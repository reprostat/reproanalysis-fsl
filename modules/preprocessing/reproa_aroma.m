function rap = reproa_aroma(rap, command, subj, run)
% denoise an EPI using FSL's ICA-AROMA
%
% This requires AROMA-ICA to be installed and added as an aa toolbox.
% AROMA requires Python (preferably 3.6) and pip. reproa also requires it to be installed in a conda environment, as specified in the parameterset.
%
% The process to install is
%   $ cd $HOME/tools
%   $ git clone https://github.com/tiborauer/ICA-AROMA.git
%   $ conda create -n aroma python=3.6
%   $ conda activate aroma
%   $ pip install -r ICA-AROMA/requirements.txt
%
%   Add the corresponding entry to your parameterset
%     <toolbox desc="Toolbox with implemented interface in extrafunctions/toolboxes" ui="custom">
%         <name desc="Name corresponding to the name of the interface without the 'Class' suffix" ui="text">aroma</name>
%         <dir ui="dir">$HOME/tools/ICA-AROMA</dir>
%         <extraparameters>
%             <condaEnvironment>aroma</condaEnvironment>
%         </extraparameters>
%     </toolbox>
%
% Note this is unsupervised denoising, which as a rule does not perform as well as using ICA with training data, such as FSL FIX.
%
% Based on https://github.com/automaticanalysis/automaticanalysis/blob/master/aa_modules/aamod_AROMA_denoise.m by Michael S. Jones (@jones-michael-s)

switch command

    case 'report'

    case 'doit'
        %% Init
        localPath = getPathByDomain(rap,'fmrirun', [subj run]);
        aromaPath = fullfile(localPath,'AROMA');
        if exist(aromaPath,'dir'), dirRemove(aromaPath); end

        global reproacache
        AROMA = reproacache('toolbox.aroma');
        fnAROMA = fullfile(AROMA.toolPath,'ICA_AROMA.py');

        % Highpass filter (for data and noise)
        fnFmri = char(getFileByStream(rap, 'fmrirun',[subj run], 'fmri'));
        V = spm_vol(fnFmri);
        load(char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri_header')),'header');
        K = spm_filter(struct('RT',header{1}.volumeTR,...
                              'row',1:numel(V),...
                              'HParam',getSetting(rap,'highpassfilter')));

        %% Process noise parameters
        % motion correction parameters, we expect 6 DOF as text file
        fnMoCo = char(getFileByStream(rap, 'fmrirun',[subj run], 'movementparameters'));
        movMat = getSetting(rap,'includemovementparameters');
        rp = load(fnMoCo); %
        movRegFinal = zeros(size(rp,1),0);
        for o = 1:size(movMat,1)
            for d = 1:size(movMat,2)
                if movMat(o,d)
                    movReg = rp.^o;
                    for i = 2:d
                        if d < size(movMat,2) % gradient/derivative
                            movReg = gradient(movReg);
                        else % spin history
                            movReg = [zeros(1,size(movReg,2)); diff(movReg)];
                        end
                    end
                    movRegFinal = [movRegFinal movReg];
                end
            end
        end

        % segmentation signals, we expect tsv according to BIDS
        if ~isempty(getSetting(rap,'includesegmentationsignal'))
            streamSignal = getFileByStream(rap,'fmrirun',[subj run],'segmentationsignal');
            for sig = getSetting(rap,'includesegmentationsignal')
                sigSpec = strsplit(sig{1},'-');
                MSig = importdata(streamSignal.(sigSpec{1}){1},'\t',1);
                movRegFinal = [movRegFinal MSig.data(:,strcmp(MSig.colheaders,sigSpec{2}))];
            end
        end

        % filter and save
        fnNoise = fullfile(localPath,'noise');
        dlmwrite(fnNoise,spm_filter(K,movRegFinal),'  ','precision',12);

        %% Process structural
        fnStructural = char(getFileByStream(rap, 'subject',subj, 'structural'));

        % this is the MNI structural reference
        % using the 1mm reference takes >20x longer to run than using 2mm
        fnMNI = fullfile(rap.directoryconventions.fsldir,['data/standard/MNI152_T1_' getSetting(rap,'normalisationresolution')]);

        % brain extraction
        if hasStream(rap, 'subject',subj, 'structural_brain')
            fnStructuralBrain = char(getFileByStream(rap, 'subject',subj, 'structural_brain'));
        else
            fnStructuralBrain = fullfile(localPath,'struct_brain');
            logging.info('Brain extraction on %s...', fnStructural);
            runFslCommand(rap, sprintf('bet %s %s', fnStructural, fnStructuralBrain));
        end

        % pre-register (FLIRT) structural to MNI
        fnMStruct2MNI = fullfile(localPath,'struct2MNI.mat');
        logging.info('Running FLIRT on %s...', fnStructuralBrain);
        runFslCommand(rap, sprintf('flirt -in %s -ref %s -omat %s -dof 12', fnStructuralBrain, spm_file(fnMNI,'suffix','_brain'), fnMStruct2MNI));

        % register (FNIRT) structural to MNI
        fnWStruct2MNI = fullfile(localPath,'struct2MNI.nii');
        logging.info('Running FNIRT on %s...', fnStructural);
        runFslCommand(rap,sprintf('fnirt --in=%s --ref=%s --refmask=%s --aff=%s --cout=%s --config=T1_2_MNI152_2mm',...
                                  fnStructural, fnMNI, spm_file(fnMNI,'suffix','_brain_mask_dil'), fnMStruct2MNI, fnWStruct2MNI));

        %% Process fMRI
        % highpass filter (and re-add mean)
        fnFmri = spm_file(fnFmri,'prefix','hp_');
        [R,C,P] = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
        Y = reshape(spm_filter(K,...
                               spm_get_data(V,...
                                            [R(:)';C(:)';P(:)']))',...
                    V(1).dim(1),V(1).dim(2),V(1).dim(3),[]);

        % - write 4D NIfTI
        dim = size(Y);
        N      = nifti;
        N.dat  = file_array(fnFmri,dim,V(1).dt);
        N.mat  = V(1).private.mat;
        N.mat0 = V(1).private.mat;
        N.descrip     = 'filtered';
        create(N);

        dim = [dim 1];
        for i = 1:prod(dim(4:end))
            N.dat(:,:,:,i) = Y(:,:,:,i);
            spm_get_space([N.dat.fname ',' num2str(i)], V(1).mat);
        end
        N.dat = reshape(N.dat,dim);

        % create meanfMRI for registration
        fnMeanFmri = spm_file(fnFmri,'prefix','mean_');
        runFslCommand(rap, sprintf('fslmaths %s -Tmean %s', fnFmri, fnMeanFmri));

        % rigid-body (6DOF) register (FLIRT) functional to structural
        fnMFunc2Struct = fullfile(localPath,'func2struct.mat');
        logging.info('Running FLIRT on %s...', fnStructuralBrain);
        runFslCommand(rap, sprintf('flirt -in %s -ref %s -omat %s -dof 6', fnMeanFmri, fnStructuralBrain, fnMFunc2Struct));

        %% Run AROMA
        switch getSetting(rap,'componentregression')
            case 'partial'
                parDenoise = 'nonaggr';
                pttrOut = '^[w]?denoised_func_data_nonaggr..*';
            case 'full'
                parDenoise = 'aggr';
                pttrOut = '^[w]?denoised_func_data_aggr..*';
            case 'both'
                parDenoise = 'both';
                pttrOut = '^[w]?denoised_func_data_.*aggr..*';
        end
        fsloutputtype0 = rap.directoryconventions.fsloutputtype;
        rap.directoryconventions.fsloutputtype = 'NIFTI2_GZ';
        logging.info('Running AROMA on %s...', fnFmri);
        runPyCommand(rap,...
                     sprintf('python %s -i %s -o %s -tr %1.3f -mc %s -a %s -w %s -m %s -dim %d -den %s',...
                             fnAROMA, fnFmri, aromaPath, header{1}.volumeTR, fnNoise, fnMFunc2Struct, fnWStruct2MNI,...
                             char(getFileByStream(rap,'fmrirun',[subj run],'brainmask')),getSetting(rap,'numberofdimensions'),parDenoise),...
                     AROMA.condaEnvironment,'runFsl',true);
        rap.directoryconventions.fsloutputtype = fsloutputtype0;

        %% Output
        % output will go in <aromaPath>/denoised_func_data_nonaggr.nii[.gz]
        for f = spm_file(spm_file(cellstr(spm_select('FPList',aromaPath,pttrOut)),'ext',''),'ext','')'
            if getSetting(rap,'writenormalised')
                logging.info('Applying FNIRT on %s...', f{1});
                runFslCommand(rap,sprintf('applywarp --in=%s --ref=%s --warp=%s --out=%s',...
                                          f{1}, fnMNI, fnWStruct2MNI, spm_file(f{1},'path',localPath,'prefix','w')));
            else
                if ~endsWith(rap.directoryconventions.fsloutputtype,'GZ'), gunzip([f{1} '.nii.gz'],localPath); delete([f{1} '.nii.gz']);
                else, movefile(f{1},localPath);
                end
            end
        end
        out = cellstr(spm_select('FPList',localPath,pttrOut));
        if numel(out) > 1
            out = reshape([{'full' 'partial'}; out'],1,[]);
            out = struct(out{:});
        end
        putFileByStream(rap,'fmrirun',[subj run],'fmri',out);

        %% Diagnostics
        visFig = 'off';

        % read classification
        fnClass = fullfile(aromaPath,'classified_motion_ICs.txt');
        iCBad = dlmread(fnClass);

        % get structural
        bg = char(getFileByStream(rap,'subject',subj,'structural_brain'));

        offX = 0; offY = 0;
        for iC = 1:size(spm_select('List',fullfile(aromaPath,'melodic.ica','report'),'^t[0-9]*\.txt'),1)
            logging.info('Processing IC: %d\n',iC);

            % Spatial component
            fnStat = sprintf(fullfile(aromaPath,'melodic.ica','stats','thresh_zstat%d.nii.gz'),iC);
            Y = spm_read_vols(spm_vol(fnStat));
            fig = spm_figure('GetWin','Graphics');
            set(fig,'visible',visFig);
            spm_check_registration(bg);
            spm_orthviews('addtruecolourimage',1,fnStat,[0.5 0.5 0.5; hot],0.6,max(Y(:)),min(Y(Y>0)));
            spm_orthviews('Reposition',[0 0 20]); % assuming neck -> brain is a bit higher
            global st
            slicesSummary = cell(1,3);
            spm_orthviews('Xhairs','off');
            for a = 1:3
                fr = getframe(st.vols{1}.ax{a}.ax);
                slicesSummary{a} = fr.cdata;
            end
            close(fig);

            img = slicesSummary{3};
            img(1:size(slicesSummary{2},1),end+1:end+size(slicesSummary{2},2),:) = slicesSummary{2};
            img(1:size(slicesSummary{1},1),end+1:end+size(slicesSummary{1},2),:) = slicesSummary{1};
            img = img(11:end-10,11:end-10,:);

            % - remove background
            [nVal, val] = hist(img(:),0:255);
            [~,oVal] = sort(nVal);
            for cBg = val(oVal(end-1:end))
                img(repmat((img(:,:,1) == cBg) & (img(:,:,2) == cBg) & (img(:,:,3) == cBg),1,1,3)) = 255;
            end

            % temporal component
            fnTs = sprintf(fullfile(aromaPath,'melodic.ica','report','t%d.txt'),iC);
            fig = figure('visible',visFig);
            set(fig, 'Position',[1 1 size(img,1)*4 size(img,1)])
            ts = dlmread(fnTs);
            if ismember(iC,iCBad)
                plot(1:numel(ts),ts,'r','LineWidth',1.5);
            else
                plot(1:numel(ts),ts,'g','LineWidth',1.5);
            end
            set(gca,'xaxislocation','origin','Units','pixels','Position',[1 1 size(img,1)*4 size(img,1)])
            legend(sprintf('IC: %d',iC),'location','northwest')
            set(gca,'FontSize',24);
            fr = getframe(gca);
            close(fig);

            img = [img fr.cdata];

            imgAll(offY+1:offY+size(img,1),offX+1:offX+size(img,2),:) = img;
            if mod(iC,2) % odd
                offX = 0;
                offY = offY + size(img,1);
            else % even
                offX = offX + size(img,2);
            end
        end

        imwrite(imgAll,fullfile(localPath,sprintf('diagnostic_%s.jpg',rap.tasklist.currenttask.name)));

        %% Cleanup
        if getSetting(rap,'cleanup')
            delete(fnFmri);
            delete(fnMeanFmri);
            dirRemove(fullfile(aromaPath,'melodic.ica'));
        end

        % aside: if you delete DENOISE, it would remove about about 99%
        % of unnecessary files (which might be more than 1 GB!). Maybe make
        % this an option. Could also look into garbage collection.
        % N.B.: DENOISE/ICA_AROMA_component_assessment.pdf might be worth
        % keeping

    case 'checkrequirements'
        % Check environment
        % - FSL
        if isempty(rap.directoryconventions.fsldir) ...
                || ~exist(rap.directoryconventions.fsldir,'dir')
           logging.error('%s:rap.directoryconventions.fsldir is not set or do not exists.',mfilename);
        end
        runFslCommand(rap,'flirt -version');
        [s, w]= runFslCommand(rap,'which imcp',{},'quiet',true);
        if ~strcmp(deblank(w),fullfile(rap.directoryconventions.fsldir,'share/fsl/bin/imcp'))
            logging.error('FSL''s python environment is not set up correctly.');
        end

        % - AROMA
        global reproacache
        if ~reproacache.isKey('toolbox.aroma'), logging.error('Toolbox ICA-AROMA is not defined.'); end
        AROMA = reproacache('toolbox.aroma'); fnAROMA = fullfile(AROMA.toolPath,'ICA_AROMA.py');
        if ~exist(fnAROMA,'file'), logging.error('%s: ICA-AROMA executable %s not found.', mfilename, fnAROMA);end

        % - conda environment
        [s,w] = runPyCommand(rap,'conda info -e','','quiet',true);
        if ~any(strcmp(regexp(w,'(?<=\n)[^\ ]*','match'),AROMA.condaEnvironment))
            logging.error('%s: conda environment %s not found.', mfilename, AROMA.condaEnvironment);
        end
end

end
