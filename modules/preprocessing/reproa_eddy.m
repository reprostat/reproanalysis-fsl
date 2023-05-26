function rap = reproa_eddy(rap,command,subj,run)

switch command
    case 'doit'
        sliceAxes = 'ijk';
        pfxEddy = 'eddy';

        % Detect CUDA
        [s, w] = shell('nvidia-smi -q | grep CUDA','ignoreerror',true);
        verCUDA = ~s && ~isempty(w);
        if verCUDA, verCUDA = str2double(char(regexp(w,'(?<=:)[0-9\ \.]*','match'))); end

        % Output directory
        localFolder = getPathByDomain(rap,rap.tasklist.currenttask.domain,[subj run]);

        % Data
        fnData = getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri','streamType','input');

        % Mask
        fnFM = char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'dualpefieldmap'));
        runFslCommand(rap,['fslmaths ' fnFM ' -Tmean ' strrep(fnFM,'fieldmap.nii','mean.nii')]);
        runFslCommand(rap,['bet ' strrep(fnFM,'fieldmap.nii','mean.nii') ' ' strrep(fnFM,'fieldmap.nii','brain.nii') ' -n -m -f 0.2']);
        fnMask = strrep(fnFM,'fieldmap.nii','brain_mask.nii');

        % Acqisition parameters
        fnAcq = getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'acquisitionparameters');
        tableTopup = dlmread(fnAcq{1});

        % Index
        fnAcqInd = fullfile(localFolder,'index.txt');
        nVol = numel(spm_vol(fnData{1}));
        load(char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri_header')),'header');
        indTable = cellfun(@(p) ...
                           min(find(arrayfun(@(r) ...
                                             isequal((sliceAxes == p(1))*str2double(fliplr(regexprep(p,['[' sliceAxes ']'],'1'))),tableTopup(r,1:3)),...
                                             1:size(tableTopup,1)))),...
                           cellstr(header{1}.PhaseEncodingDirection));
        dlmwrite(fnAcqInd,indTable*ones(1,nVol),'delimiter',' ');

        % Topup
        fnTopup = regexprep(char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fieldcoefficients')),'_fieldcoef.nii.*','');

        % b
        if hasStream(rap,rap.tasklist.currenttask.domain,[subj run],'bvals')
            fnBVals = char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'bvals'));
        else
            fnBVals = fullfile(localFolder,'bvals');
            dlmwrite(fnBVals,zeros(1,nVol),'delimiter',' ');
        end
        if hasStream(rap,rap.tasklist.currenttask.domain,[subj run],'bvecs')
            fnBVecs = char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'bvecs'));
        else
            fnBVecs = fullfile(localFolder,'bvecs');
            dlmwrite(fnBVecs,repmat([1;0;0],1,nVol),'delimiter',' ');
        end

        % eddy parameters
        strQCParam = '';
        % - binary
        switch getSetting(rap,'wheretoprocess')
            case 'auto'
                strEddy = 'eddy';
            case 'cpu'
                strEddy = 'eddy_cpu';
            case 'gpu'
                binaries = cellstr(spm_select('List',fullfile(rap.directoryconventions.fsldir,'bin'),'^eddy_cuda.*'));
                switch numel(binaries)
                    case 0
                        logging.error('No CUDA version of Eddy is detected in %s',fullfile(rap.directoryconventions.fsldir,'bin'));
                    case 1
                        strEddy = binaries{1};
                    otherwise
                        verEddy = sort(cellfun(@(b) sscanf(b,'eddy_cuda%f'), binaries));
                        strEddy = sprintf('eddy_cuda%1.1f', verEddy(find(verEddy <= verCUDA,1,'last')));
                end
        end
        % - mode
        switch getSetting(rap,'mode')
            case 'minimum'
                strEddyParam = '-fwhm=0';
            case 'medium'
                strEddyParam = '--niter=7 --fwhm=10,0,0,0,0,0,0';
            case 'extensive'
                strEddyParam = '--niter=10 --fwhm=10,10,5,5,0,0,0,0,0,0';
        end
        % - movements
        if getSetting(rap,'susceptibilitybymovement'), strEddyParam = [strEddyParam ' --estimate_move_by_susceptibility --mbs_niter=20 --mbs_lambda=5 --mbs_ksp=5']; end
        if getSetting(rap,'slicetovolume')
            doStV = true;
            if strcmp(strEddy,'eddy_cpu')
                logging.warning('slicetovolume correction is availabe only in the CUDA version');
                doStV = false;
            end
            if doStV && ~verCUDA
                logging.warning('no CUDA is detected');
                doStV = false;
            end
            if doStV && ~isfield(header{1},'sliceorder')
                logging.warning('sliceorder is not detected');
                doStV = false;
            end
            if ~doStV, logging.warning('\t-> slicetovolume correction will be skipped');
            else
                fnSlSpec = fullfile(localFolder,'slspec.txt');
                slspec = unique(header{1}.sliceorder,'stable')-1;
                mbf = numel(header{1}.sliceorder)/numel(slspec);
                for b = 2:mbf
                    slspec(:,b) = slspec(:,1)+(b-1)*size(slspec,1);
                end
                dlmwrite(fnSlSpec,slspec,'delimiter',' ');
                strEddyParam = [strEddyParam ' --slspec=' fnSlSpec ' --mporder=8 --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear'];
                strQCParam = [strQCParam ' --slspec=' fnSlSpec];
            end
        end
        % - data
        if strcmp(rap.tasklist.currenttask.domain,'fmrirun'), strEddyParam = [strEddyParam ' --b0_only']; end

        % Now run eddy
        runFslCommand(rap,sprintf('%s --imain=%s --mask=%s --acqp=%s --index=%s --bvecs=%s --bvals=%s --topup=%s --out=%s %s --repol --data_is_shelled --verbose', ...
                                  strEddy,fnData{1},fnMask,fnAcq{1},fnAcqInd,fnBVecs,fnBVals,fnTopup,fullfile(localFolder,pfxEddy),strEddyParam));
        fnFmri = spm_select('FPList',localFolder, ['^' pfxEddy '\.[^(eddy)]']);

        % Generate mean fmri
        fnMeanFmri = spm_file(fnFmri(1,:),'prefix','mean_');
        runFslCommand(rap,sprintf('fslmaths %s -Tmean %s',fnFmri,fnMeanFmri));

        % Outputs
        putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri',fnFmri);
        putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'meanfmri',fnMeanFmri);
        putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'brainmask',fnMask);
        putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'movementparameters',fullfile(localFolder,[pfxEddy '.eddy_parameters']));

        % QA (for DWI only)
        if hasStream(rap,rap.tasklist.currenttask.domain,[subj run],'bvals')
            fnFM = getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'dualpefieldmap');
            runFslCommand(rap,sprintf('eddy_quad %s --mask=%s --eddyParams=%s --eddyIdx=%s --bvals=%s --field=%s %s',...
                                      fullfile(localFolder,pfxEddy),fnMask,fnAcq{1},fnAcqInd,fnBVals,fnFM{1},strQCParam(2:end)));
        end
end
end

