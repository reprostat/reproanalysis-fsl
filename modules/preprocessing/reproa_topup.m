function rap = reproa_topup(rap,command,subj,run)

    switch command

        case 'doit'
            sliceAxes = 'ijk';
            pfxTopup = 'topup';

            % Output directory
            localFolder = getPathByDomain(rap,rap.tasklist.currenttask.domain,[subj run]);

            % Get inputs
            fnAll = getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'dualpefieldmap','streamType','input');
            load(char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'dualpefieldmap_header','streamType','input')),'dcmhdr');

            % Table
            % - PE
            tableTopup = cellfun(@(p) (sliceAxes == p(1))*str2double(fliplr(regexprep(p,'[ijk]','1'))), dcmhdr{1}.PhaseEncodingDirection, 'UniformOutput',false);
            tableTopup = cat(1,tableTopup{:});
            % - Total Readout Time for FSL
            if isfield(dcmhdr{1},'TotalReadoutTime') % directly from the header
                tableTopup(:,4) = dcmhdr{1}.TotalReadoutTime;
            elseif isfield(dcmhdr{1},'EffectiveEchoSpacing') % based on Echo Spacing
                headerFields = fieldnames(dcmhdr{1});
                if any(strcmpi(headerFields,'NumberOfPhaseEncodingSteps'))
                    nPE = dcmhdr{1}.(headerFields{strcmpi(headerFields,'NumberOfPhaseEncodingSteps')});
                else
                    V = spm_vol(fnAll{1});
                    nPE = V(1).dim(logical(abs(tableTopup(1,:))));
                end
            else
                logging.error('TotalReadoutTime or EffectiveEchoSpacing MUST be specified in the fieldmap header');
            end
            % - match table rows with volumes
            tableTopup = [repmat(tableTopup(1,:),numel(spm_vol(fnAll{1})),1);...
                          repmat(tableTopup(2,:),numel(spm_vol(fnAll{2})),1)];
            % - save
            fnAcq = fullfile(localFolder,'acquisition_parameters.txt');
            fid = fopen(fnAcq,'w');
            for ind = 1:size(tableTopup,1), fprintf(fid,'%d %d %d %1.6f\n',tableTopup(ind,:)); end
            fclose(fid);
            putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'acquisitionparameters',fnAcq);

            % Create merged file
            fnAllPE = fullfile(localFolder,'allpe.nii');
            runFslCommand(rap,['fslmerge -t ' fnAllPE ' ' strjoin(fnAll,' ')]);

            % Estimate topup
            runFslCommand(rap,sprintf('topup --imain=%s --datain=%s --config=%s/etc/flirtsch/b02b0.cnf  --out=%s --iout=%s',...
                                      fnAllPE, fnAcq, rap.directoryconventions.fsldir, fullfile(localFolder,pfxTopup), fullfile(localFolder,[pfxTopup '_fieldmap'])));

            putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fieldcoefficients',spm_select('FPList',localFolder, ['^' pfxTopup '_fieldcoef.*']));
            putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'movementparameters',fullfile(localFolder, [pfxTopup '_movpar.txt']));
            putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'dualpefieldmap',spm_select('FPList',localFolder, ['^' pfxTopup '_fieldmap.*']));
 

            % Apply topup
            if getSetting(rap,'applytopup')
                % - data
                fnData = getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri','streamType','input'); fnData = strjoin(fnData,',');
                load(char(getFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri_header','streamType','input')),'header');

                % - index
                indTable = cellfun(@(p) ...
                                   min(find(arrayfun(@(r) ...
                                                     isequal((sliceAxes == p(1))*str2double(fliplr(regexprep(p,'[ijk]','1'))),tableTopup(r,1:3)),...
                                                     1:size(tableTopup,1)))),...
                                   cellstr(header{1}.PhaseEncodingDirection));
                
                % - run
                fnOut = spm_file(fnData,'prefix',[pfxTopup '_']);
                try
                    runFslCommand(rap,sprintf('applytopup --imain=%s  --datain=%s  --inindex=%s --topup=%s --out=%s --datatype=preserve',...
                                              fnData,fnAcq,strrep(num2str(indTable),'  ',','),fullfile(localFolder,pfxTopup),fnOut));
                catch % (default) least-squares restoration method often fails -> try jacobian method
                    runFslCommand(rap,sprintf('applytopup --imain=%s  --datain=%s  --inindex=%s --topup=%s --method=jac --out=%s --datatype=preserve',...
                                              fnData,fnAcq,strrep(num2str(indTable),'  ',','),fullfile(localFolder,pfxTopup),fnOut));
                end
                putFileByStream(rap,rap.tasklist.currenttask.domain,[subj run],'fmri',fnOut);
            end           
            
        case 'checkrequirements'
            if ~getSetting(rap,'applytopup') % Remove (optional) input and output stream not to be used and created
                for input = setdiff({rap.tasklist.currenttask.inputstreams.name},{'dualpefieldmap' 'dualpefieldmap_header'})
                    rap = renameStream(rap,rap.tasklist.currenttask.name,'input',input{1},[]);
                    logging.info('REMOVED: %s input stream: %s', rap.tasklist.currenttask.name,input{1});
                    if contains(input{1},'header'), continue; end
                    rap = renameStream(rap,rap.tasklist.currenttask.name,'output',input{1},[]);
                    logging.info('REMOVED: %s output stream: %s', rap.tasklist.currenttask.name,input{1});
                end
            end
    end
end

