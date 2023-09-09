function rap = reproa_bet(rap,command,varargin)

switch command
    case 'report'
        registrationReport(rap,varargin{:});

    case 'doit'
        % Get inputs
        streamName = rap.tasklist.currenttask.inputstreams(1).name;
        imgInput = char(getFileByStream(rap,rap.tasklist.currenttask.domain,[varargin{:}],streamName));

        imgOutput = spm_file(imgInput,'prefix','bet_');

        fslcommand = sprintf('bet %s %s -f %f -v',imgInput,imgOutput, getSetting(rap,'f'));

        doExtra = '';
        % N.B.1: mutually exclusive this is the order of preference
        % N.B.2: ensure that the name of the option start with "do[the parameter]" (e.g. dofmri -> F)
        for doOption = {'dofmri' 'dorobust' 'dobiasandneck'}
            if ~isempty(getSetting(rap,doOption{1})) && getSetting(rap,doOption{1})
                doExtra = doOption{1};
                break;
            end
        end
        if ~isempty(doExtra)
            [~,~,attr] = getSetting(rap,doExtra);
            logging.info(attr.desc);
            fslcommand = [fslcommand ' -' upper(doExtra(3))];
            switch doExtra
                case 'dofmri'
                    rap.tasklist.currenttask.settings.f = 0.3; % for the 2nd pass
            end
        end

        [~, w] = runFslCommand(rap, fslcommand);

        % This outputs last radius from recursive command...
        strRad = regexp(w,'(?<=radius\ )[0-9\.]*','match'); strRad = strRad{end};

        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(imgOutput));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];

        logging.info('\t...calculated COG (vox): %0.4f %0.4f %0.4f and radius (mm): %s', ...
                     COG(1), COG(2), COG(3), strRad)

        logging.info('2nd BET pass extracting brain masks')
        % Run BET [-A Now let's get the brain masks and meshes!!]
        [~, w] = runFslCommand(rap, ...
                               sprintf('bet %s %s -f %f -c %0.4f %0.4f %0.4f -r %s -m -v -A',imgInput,imgOutput, ...
                                       getSetting(rap,'f'), COG(1), COG(2), COG(3), strRad)...
                               );

        % This outputs last radius from recursive command...
        strRad = regexp(w,'(?<=radius\ )[0-9\.]*','match'); strRad = strRad{end};

        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(imgOutput));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];

        logging.info('\t...final COG (vox): %0.4f %0.4f %0.4f and radius (mm): %s', ...
                     COG(1), COG(2), COG(3), strRad)

        %% DESCRIBE OUTPUTS!
        outstreams = {rap.tasklist.currenttask.outputstreams.name};
        putFileByStream(rap,rap.tasklist.currenttask.domain,[varargin{:}],...
                        outstreams{~endsWith(outstreams,'mask') & ~endsWith(outstreams,'mesh')},imgOutput)
        putFileByStream(rap,rap.tasklist.currenttask.domain,[varargin{:}],...
                        outstreams{endsWith(outstreams,'mask')},spm_file(imgOutput,'suffix','_mask'))
        putFileByStream(rap,rap.tasklist.currenttask.domain,[varargin{:}],...
                        outstreams{endsWith(outstreams,'mesh')},spm_file(imgOutput,'suffix','_mesh','ext','.vtk'))


        %% DIAGNOSTIC IMAGE
        fnMesh = cellstr(spm_select('FPListRec',getPathByDomain(rap,rap.tasklist.currenttask.domain,[varargin{:}]),'.*mesh\.nii'));
        registrationCheck(rap,rap.tasklist.currenttask.domain,[varargin{:}],imgInput,imgOutput,fnMesh{:},'mode','combined')
end
end

