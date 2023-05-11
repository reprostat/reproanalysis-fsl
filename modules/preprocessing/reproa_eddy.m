% mask
%   fslmaths topup_fieldmap -Tmean topup_mean
%   bet topup_mean topup_brain -n -m -f 0.2
%   topup_brain_mask.nii
% index
%    dlmwrite('index.txt',ones(1,nVol),'delimiter',' ');
% bvals
%    dlmwrite('bvals',zeros(1,nVol),'delimiter',' ');
% bvecs
%    dlmwrite('bvecs',repmat([1;0;0],1,nVol),'delimiter',' ');
% slspec
%    slspec = unique(header{1}.sliceorder,'stable')-1;
%    mbf = numel(header{1}.sliceorder)/numel(slspec);
%    for b = 2:mbf
%        slspec(:,b) = slspec(:,1)+(b-1)*size(slspec,1);
%    end
%    dlmwrite('slspec.txt',slspec,'delimiter',' ');
%
% eddy_cuda10.2 \
%     --imain=$WDIR/sub-032301_ses-01_task-rest_acq-AP_run-01_bold.nii \
%     --mask=$WDIR/topup_mask_mask.nii \
%     --acqp=$WDIR/acquisition_parameters.txt \
%     --index=$WDIR/index.txt \
%     --bvecs=$WDIR/bvecs \
%     --bvals=$WDIR/bvals \
%     --topup=$WDIR/topup \
%     --out=$WDIR/eddy-cuda-full \
%     --niter=10 --fwhm=10,10,5,5,0,0,0,0,0,0 --repol \
%     --slspec=$WDIR/slspec.txt --mporder=8 --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear \
%     --estimate_move_by_susceptibility --mbs_niter=20 --mbs_lambda=5 --mbs_ksp=5 \
%     --b0_only --data_is_shelled --verbose 
function rap = reproa_eddy(rap,command,subj,run)

switch command
    case 'doit'
        
        % Get nii filenames from stream
        pelabels=[]; % for each volume, which acquisition does it come from
        bvals=[];
        bvecs=[];
        allfns='';
        for peind=1:2
            diffinput=aas_getfiles_bystream(rap,'diffusion_session_phaseencode_direction',[subj run peind],'diffusion_data');
            allfns=[allfns ' ' diffinput];
            V=spm_vol(diffinput);
            pelabels=[pelabels;peind*ones(length(V),1)];
            bvals=[bvals load(aas_getfiles_bystream(rap,'diffusion_session_phaseencode_direction',[subj run peind],'bvals'))];
            bvecs=[bvecs load(aas_getfiles_bystream(rap,'diffusion_session_phaseencode_direction',[subj run peind],'bvecs'))];
        end;
        % Output directory
        dsess=aas_getpath_bydomain(rap,'diffusion_session',[subj run]);
         
        % Collect other files eddy needs
        BETmask=aas_getfiles_bystream(rap,'diffusion_session',[subj run],'BETmask');
        BETmask=BETmask(1,:); % I hate having to do this filtering, one day should refactor BET outputs
        
        acqparms=aas_getfiles_bystream(rap,'diffusion_session',[subj run],'topup_acquisition_parameters');
        topup_output_fieldcoef=aas_getfiles_bystream(rap,'diffusion_session',[subj run],'topup_output_fieldcoef');        
        [pth nme ext]=fileparts(topup_output_fieldcoef);
        topup_output_rootfn=fullfile(pth,nme(1:end-10));
        
        % Write phase encode labels (rows in acq parms)
        pelabelsfn=fullfile(dsess,'pelabels.txt');  
        fid=fopen(pelabelsfn,'w');
        fprintf(fid,'%d\n',pelabels);
        fclose(fid);
        
        % Write bvals
        bvalsfn=fullfile(dsess,'bvals');
        fid=fopen(bvalsfn,'w');
        fprintf(fid,'%f ',bvals);
        fprintf(fid,'\n');
        fclose(fid);
        
        % Write bvecs
        bvecsfn=fullfile(dsess,'bvecs');
        fid=fopen(bvecsfn,'w');
        fprintf(fid,'%f ',bvecs(1,:));
        fprintf(fid,'\n');
        fprintf(fid,'%f ',bvecs(2,:));
        fprintf(fid,'\n');
        fprintf(fid,'%f ',bvecs(3,:));
        fprintf(fid,'\n');
        fclose(fid);
        
        
        % Create merged file
        mergedfn=fullfile(dsess,'allpe.nii');
        cmd=['fslmerge -t ' mergedfn ' ' allfns];
        [s w]=aas_runfslcommand(rap,cmd);
        
        % Now run eddy
        outfn=fullfile(dsess,'eddy_output');       
        cmd=[sprintf('eddy --imain=%s --mask=%s --acqp=%s --index=%s ',mergedfn,BETmask,acqparms,pelabelsfn)...
             sprintf('--bvecs=%s --bvals=%s --topup=%s ',bvecsfn,bvalsfn,topup_output_rootfn)...
             sprintf('--out=%s',outfn)];
         
        aas_log(rap,false,sprintf('Running %s',cmd));        
         [s w]=aas_runfslcommand(rap,cmd);
        if s
             aas_log(rap,true,sprintf('Error %s',w));
        end;
        
        output_list={'eddy_movement_rms','eddy_outlier_map','eddy_outlier_n_stdev_map','eddy_outlier_report','eddy_parameters','eddy_post_eddy_shell_alignment_parameters','eddy_rotated_bvecs','nii'};
        
        % Describe outputs
        for outputind=1:length(output_list)
            rap=aas_desc_outputs(rap,'diffusion_session',[subj run],'diffusion_data',[outfn '.' output_list{outputind}]);
        end;
       
        rap=aas_desc_outputs(rap,'diffusion_session',[subj run],'bvals',bvalsfn);
        rap=aas_desc_outputs(rap,'diffusion_session',[subj run],'bvecs',bvecsfn);
        
end
end

