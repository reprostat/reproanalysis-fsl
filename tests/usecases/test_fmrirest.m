function test_fmrirest(rap)
    rap.options.parallelresources.walltime = 6;

    rap.tasksettings.reproa_fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap.tasksettings.reproa_eddy_fmri.mode = 'extensive';

    rap.tasksettings.reproa_fromnifti_structural.sfxformodality = 'T1w:T2w';

    rap.tasksettings.reproa_coregextended.target = 'meanfmri';
    rap.tasksettings.reproa_coregextended.reorienttotemplate = 0;

    rap.tasksettings.reproa_coregextended_t2.reorienttotemplate = 0;

    rap.tasksettings.reproa_segment.writenormalised.method = 'none';
    rap.tasksettings.reproa_segment.writecombined = [0.5 0.1 0 0 0 0]; % ~ BET (GM and WM only)

    rap = renameStream(rap,'reproa_mask_segmentations_00001','input','reference','meanfmri');
    rap = renameStream(rap,'reproa_mask_segmentations_00001','input','segmentations','native_segmentations');
    rap.tasksettings.reproa_mask_segmentations.threshold = 'exclusive';

    rap.tasksettings.reproa_smooth_fmri.FWHM = 4; % ~1.5 voxel size

    for b = 1:2 % AROMA
        rap = renameStream(rap,sprintf('reproa_aroma_%05d',b),'input','structural',...
                           'reproa_coregextended_00001.structural');
        rap = renameStream(rap,sprintf('reproa_aroma_%05d',b),'input','structural_brain',...
                           'reproa_segment_00001.structural');
        rap.tasksettings.reproa_aroma(b).highpassfilter = 100; % in seconds
        rap.tasksettings.reproa_aroma(b).includemovementparameters = [1 1 0; 1 1 0];
        rap.tasksettings.reproa_aroma(b).numberofdimensions = 100;
        rap.tasksettings.reproa_aroma(b).componentregression = 'both';
    end
    % - AROMA + aCompCorr (WM+CSF)
    rap.tasksettings.reproa_aroma(2).includesegmentationsignal = ...
        {'components-WMc1' 'components-WMc2' 'components-WMc3' 'components-WMc4' 'components-WMc5' ...
         'components-CSFc1' 'components-CSFc2' 'components-CSFc3' 'components-CSFc4' 'components-CSFc5'};

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
