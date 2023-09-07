function test_LEMON_MRI_fmri(rap)
    rap.options.parallelresources.walltime = 6;

    rap.tasksettings.reproa_fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap.tasksettings.reproa_eddy_fmri.mode = 'extensive';

    rap.tasksettings.reproa_coregextended(1).target = 'structural';
    rap = renameStream(rap,'reproa_coregextended_00001','input','append','brainmask');
    rap.tasksettings.reproa_coregextended(2).target = 'meanfmri';

    rap = renameStream(rap,'reproa_aroma_00001','input','structural','reproa_fromnifti_structural_00001.structural');
    rap = renameStream(rap,'reproa_aroma_00002','input','structural','reproa_coregextended_00002.structural');

    for b = 1:2
        rap.tasksettings.reproa_coregextended(b).reorienttotemplate = 0;

        rap.tasksettings.reproa_coregextended_t2(b).reorienttotemplate = 0;

        rap.tasksettings.reproa_segment(b).writenormalised.method = 'none';
        rap.tasksettings.reproa_segment(b).writecombined = [0.5 0.1 0 0 0 0]; % ~ BET (GM and WM only)

        rap = renameStream(rap,sprintf('reproa_mask_segmentation_%05d',b),'input','reference','meanfmri');
        rap = renameStream(rap,sprintf('reproa_mask_segmentation_%05d',b),'input','segmentations','native_segmentations');
        rap.tasksettings.reproa_mask_segmentation(b).threshold = 'exclusive';

        rap.tasksettings.reproa_smooth_fmri(b).FWHM = 4; % ~1.5 voxel size

        rap = renameStream(rap,sprintf('reproa_aroma_%05d',b),'input','structural_brain',...
                           sprintf('reproa_segment_%05d.structural',b));
        rap.tasksettings.reproa_aroma(b).highpassfilter = 100; % in seconds
        rap.tasksettings.reproa_aroma(b).includemovementparameters = [1 1 0; 1 1 0];
        rap.tasksettings.reproa_aroma(b).includesegmentationsignal = ...
            {'components-WMc1' 'components-WMc2' 'components-WMc3' 'components-WMc4' 'components-WMc5' ...
             'components-CSFc1' 'components-CSFc2' 'components-CSFc3' 'components-CSFc4' 'components-CSFc5'};
        rap.tasksettings.reproa_aroma(b).numberofdimensions = 100;
        rap.tasksettings.reproa_aroma(b).componentregression = 'both';
    end

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
