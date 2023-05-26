function test_LEMON_MRI_fmri(rap)
    rap.options.parallelresources.walltime = 6;

    rap.tasksettings.reproa_fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap.tasksettings.reproa_eddy_fmri.mode = 'extensive';

    rap.tasksettings.reproa_coregextended(1).target = 'structural';
    rap.tasksettings.reproa_coregextended(2).target = 'meanfmri';

    for b = 1:2
        rap.tasksettings.reproa_coregextended(b).reorienttotemplate = 0;

        rap.tasksettings.reproa_coregextended_t2(b).reorienttotemplate = 0;

        rap.tasksettings.reproa_segment(b).writenormalised.method = 'none';

        rap = renameStream(rap,sprintf('reproa_mask_segmentation_%05d',b),'input','reference','meanfmri');
        rap = renameStream(rap,sprintf('reproa_mask_segmentation_%05d',b),'input','segmentations','native_segmentations');
        rap.tasksettings.reproa_mask_segmentation(b).threshold = 'exclusive';
    end

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
