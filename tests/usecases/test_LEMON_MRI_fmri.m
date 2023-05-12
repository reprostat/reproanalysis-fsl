function test_LEMON_MRI_fmri(rap)
    rap.tasksettings.reproa_fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap.tasksettings.reproa_eddy.mode = 'extensive';

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
