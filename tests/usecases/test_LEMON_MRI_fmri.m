function test_LEMON_MRI_fmri(rap)
    rap.tasksettings.fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
