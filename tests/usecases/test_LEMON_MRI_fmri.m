function test_LEMON_MRI_fmri(rap)
    rap.tasksettings.reproa_fromnifti_fieldmap.pattern = 'acq-SEfmapBOLDpre';

    rap.tasksettings.reproa_eddy_fmri.mode = 'extensive';

    rap.tasksettings.reproa_coregextended(1).target = 'structural';
    rap.tasksettings.reproa_coregextended(1).reorienttotemplate = 0;
    rap.tasksettings.reproa_coregextended(2).target = 'meanfmri';
    rap.tasksettings.reproa_coregextended(2).reorienttotemplate = 0;

%    rap.tasksettings.reproa_segment.writenormalised.method = 'none';

    rap = processBIDS(rap);

    processWorkflow(rap);

    reportWorkflow(rap);
end
