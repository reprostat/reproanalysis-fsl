# reproanalysis-fsl
Reproanalysis Extension for integrating FSL

## Requirements
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
- [ICA-AROMA](https://github.com/maartenmennes/ICA-AROMA) - required for the _reproa_aroma_ module

## Supported use-cases
- fmriprep-like preprocessing ([task list](https://github.com/reprostat/reproanalysis-fsl/blob/main/examples/FSL_rest.xml) and [user script](https://github.com/reprostat/reproanalysis-fsl/blob/main/tests/usecases/test_fmrirest.m)), including [_topup_](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup), [_eddy_](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy), [_ICA-AROMA_](https://github.com/maartenmennes/ICA-AROMA) with or without _aCompCorr_
