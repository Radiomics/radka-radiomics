RADIOMICS
==========

Install
==========

Installation libraries
------------------------

```
pip install javabridge
pip install python-weka-wrapper3
pip install msgpack-python
pip install pyradiomics
pip install pydicom
pip install opencv-python
pip install coloredlogs
```

PyCharm
------------------------
In order to properly configure this project in PyCharm you need to
set the `src` folder and the `imagevisualizer` and `preprocUtils` submodules as
**Sources Root**. Simply right click on the folders and
then go to `Mark directory as` --> `Sources root`. 


PROGRAMS
==========

The objective of this project is to be able to predict genomic marker from imaging. 

The data comes from T2-MRIs and DWI (ADC and BVAL). There is a total of X cases, for each of them
dwi 220 features are computed. The features being computed are the following:

  * 3 image types: t2, adc, and bval
  * 4 statistics: mean, std, kurtosis, skewness
  * 6 textures: intensity, correlation, energy, entropy
  * 5 percentiles: 10, 25, 50, 75, 90
  * 3 tissues: lesion, natural appearing peripheral zone,  natural appearing transition zone 
  * 4 extra features: prostate volume, lesion volume (ROI_vol), and ROI/Vol


In order to obtain the final 'table' with the radiomics features appended to each clinical feature from the patients
we need to run the following files: 

0_ComputeRadiomicsProstateMasks.py
-----------------------------------


1_Read_Preprocess_Data.py
------------

This program does the following processes:

1. _mergeRadiomicsByBiopsy_ Combines the radiomics table
into one row for each biopsy (Lesion, NAP and NATZ)
2. _mergeNormalition_ Merges the normalization info (Femur, Fat and Urine)
    with the clinical data 
3. _joinTablesByID_ 1) Fixes 15 Symphony (multiply by 10 and 2)
    Concatenates the clinical data with the radiomics data
4. _mainNormalizeIntensities_ Normalizes the data taking into
account (bone, fat, and urine)
