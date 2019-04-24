#!/bin/bash

# ------------------------------- Patient 000 -----------------------------
# T2
SOURCE=/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/PX0000/ProstateX-0000_ProstateX-0000_MR_2011-07-07_114731_MR.prostaat.kanker.detectie.WDS.mc.MCAPRODETW_t2.tse.tra_n19__00000
MASK=/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/PX0000/ProstateX-0000_ProstateX-0000_MR_2011-07-07_114731_MR.prostaat.kanker.detectie.WDS.mc.MCAPRODETW_T2.Mask_n19__00000
OUTPUT=/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/OUTPUT/T2
./build/BatchTextureAnalyzerExec $SOURCE $MASK $OUTPUT  2 16
