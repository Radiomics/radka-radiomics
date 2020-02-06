from os.path import join
import os
from preprocRadiomics.preprocRadiomics import *
from preprocRadiomics.preprocMain import *
from preprocRadiomics.nomalizeData import *
from Normalization_Tests import mainNormalizeIntensities
from config.MainConfig import get_config_preproc_data
from constants.ConfigParams import ParamsPreprocData

__author__='Olmo S. Zavala Romero'

if __name__ == '__main__':
    # This code does the following:
    # 1: Merges the Radiomics file putting the NAPZ, NATZ and GM radiomics into a single row

    config = get_config_preproc_data()

    # ------------------- Merging radiomics -------------------
    input_folder = config[ParamsPreprocData.input_folder]
    output_folder = config[ParamsPreprocData.output_folder]
    input_radiomics = config[ParamsPreprocData.input_radiomics_file]
    input_biopsy_clinical_info = config[ParamsPreprocData.input_clinical_biopsy_info_file]
    input_mean_intensities = config[ParamsPreprocData.input_mean_ctr_int_values_file]

    output_radiomics_by_patient = config[ParamsPreprocData.output_radiomics_by_patient_file]
    output_clinical_and_intensities = config[ParamsPreprocData.output_merged_clinical_and_mean_ctr_values_file]
    output_merged_clinical_radiomics = config[ParamsPreprocData.output_merged_clinical_and_radiomics]
    output_merged_and_normalized = config[ParamsPreprocData.output_normalized_clinical_and_radiomics]
    # This is the input to create all the subfiles with filtered columns
    input_temp_merged_clinical_radiomics = output_merged_clinical_radiomics
    NAN_VALUE = '?'

    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)

    # -------------------- Merging Radiomics -------------------
    print('****** Merging Radiomics table, output: {}...'.format(output_radiomics_by_patient))
    mergeRadiomicsByBiopsy(input_radiomics, output_radiomics_by_patient)
    print('Done!')

    # -------------------- Merging with normalization columns (femur, muscle, fat)-------------------
    print('****** Joining main with normalization')
    print(F'{input_biopsy_clinical_info} with {input_mean_intensities}, output: {output_clinical_and_intensities}')
    mergeNormalization(input_biopsy_clinical_info, input_mean_intensities, output_clinical_and_intensities)
    print('Done!')

    # -------------------- Merging to MASTER TABLE -------------------
    print('****** Joining with master table, output: {}...'.format(output_merged_clinical_radiomics))
    joinTablesByID(output_clinical_and_intensities, output_radiomics_by_patient, output_merged_clinical_radiomics)
    print('Done!')

    # %% -------------------- Normalizing MASTER TABLE -------------------
    print('****** Normalizing master table in, output: {}...'.format(output_merged_and_normalized))
    mainNormalizeIntensities(output_merged_clinical_radiomics, output_merged_and_normalized)

    print('Done!')

    exit()

    # norm_types = ['PZTZ', 'None', 'Muscle']
    norm_types = ['None']

    # ------------------- Cleaning columns -------------------
    print('Cleaning columns....')
    mergedData = pd.read_csv(input_temp_merged_clinical_radiomics,  index_col=0)
    allColumns = mergedData.keys()

    # ------- Select variables that we want to keep ------
    keep = getClinicalVariables() + getImagingVariables(allColumns) + getOtherVariables() # Default
    # keep =  getClinicalVariables()
    # keep =  getImagingVariables(allColumns)
    # keep =  getClinicalVariables() + getImagingVariables(allColumns)

    keepIdx= [np.where(allColumns == x)[0][0]  for x in keep]

    # classes = [ 'GS_um', 'total_gleason_score', 'NCCN(form)', 'Decipher', 'SUM (5 tier)', 'SUM (3 tier)' ]
    # titles = [ 'UM', 'GenDX_GG', 'NCCN', 'DecipherNCCN', '6-Tier', '3-Tier']
    classes = ['SUM (3 tier)' ]
    titles = ['3-Tier']

    removeColumnsIdxOriginal = np.arange(0,len(allColumns)) # Adding all columns
    # Remove the 'keep' columns in the list of columns to be deleted
    removeColumnsIdxOriginal = removeColumnsIdxOriginal[ np.logical_not(np.isin(removeColumnsIdxOriginal, keepIdx)) ]
    # Remove everything after 'GM_t2_con_Kur' (the empty muscles)
    gm_t2_con_idx = np.where(allColumns == 'GM_t2_con_Kur')[0][0]
    print('Deleting all columns after GM_t2_con_Kur: ',gm_t2_con_idx)
    removeColumnsIdxOriginal = np.append(removeColumnsIdxOriginal, np.arange(gm_t2_con_idx, mergedData.shape[1]))

    makeHeatMaps = False

    for onlyLesion in [True, False]:
        print('-------------------------- Only lesion {} -----------------------'.format(onlyLesion))
        # This loop iterates over the normalization types of the intensities
        for norm_type in norm_types:
            print('########################## Normalization {} #######################'.format(norm_type))
            if norm_type == 'PZTZ':
                titleNorm = 'Normalize (PZ+TZ)/2'
                pre='NormByPZ_and_TZ'
            if norm_type == 'None':
                titleNorm = 'No normalization'
                pre='NO_norm'
            if norm_type == 'Muscle':
                titleNorm = 'Normalize Muscle'
                pre='NormByMuscle'
            if onlyLesion:
                pre = 'OnlyLesion_{}'.format(pre)
            else:
                pre = 'All_{}'.format(pre)

            outFiles=[join(output_folder,'Normal_{}_{}'.format(pre,title)) for title in titles]

            # Iterate over all the final files
            for idx, fileName in enumerate(outFiles):
                t1 = np.where(allColumns == classes[idx])[0][0]
                t2 = np.where(removeColumnsIdxOriginal == t1)[0][0]
                removeColumnsIdx = np.delete(removeColumnsIdxOriginal, t2) # Add the class column
                removeColumnsNames= [allColumns[x] for x in removeColumnsIdx]

                # Removing all columns that have at least 3 nan
                # allData = allData.dropna(axis='columns', thresh=3)
                allData = setClassCol(mergedData, classes[idx])

                cleanedData = normalizeData(allData, norm_type)

                # Making the final files by merging GSC depending on the source
                for merged in (True,False):
                    data, finalFileName = mergeGSC(cleanedData, merged, titles[idx], fileName)
                    print('*****{} ****'.format(finalFileName))
                    print('Removing columns: ',removeColumnsNames)
                    data= data.drop(removeColumnsNames, axis=1)
                    print('Reseting index')
                    data= data.reset_index(drop=True)
                    # data = removeMuscle(data)
                    data = keepOnlyLesion(data, onlyLesion)
                    data.to_csv('{}.csv'.format(finalFileName), na_rep=NAN_VALUE, index=False)

                # if makeHeatMaps:
                #     if (idx == 0) and (merged): # Only make maps for 'merged' and 'UM' (should be the same)
                #         if onlyLesion:
                #             plotHeatmaps(cleanedData, output_folder, pre, figsize=35)
                #         else:
                #             plotHeatmaps(cleanedData, output_folder, pre, figsize=100)

