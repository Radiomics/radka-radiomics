from constants.ConfigParams import *
import logging

def config_compute_texture_features():
    _root_folder = 'S:/DATA_MRI/RADIOMICS_and_MRI_Normalization'
    cur_config = {ParamsComputeRadiomics.input_folder: F'{_root_folder}/Preproc',
                  ParamsComputeRadiomics.output_folder: F'{_root_folder}/RADIOMICS_ISAAC_test2',
                  ParamsComputeRadiomics.output_imgs_folder: F'{_root_folder}/RADIOMICS_ISAAC_IMGS_del2',
                  ParamsComputeRadiomics.img_names: ['img_t2_tra.nrrd', 'img_adc.nrrd', 'img_bval.nrrd'],
                  ParamsComputeRadiomics.img_output_names:  ['t2', 'adc', 'bval'],  # name of the output image
                  ParamsComputeRadiomics.ctr_names: ['ctr_pro.nrrd','ctr_NAPZ.nrrd', 'ctr_NATZ.nrrd'],
                  # The second parameter is the output name that will be used when creating the CSV file
                  ParamsComputeRadiomics.radiomics_textures: {RadiomicFeatures.CONTRAST: 'con',
                                                              RadiomicFeatures.CORRELATION: 'cor',
                                                              RadiomicFeatures.ENTROPY: 'ene',
                                                              RadiomicFeatures.ENERGY: 'ent',
                                                              RadiomicFeatures.HOMOGENEITY: 'hom'},
                  RadiomicsParams.BIN_WIDTH: 128,
                  # It can be INFO or DEBUG. Check coloredlogs documentation
                  ParamsComputeRadiomics.logging_level: 'INFO',
                  # These parameters are together. Indicate if we want to generate intermediate images
                  ParamsComputeRadiomics.generate_images: True,
                  ParamsComputeRadiomics.show_images: False
                  }

    return cur_config


def get_config_preproc_data():
    _root_folder = '/media/osz1/DATA/Dropbox/UMIAMI/WorkUM/DataAnalisis/data/'
    _input_folder = F'{_root_folder}/INPUT'
    _output_folder = F'{_root_folder}/OUTPUT'
    _output_temp_folder = F'{_root_folder}/OUTPUT_TEMP'
    cur_config = {
        ParamsPreprocData.input_folder: _input_folder,
        ParamsPreprocData.output_folder:  _output_temp_folder,
        ParamsPreprocData.output_temp_folder: _output_folder,
        ParamsPreprocData.input_radiomics_file: F'{_input_folder}/Radiomics.csv',
        ParamsPreprocData.input_clinical_biopsy_info_file: F'{_input_folder}/Clinical_Variables.csv',
        ParamsPreprocData.input_mean_ctr_int_values_file: F'{_input_folder}/Mean_ctr_values.csv',

        ParamsPreprocData.output_radiomics_by_patient_file: F'{_output_temp_folder}/TempRadiomicsByPatient.csv',
        ParamsPreprocData.output_merged_clinical_and_mean_ctr_values_file: F'{_output_temp_folder}/TempClinicalAndMeanIntensities.csv',
        ParamsPreprocData.output_merged_clinical_and_radiomics: F'{_output_folder}/Merged_Clinical_Radiomics.csv',
        ParamsPreprocData.output_normalized_clinical_and_radiomics: F'{_output_folder}/Normalized_Merged_Clinical_Radiomics.csv',
    }

    return cur_config
