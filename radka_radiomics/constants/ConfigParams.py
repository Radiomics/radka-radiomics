from enum import Enum


# ------------- These are for the First step, computing radiomics --------------
class ParamsComputeRadiomics(Enum):
    input_folder = 1  # Input folder (preprocRadiomics)
    output_folder = 2  # Where to save all the radiomics files in nrrd format (also the CSV file)
    output_imgs_folder = 3  # Where to save intermediate images
    img_names = 4  # Name of the images where we want to compute the radiomics
    img_output_names = 5  # What are the prefix name of each of the images analized
    ctr_names = 6  # Which contours to analyze (besied the lesions)
    radiomics_textures = 7  # Name of each texture feature to compute
    logging_level = 8  # Indicate the level of logging we are using

    generate_images = 9  # Boolean that indicates if we want to generate intermediate images
    show_images = 10  # Indicates if we want to display images (PyCharm)

class RadiomicFeatures(Enum):
    CONTRAST = 'Contrast'
    CORRELATION = 'Correlation'
    ENTROPY = 'JointEntropy'
    ENERGY = 'JointEnergy'
    HOMOGENEITY = 'Id'


class RadiomicsParams(Enum):
    BIN_WIDTH = 1


# ------------- These are for the Second step, normalizing, merging etc, --------------
class ParamsPreprocData(Enum):
    input_folder = 1
    output_folder = 2
    output_temp_folder = 20
    input_radiomics_file = 3  # This file s obtained on the previous step (has al the radiomics statistics)
    input_clinical_biopsy_info_file = 4  # This file contains the clinical and genetic information
    # This file has the mean values for the ctrs of interest to perform tissue normalization
    input_mean_ctr_int_values_file = 5

    output_radiomics_by_patient_file = 6  # Has the radiomics values organized with a single row by patient
    # Contains all the clinical and mean intensity ctr information merged
    output_merged_clinical_and_mean_ctr_values_file = 7
    output_merged_clinical_and_radiomics = 8  # Contains the final radiomics and clinical variables wo normalization
    output_normalized_clinical_and_radiomics = 9
