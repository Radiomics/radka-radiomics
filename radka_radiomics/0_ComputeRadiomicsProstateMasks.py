import numpy as np
import re
import time
import SimpleITK as sitk
from img_viz.medical import *
from os.path import join
from os import listdir
from radiomics import featureextractor
import radiomics
from radiomics.imageoperations import *
from scipy.stats import kurtosis, skew
from preprocUtils.preproc.UtilsPreproc import resample_to_reference_itk, get_mask_volume_itk
from preprocUtils.preproc.UtilsItk import copyItkImage, save_image
import pandas as pd
import logging
import coloredlogs
from config.MainConfig import config_compute_texture_features
from constants.ConfigParams import *

# Gets the configuration from the user
config = config_compute_texture_features()
logger = logging.getLogger('OZ')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
logger.addHandler(ch)
fmt = '%(asctime)s,%(msecs)03d %(name)s %(levelname)s %(message)s'
coloredlogs.install(level=config[ParamsComputeRadiomics.logging_level], logger=logger, fmt=fmt)


def main():
    log_file = 'LOG.txt'
    handler = logging.FileHandler(filename=log_file, mode='w')  # overwrites log_files from previous runs. Change mode to 'a' to append.
    formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")  # format string for log messages
    handler.setFormatter(formatter)
    radiomics.logger.addHandler(handler)

    # Control the amount of logging stored by setting the level of the logger. N.B. if the level is higher than the
    # Verbositiy level, the logger level will also determine the amount of information logger.infoed to the output
    radiomics.logger.setLevel(logging.INFO)

    createRadiomicsTexturesByContour(config)


def read_imgs_and_ctrs_from_folder(input_folder, c_folder, img_names, ctr_names_input, process_all_lesions,
                                   interpolator=sitk.sitkNearestNeighbor):
    '''
    Reads all the imgs and contours desired from the current folder
    :param input_folder:  Root folder where all the cases are
    :param c_folder:   Current case folder
    :param img_names:   Which image names to read ('img_t2_tra.nrrd' for example)
    :param ctr_names_input: Which cotr names to read ('ctr_pro.nrrd' for example)
    :param process_all_lesions:  Boolean Indicates if we try to read lesions (starting with G or M)
    :param interpolator:
    :return:
    '''
    np_imgs = []
    itk_imgs = []
    np_ctrs = []
    itk_ctrs = []
    ctr_names_by_file = []
    logger.info(' Reading images  and contours....')
    all_files = listdir(join(input_folder,c_folder))

    # Reading all contour files
    ctr_names = ctr_names_input.copy()
    if process_all_lesions:
        logger.info(' Reading all the lesions (all files that start with ctr_G)')
        for c_file in all_files: # Adding contours
            # is_found = re.search('^ctr_M', c_file) or re.search('^ctr_G', c_file)
            is_found = re.search('^ctr_G', c_file)
            if is_found:
                ctr_names.append(c_file)

    # Adding all images (t2, adc, bval)
    for c_img_name in img_names:
        try:
            file_name = join(input_folder, c_folder, c_img_name)
            logger.debug(F'\tReading {file_name}')
            c_img_type_itk = sitk.ReadImage(file_name) # Reads the image
            np_imgs.append(sitk.GetArrayFromImage(c_img_type_itk)) # Gets its numpy array
            itk_imgs.append(c_img_type_itk)

            for cur_ctr_name in ctr_names:
                file_name = join(input_folder, c_folder, cur_ctr_name)
                c_ctr_itk = sitk.ReadImage(file_name) # Reads the image
                logger.debug(F'\tReading {file_name} and resampling ctr {c_ctr_itk.GetSize()} to img {c_img_type_itk.GetSize()}')
                temp_ctr = resample_to_reference_itk(c_ctr_itk, c_img_type_itk, interpolator, 0)
                this_ctr_vol = get_mask_volume_itk(temp_ctr)
                if this_ctr_vol == 0:
                    logger.warning(F' WARNING: contour {c_img_name}--{cur_ctr_name} has volume 0 (resampling {c_img_type_itk.GetSize()} with {c_ctr_itk.GetSize()}!!!')
                itk_ctrs.append(sitk.Cast(temp_ctr, sitk.sitkUInt8))
                np_ctrs.append(sitk.GetArrayFromImage(temp_ctr))
                ctr_names_by_file.append(cur_ctr_name)

        except Exception as e:
            logger.error(F'Failed reading data for {c_folder} and {c_img_name} : Error: {e}')
    logger.info('** Done Reading!')
    return itk_imgs, itk_ctrs, np_imgs, np_ctrs, ctr_names, ctr_names_by_file


def createRadiomicsTexturesByContour(config):
    """
    Iterates over all the Cases in the input folder and computes the requested radiomic textures
    :param input_folder:
    :param output_radiomics_folder:
    :param output_images_folder:
    :param img_names:
    :param out_names:
    :param ctr_names_input:
    :param radiomics_textures:
    :param radiomics_names:
    :return:
    """
    input_folder = config[ParamsComputeRadiomics.input_folder]
    output_radiomics_folder = config[ParamsComputeRadiomics.output_folder]
    output_images_folder = config[ParamsComputeRadiomics.output_imgs_folder]

    img_names = config[ParamsComputeRadiomics.img_names]
    out_names = config[ParamsComputeRadiomics.img_output_names]
    ctr_names_input = config[ParamsComputeRadiomics.ctr_names]

    generate_images = config[ParamsComputeRadiomics.generate_images]
    process_all_lesions = True # If we want to compute the radiomics for all the lesions

    # Gets the name of each of the texture features
    radiomics_textures = [name.value for name in list(config[ParamsComputeRadiomics.radiomics_textures].keys())]
    radiomics_names = list(config[ParamsComputeRadiomics.radiomics_textures].values())

    # Normalize to 0 to 255
    # params = {'binWidth': 16}
    # distances 2
    params = {'binWidth': config[RadiomicsParams.BIN_WIDTH],
              'kernelRadius': 2}

    viz_obj = MedicalImageVisualizer(output_folder=output_images_folder,
                                     disp_images=config[ParamsComputeRadiomics.show_images])

    # This is the order on how the columns will be stored on the output file
    # For future projects, you need to be carefull that the name of the columns here, correspond to the column
    # names being generated by the program
    columns_order = ['Patient', 'Contour', 'Pros_Vol', 'PZ_Vol', 'ROI_Vol', 't2_int_10', 't2_int_25', 't2_int_50', 't2_int_mean', 't2_int_75', 't2_int_90', 't2_int_SD', 't2_int_Kur', 't2_int_Ske', 'adc_int_10', 'adc_int_25', 'adc_int_50',
                     'adc_int_mean', 'adc_int_75', 'adc_int_90', 'adc_int_SD', 'adc_int_Kur', 'adc_int_Ske', 'b_int_10', 'b_int_25', 'b_int_50', 'b_int_mean', 'b_int_75', 'b_int_90', 'b_int_SD', 'b_int_Kur', 'b_int_Ske', 't2_con_10', 't2_con_25', 't2_con_50', 't2_con_mean', 't2_con_75',
                     't2_con_90', 't2_con_SD', 't2_con_Kur', 't2_con_Ske', 'adc_con_10', 'adc_con_25', 'adc_con_50', 'adc_con_mean', 'adc_con_75', 'adc_con_90', 'adc_con_SD', 'adc_con_Kur', 'adc_con_Ske', 'b_con_10', 'b_con_25', 'b_con_50', 'b_con_mean', 'b_con_75', 'b_con_90', 'b_con_SD', 'b_con_Kur', 'b_con_Ske', 't2_cor_10', 't2_cor_25', 't2_cor_50', 't2_cor_mean',
                     't2_cor_75', 't2_cor_90', 't2_cor_SD', 't2_cor_Kur', 't2_cor_Ske', 'adc_cor_10', 'adc_cor_25', 'adc_cor_50', 'adc_cor_mean', 'adc_cor_75', 'adc_cor_90', 'adc_cor_SD', 'adc_cor_Kur', 'adc_cor_Ske', 'b_cor_10', 'b_cor_25', 'b_cor_50', 'b_cor_mean', 'b_cor_75', 'b_cor_90', 'b_cor_SD', 'b_cor_Kur',
                     'b_cor_Ske', 't2_ene_10', 't2_ene_25', 't2_ene_50', 't2_ene_mean', 't2_ene_75', 't2_ene_90', 't2_ene_SD', 't2_ene_Kur', 't2_ene_Ske', 'adc_ene_10', 'adc_ene_25', 'adc_ene_50', 'adc_ene_mean', 'adc_ene_75', 'adc_ene_90', 'adc_ene_SD', 'adc_ene_Kur', 'adc_ene_Ske', 'b_ene_10', 'b_ene_25', 'b_ene_50', 'b_ene_mean', 'b_ene_75', 'b_ene_90', 'b_ene_SD', 'b_ene_Kur', 'b_ene_Ske',
                     't2_ent_10', 't2_ent_25', 't2_ent_50', 't2_ent_mean', 't2_ent_75', 't2_ent_90', 't2_ent_SD', 't2_ent_Kur', 't2_ent_Ske', 'adc_ent_10', 'adc_ent_25', 'adc_ent_50', 'adc_ent_mean', 'adc_ent_75', 'adc_ent_90', 'adc_ent_SD', 'adc_ent_Kur', 'adc_ent_Ske', 'b_ent_10', 'b_ent_25', 'b_ent_50', 'b_ent_mean',
                     'b_ent_75', 'b_ent_90', 'b_ent_SD', 'b_ent_Kur', 'b_ent_Ske', 't2_hom_10', 't2_hom_25', 't2_hom_50', 't2_hom_mean', 't2_hom_75', 't2_hom_90', 't2_hom_SD', 't2_hom_Kur', 't2_hom_Ske', 'adc_hom_10', 'adc_hom_25', 'adc_hom_50', 'adc_hom_mean', 'adc_hom_75', 'adc_hom_90', 'adc_hom_SD', 'adc_hom_Kur', 'adc_hom_Ske', 'b_hom_10', 'b_hom_25', 'b_hom_50', 'b_hom_mean', 'b_hom_75', 'b_hom_90', 'b_hom_SD', 'b_hom_Kur',
                     'b_hom_Ske']

    extractor = featureextractor.RadiomicsFeatureExtractor(**params)
    extractor.disableAllFeatures()
    extractor.enableFeaturesByName(glcm=radiomics_textures)
    logger.info(F'Enabled features: {extractor.enabledFeatures}')

    all_patients = listdir(input_folder)
    all_patients.sort()


    all_data = pd.DataFrame()
    # Iterates over all the cases
    for c_folder in all_patients:
        try:
            logger.info('*** {} ***'.format(c_folder))
            # Reading the images and contours for this case
            itk_imgs, itk_ctrs, np_imgs, np_ctrs, ctr_names, ctr_names_by_file = \
                read_imgs_and_ctrs_from_folder(input_folder, c_folder, img_names, ctr_names_input, process_all_lesions)

            # It always computes the prostate and pz volumes
            file_name = join(input_folder, c_folder, 'ctr_pro.nrrd')
            prostate_itk = sitk.ReadImage(file_name) # Reads the image
            prostate_vol = get_mask_volume_itk(prostate_itk)
            file_name = join(input_folder, c_folder, 'ctr_pz.nrrd')
            pz_itk = sitk.ReadImage(file_name) # Reads the image
            pz_vol = get_mask_volume_itk(pz_itk)

            tot_ctrs = len(ctr_names)
            # Iterates over all the contours
            for c_ctr_name_idx in range(tot_ctrs): # Iterate over all the contours (Prostate, PZ, Biopsies, etc)
                cur_row = {}
                cur_row['Patient'] = c_folder
                cur_row['Pros_Vol'] = prostate_vol
                cur_row['PZ_Vol'] = pz_vol
                # From the ctr_names it removes ctr_ and .nrrd, what is left is used for the column names
                cur_row['Contour'] = ctr_names[c_ctr_name_idx].upper().replace('CTR_','').replace('.NRRD','')

                # Iterate over images (T2, ADC, BVAL)
                for c_img_idx, c_img_itk in enumerate(itk_imgs):
                    try:
                        # TODO hardcoded name to adc, set is as a parameter
                        if img_names[c_img_idx].lower().find('adc') != -1:
                            c_img_np = sitk.GetArrayFromImage(c_img_itk)
                            real_min_adc_val = np.amin(c_img_np)
                            min_adc_val = max(0,real_min_adc_val)
                            # If the minimum value inside the image is less than 0 we fix the adc intensities
                            if min_adc_val <= 0:
                                # Everything that is lower than 0, will be replaced with the value at the
                                # 5th percentile of the remaining pixels. Those above min_value
                                new_low_value = np.percentile(c_img_np[c_img_np > min_adc_val], 5)
                                logger.debug(F'Fixing ADC. Replacing low values < min(0,{real_min_adc_val}) with the 5th percentile {new_low_value}')
                                c_img_np[c_img_np <= min_adc_val] = new_low_value
                                c_img_itk = copyItkImage(c_img_itk, c_img_np)

                        # Because we resampled every contour to the image, they are 'repeated'
                        c_ctr_img_idx = c_img_idx*tot_ctrs+c_ctr_name_idx
                        c_ctr_itk = itk_ctrs[c_ctr_img_idx]
                        roi_volume = get_mask_volume_itk(c_ctr_itk)
                        cur_row['ROI_Vol'] = roi_volume

                        # logger.debug(F'\r\t ---Processing: {img_names[c_img_idx]}{c_img_itk.GetSize()}{c_img_itk.GetPixelIDTypeAsString()}'
                        logger.debug(F'\r\t ---Processing: {img_names[c_img_idx]}{c_img_itk.GetSize()} -- {ctr_names_by_file[c_ctr_img_idx]}{c_ctr_itk.GetSize()}'
                                     F'  IDs:{c_img_idx} and {c_img_idx*tot_ctrs+c_ctr_name_idx}---')

                        # utilsviz.drawSeriesItk(c_img_itk, slices='middle', title=img_names[c_img_idx], contours=[c_ctr_itk])
                        logger.debug('\t Processing intensities (cropped mask and img)....')
                        # tit = F'{img_names[c_img_idx]} -- {ctr_names_by_file[c_ctr_img_idx]}'
                        # utilsviz.view_results = True
                        # utilsviz.drawSeriesItk(c_img_itk, slices='all', title=tit, contours=[c_ctr_itk])
                        bbox = checkMask(c_img_itk, c_ctr_itk)
                        img_crop_itk, ctr_crop_itk = cropToTumorMask(c_img_itk, c_ctr_itk, bbox[0])   # Mask image with current contour
                        file_name = '{}_{}_{}'.format(out_names[c_img_idx], 'int', ctr_names[c_ctr_name_idx])

                        # Saves the masked ctr of the image (intensities)
                        img_crop_np = sitk.GetArrayFromImage(img_crop_itk)
                        ctr_crop_np = sitk.GetArrayFromImage(ctr_crop_itk)
                        # Saving the cropped contour with original intensities
                        img_masked = sitk.GetImageFromArray(ctr_crop_np*img_crop_np)
                        save_image(img_masked, join(output_radiomics_folder, c_folder), file_name)

                        # Do not compute the textures for the prostate contour
                        # if ctr_names[c_ctr_name_idx] == 'ctr_pro.nrrd':
                        #     logger.debug('Not computing textures for the prostate!')
                        #     continue

                        logger.debug('Processing textures....')
                        t = time.time()
                        result = extractor.execute(img_crop_itk, ctr_crop_itk, voxelBased=True,)
                        logger.debug(F' Done! Total time: {(time.time() - t)}')

                        logger.debug(F'\t Computing statistic for intensities...')
                        cur_row = computeStatistics(sitk.GetArrayFromImage(img_crop_itk),
                                                    sitk.GetArrayFromImage(ctr_crop_itk), img_name=img_names[c_img_idx],
                                                    texture_name='int', cur_row_dict=cur_row)

                        # Iterate over all the Radiomic features
                        logger.debug('\t Computing statistics for all textures, and saving files....')
                        for c_texture_idx,c_radio in enumerate(radiomics_textures): # Adding all images
                            # logger.debug(F'\t\t Texture {c_radio}')
                            # logger.debug('\t\t ___ {} ___'.format(c_radio))
                            # logger.debug(result.keys())
                            res_img = result['original_glcm_{}'.format(c_radio)]  # Retrieves the proper texture feature
                            file_name = F'{out_names[c_img_idx]}_{radiomics_names[c_texture_idx]}_{ctr_names[c_ctr_name_idx]}'

                            if generate_images:
                                viz_obj.output_folder = join(output_images_folder, c_folder)
                                file_prefix = file_name.replace('.nrrd','')
                                viz_obj.plot_img_and_ctrs_itk(img_itk=res_img, slices=SliceMode.MIDDLE,
                                                              title=file_prefix,
                                                              file_name_prefix=file_prefix)

                            cur_row = computeStatistics(sitk.GetArrayFromImage(img_crop_itk),
                                                        sitk.GetArrayFromImage(ctr_crop_itk), img_name=img_names[c_img_idx],
                                                        texture_name=radiomics_names[c_texture_idx], cur_row_dict=cur_row)

                            save_image(res_img, join(output_radiomics_folder, c_folder), file_name)

                        # logger.info(F'\t\t Statistics for {ctr_names[c_ctr_name_idx]}: {cur_row}')
                        logger.debug('\t Done!....')

                    except Exception as e:
                        logger.error(F' ############# ERROR ({img_names[c_img_idx]}--{ctr_names_by_file[c_ctr_img_idx]} ){e} ########### ')

                if ctr_names[c_ctr_name_idx] == 'ctr_pro.nrrd': # Skip everything when the contour is the prostate
                    logger.info('This is the prostate contour, not saving the results into the csv file')
                    continue

                all_data = all_data.append(cur_row, ignore_index=True)
                all_data = all_data[columns_order]
                all_data.to_csv(join(output_radiomics_folder,'Radiomics.csv'), index=False)
                # logger.debug(F"{all_data[['Patient','Contour']].values}")

        except Exception as e:
            logger.info(F'\t ### ERROR Patient {c_folder} failed: {e}')
    exit()


def computeStatistics(img_np, mask_np, img_name='', texture_name='', cur_row_dict={}):
    mask_indexes = mask_np > 0
    img_values = img_np[mask_indexes]

    final_col_name = ''
    if img_name.find('t2') != -1:
        final_col_name += 't2_'
    if img_name.find('adc') != -1:
        final_col_name += 'adc_'
    if img_name.find('bval') != -1:
        final_col_name += 'b_'

    stats = { F'{final_col_name}{texture_name}_10': np.percentile(img_values, 10),
              F'{final_col_name}{texture_name}_25': np.percentile(img_values, 25),
              F'{final_col_name}{texture_name}_50': np.percentile(img_values, 50),
              F'{final_col_name}{texture_name}_75': np.percentile(img_values, 75),
              F'{final_col_name}{texture_name}_90': np.percentile(img_values, 90),
              F'{final_col_name}{texture_name}_mean': np.mean(img_values),
              F'{final_col_name}{texture_name}_SD': np.std(img_values),
              F'{final_col_name}{texture_name}_Kur': kurtosis(img_values),
              F'{final_col_name}{texture_name}_Ske': skew(img_values)
              }

    cur_row_dict.update(stats)
    return cur_row_dict


if __name__ == '__main__':
    main()

