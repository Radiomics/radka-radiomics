
import logging

import itk
import SimpleITK as sitk

import numpy as np
import pandas as pd

import radiomics

radiomics.setVerbosity(logging.DEBUG)


def get_itk_glcm(pyrad_glcm, min_vox, max_vox):
  """
  Simulates a GLCM as if it was calculated by
  `itkScalarImageToCooccurrenceMatrixFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToCooccurrenceMatrixFilter.hxx>`_.

  However, the frequencies of the elements are set using the PyRadiomics' caluclated GLCM, thereby ensuring that the
  input data for feature calculation is identical and any differences in the results are due to differences in the
  feature formulas themselves.

  :param pyrad_glcm: PyRadiomics' calculated GLCM for a specific offset
  :param min_vox: minimum intensity value in the ROI, pre-discretization
  :param max_vox: maximum intensity value in the ROI, pre-discretization
  :return: ITK histogram (2D) representing the GLCM for a specific offset
  """
  hist = itk.Histogram[itk.D].New(MeasurementVectorSize=2)
  hist_size = itk.Array.UL(2)
  hist_size.Fill(pyrad_glcm.shape[0])

  hist.Initialize(hist_size, (min_vox, min_vox), (max_vox, max_vox))

  for i in range(pyrad_glcm.shape[0]):
    for j in range(pyrad_glcm.shape[1]):
      hist.SetFrequencyOfIndex((i, j), int(pyrad_glcm[i, j]))

  return hist


def get_feature_values(itk_glcm):
  """
  Calculate ITK feature values, based on the provided ITK GLCM Matrix.

  :param itk_glcm: ITK histogram (2D) representing the GLCM for a specific offset
  :return: dict containing features calculated by ITK (keys are the corresponding PyRadiomics feature names)
  """
  itk_feature_extractor = itk.HistogramToTextureFeaturesFilter[itk_glcm].New()
  itk_feature_extractor.SetInput(itk_glcm)
  itk_feature_extractor.Update()

  f_vec = pd.Series({
    'JointEnergy': itk_feature_extractor.GetEnergy(),
    'Correlation': itk_feature_extractor.GetCorrelation(),
    'HaralickCorrelation': itk_feature_extractor.GetHaralickCorrelation(),
    'JointEntropy': itk_feature_extractor.GetEntropy(),
    'Contrast': itk_feature_extractor.GetInertia(),
    'Idm': itk_feature_extractor.GetInverseDifferenceMoment(),
    'ClusterProminence': itk_feature_extractor.GetClusterProminence(),
    'ClusterShade': itk_feature_extractor.GetClusterShade()
  })
  return f_vec


def main(testcase, bins=16):
  im_path, ma_path = radiomics.getTestCase(testcase)

  s_im = sitk.ReadImage(im_path)
  s_ma = sitk.ReadImage(ma_path)

  im_arr = sitk.GetArrayFromImage(s_im)
  ma_arr = sitk.GetArrayFromImage(s_ma) == 1

  voxels = im_arr[ma_arr]

  max_vox = int(np.max(voxels))
  min_vox = int(np.min(voxels))

  glcm_class = radiomics.getFeatureClasses()['glcm'](s_im, s_ma, binCount=bins)

  P_glcm, angles = radiomics.cMatrices.calculate_glcm(glcm_class.imageArray,
                                                      glcm_class.maskArray,
                                                      np.array(glcm_class.settings.get('distances', [1])),
                                                      glcm_class.coefficients['Ng'],
                                                      glcm_class.settings.get('force2D', False),
                                                      glcm_class.settings.get('force2Ddimension', 0))

  glcm_class.enabledFeatures = {k: True for k in
                                ['JointEnergy', 'JointEntropy', 'Idm', 'Contrast', 'Correlation', 'ClusterShade',
                                 'ClusterProminence']}

  pyrad_features = pd.Series(glcm_class.execute())

  ux, sx, uux, ssx = recurrence_mean_sd(glcm_class.coefficients['px'][0, :, 0, 12])

  print('Differences between mean and sd calculated by Numpy and recurrent method: %g and %g' % (ux - uux, sx - ssx))

  P_glcm += P_glcm.copy().transpose((0, 2, 1, 3))

  f_df = pd.DataFrame()
  for a in range(1, P_glcm.shape[3]):
    itk_glcm = get_itk_glcm(P_glcm[0, :, :, a], min_vox, max_vox)
    f_vec = get_feature_values(itk_glcm)
    f_vec.name = a
    f_df = f_df.join(f_vec, how='outer')

  python_itk_features = pd.Series()
  python_itk_features['Correlation'] = correlation_var_squared(glcm_class)[0]
  python_itk_features['HaralickCorrelation'] = haralick_correlation(glcm_class)[0]

  itk_features = f_df.mean(axis=1)
  df = pd.DataFrame({'pyrad': pyrad_features, 'itk': itk_features, 'pyrad_changed': python_itk_features})

  print(df)


def recurrence_mean_sd(px):
  ux = px[0]
  sx = 0
  for k in range(1, px.shape[0]):
    ux_min1 = ux
    sx_min1 = sx

    ux = ux_min1 + (px[k] - ux_min1) / (k + 1)
    sx = sx_min1 + (px[k] - ux_min1) * (px[k] - ux)

  sx = sx / px.shape[0]

  uux = np.mean(px)
  ssx = np.var(px)

  return ux, sx, uux, ssx


def correlation_var_squared(glcm_class):
  """
  Alternative formula for PyRadiomics Correlation.
  Differs from PyRadiomics in 2 respects:

  - Only supports calculation for symmetrical GLCM (through use of only u_i and sigma_i and squares);
    does not result in differences in calculated values.
  - Denominator in formula is the square of variance (i.e. :math:`\sigma^4`), instead of the square of standard
    deviations (i.e. variance).

  :param glcm_class: PyRadiomics GLCM calculation class instance, after having called ``execute()``
  :return: Alternative value for correlation.
  """
  i = np.array(range(glcm_class.P_glcm.shape[1]))  # 0, 1, 2, ..., n_bins - 1
  u = np.sum(glcm_class.P_glcm * i[None, :, None, None], axis=(1, 2), keepdims=True)
  s_sq = np.sum(glcm_class.P_glcm * (i[None, :, None, None] - u) ** 2, axis=(1, 2), keepdims=True)

  corm = np.sum(glcm_class.P_glcm * (i[None, :, None, None] - u) * (i[None, None, :, None] - u), (1, 2), keepdims=True)
  corr = corm / (s_sq ** 2)

  return np.nanmean(corr, (1, 2, 3))


def haralick_correlation(glcm_class):
  """
  Alternative formula for PyRadiomics Correlation. Reflects ITK's haralickCorrelation feature
  Differs from PyRadiomics in 2 respects:

  - Only supports calculation for symmetrical GLCM (through use of only u_i and sigma_i and squares);
    does not result in differences in calculated values.
  - mean and standard deviations are based only on probability, takes no account of discretized gray values.

  :param glcm_class: PyRadiomics GLCM calculation class instance, after having called ``execute()``
  :return: Alternative value for correlation.
  """
  px = glcm_class.coefficients['px']
  py = glcm_class.coefficients['py']

  ux = np.mean(px, axis=(1, 2))
  uy = np.mean(py, axis=(1, 2))

  sx = np.std(px, axis=(1, 2))
  sy = np.std(py, axis=(1, 2))

  i = np.array(range(glcm_class.P_glcm.shape[1]))  # 0, 1, 2, ..., n_bins - 1

  h_corr = np.sum(glcm_class.P_glcm * i[None, :, None, None] * i[None, None, :, None], axis=(1, 2))
  h_corr = (h_corr - (ux * uy)) / (sx * sy)
  h_corr = np.nanmean(h_corr, 1)

  return h_corr


if __name__ == '__main__':
  test_case = 'brain1'

  main(test_case)
