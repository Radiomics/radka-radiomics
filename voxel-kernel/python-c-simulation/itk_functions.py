import itk
import numpy as np
import pandas as pd

from radiomics import _cmatrices


def itk_read_image_and_mask(im_path, ma_path, itk_dtype=itk.SS):
  im = itk.imread(im_path, itk_dtype)
  ma = itk.imread(ma_path, itk_dtype)

  return im, ma


def itk_compute_min_max(im, ma):
  lsif = itk.LabelStatisticsImageFilter[im, ma].New()
  lsif.SetInput(im)
  lsif.SetLabelInput(ma)
  lsif.Update()

  min_vox = lsif.GetMinimum(1)
  max_vox = lsif.GetMaximum(1)

  return min_vox, max_vox

def itk_compute_im_min_max(im):
  mmif = itk.MinimumMaximumImageFilter[im].New()
  mmif.SetInput(im)
  mmif.Update()

  min_vox = mmif.GetMinimum()
  max_vox = mmif.GetMaximum()

  return  min_vox, max_vox


def itk_get_offsets():
  angles = _cmatrices.generate_angles(np.array([3, 3, 3]), np.array([1]), 0, 0, -1)
  itk_offsets = [tuple(int(da) for da in a[::-1]) for a in angles]
  return itk_offsets


def itk_compute_glcm(im, ma, **kwargs):
  bins = kwargs.get('binCount', 16)



  min_vox, max_vox = itk_compute_min_max(im, ma)

  itk_offsets = itk_get_offsets()

  glcms = []

  for o in itk_offsets:
    # Recreate the class to ensure the output is new too (otherwise overwrites the other histograms)
    sitm = itk.ScalarImageToCooccurrenceMatrixFilter[im].New()
    sitm.SetInput(im)
    sitm.SetMaskImage(ma)
    sitm.SetNumberOfBinsPerAxis(bins)

    sitm.SetPixelValueMinMax(int(min_vox), int(max_vox))
    sitm.SetOffset(o)
    sitm.Update()
    glcm = sitm.GetOutput()
    glcms.append(glcm)


  return glcms, itk_offsets


def itk_compute_itk_features(itk_glcms):
  itk_feature_extractor = itk.HistogramToTextureFeaturesFilter[itk_glcms[0]].New()

  f_df = pd.DataFrame()
  for idx, itk_glcm in enumerate(itk_glcms):
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

    f_vec.name = idx
    f_df = f_df.join(f_vec, how='outer')
  return f_df.mean(axis=1)


def itk_alt_compute_glcm_features(im, ma, **kwargs):
  bins = kwargs.get('binCount', 16)


  sittff = itk.ScalarImageToTextureFeaturesFilter[im].New()
  sittff.SetInput(im)

  if ma is None:
    min_vox, max_vox = itk_compute_im_min_max(im)
  else:
    min_vox, max_vox = itk_compute_min_max(im, ma)
    sittff.SetMaskImage(ma)

  sittff.SetNumberOfBinsPerAxis(bins)
  sittff.SetPixelValueMinMax(int(min_vox), int(max_vox))

  sittff.Update()

  f_vec_itk = sittff.GetFeatureMeans()


  # Default features:
  #Energy, \
  #Entropy, \
  #InverseDifferenceMoment, \
  #Inertia, \
  #ClusterShade, \
  #ClusterProminence
  f_vec = pd.Series({
    'JointEnergy': f_vec_itk.GetElement(0),
    #'Correlation': ,
    #'HaralickCorrelation': ,
    'JointEntropy': f_vec_itk.GetElement(1),
    'Idm': f_vec_itk.GetElement(2),
    'Contrast': f_vec_itk.GetElement(3),
    'ClusterShade': f_vec_itk.GetElement(4),
    'ClusterProminence': f_vec_itk.GetElement(5),
  })

  f_vec = np.array([
    f_vec_itk.GetElement(0),
    f_vec_itk.GetElement(1),
    f_vec_itk.GetElement(2),
    f_vec_itk.GetElement(3),
    f_vec_itk.GetElement(4),
    f_vec_itk.GetElement(5),
  ])

  return f_vec


def compute_radiomics(im_path, ma_path):
  import radiomics
  def _alt_binEdges(parameterValues, **kwargs):
    binCount = kwargs.get('binCount')

    # Method 1 (Closer to PyRadiomics)
    # Copy the ROI image values and add 1 to the maximum value (Simulates ITK method)
    # Then, use numpy.histogram_bin_edges to get the bin edges (consistent with PyRadiomics behaviour)
    # temp_vals = parameterValues.copy()
    # temp_vals[temp_vals == max(temp_vals)] += 1
    # edges = np.histogram_bin_edges(temp_vals, binCount)

    # Method 2 (Closer to ITK)
    # Compute the min and max in the image, and add 1 to the maximum
    # Then generated the edges by linearly arranging the values yields:
    # binCount + 1 values; each bin will be defined as [edge[edge_idx], edge[edge_idx + 1])
    minimum = min(parameterValues)
    maximum = max(parameterValues) + 1
    binEdges = np.linspace(minimum, maximum, binCount + 1)

    # Check that both methods yield the same edges
    # assert np.allclose(edges, binEdges, rtol=1e-4)

    # Show the edges generated to the user
    # print('PyRadiomics alternative bins')
    # for bin_idx in range(binCount):
    #  print((binEdges[bin_idx], binEdges[bin_idx + 1]))

    return binEdges

  radiomics.imageoperations.getBinEdges = _alt_binEdges

  from radiomics import featureextractor, glcm
  import SimpleITK as sitk

  im, ma = sitk.ReadImage(im_path), sitk.ReadImage(ma_path)

  glcm_class = glcm.RadiomicsGLCM(im, ma, binWidth=None, binCount=16)
  glcm_class.disableAllFeatures()
  glcm_class.enabledFeatures = {f: True for f in [
    'JointEnergy',
    'JointEntropy',
    'Idm',
    'Contrast',
    'ClusterShade',
    'ClusterProminence']}

  #extractor = featureextractor.RadiomicsFeatureExtractor(binWidth=None, binCount=16)
  #extractor.disableAllFeatures()
  #extractor.enableFeaturesByName(glcm=['JointEnergy', 'JointEntropy', 'Idm', 'Contrast', 'ClusterShade', 'ClusterProminence'])

  fvec = glcm_class.execute()
  return fvec


def _compute_pyrad_glcm(im_path, ma_path, binCount=16, use_alt_discretization=True):
  import radiomics
  import SimpleITK as sitk
  if use_alt_discretization:
    # Use alternative method to determine bin edges to get ITK comparable results
    def _alt_binEdges(parameterValues, **kwargs):
      binCount = kwargs.get('binCount')

      minimum = min(parameterValues)
      maximum = max(parameterValues) + 1
      binEdges = np.linspace(minimum, maximum, binCount + 1)

      return binEdges
    radiomics.imageoperations.getBinEdges = _alt_binEdges

  s_im, s_ma = sitk.ReadImage(im_path), sitk.ReadImage(ma_path)
  glcm_class = radiomics.getFeatureClasses()['glcm'](s_im, s_ma, binCount=binCount)

  P_glcm, angles = radiomics.cMatrices.calculate_glcm(glcm_class.imageArray,
                                                      glcm_class.maskArray,
                                                      np.array(glcm_class.settings.get('distances', [1])),
                                                      glcm_class.coefficients['Ng'],
                                                      glcm_class.settings.get('force2D', False),
                                                      glcm_class.settings.get('force2Ddimension', 0))
  P_glcm += P_glcm.copy().transpose((0, 2, 1, 3))

  return P_glcm, angles


def compare_glcms(itk_glcms, pyrad_glcms):
  i_glcm = []
  for a in range(len(itk_glcms)):
    # Convert to numpy array for easier handling
    glcm_size = tuple(itk_glcms[a].GetSize())
    np_glcm = np.empty(glcm_size)

    for i in range(glcm_size[0]):
      for j in range(glcm_size[1]):
        np_glcm[i, j] = itk_glcms[a].GetFrequency([i, j])
    i_glcm.append(np_glcm.copy())

  i_glcm = np.array(i_glcm).transpose((1, 2, 0))
  i_glcm = np.abs(i_glcm - pyrad_glcms[0, :, :, :])
  print(i_glcm.sum())


if __name__ == "__main__":
  from radiomics import getTestCase
  im_path, ma_path = getTestCase('brain1')
  im, ma = itk_read_image_and_mask(im_path, ma_path)

  itk_glcms, offsets = itk_compute_glcm(im, ma, binCount=16)
  p_glcms, angles = _compute_pyrad_glcm(im_path, ma_path, binCount=16)

  for a_idx, a in enumerate(angles):
    assert tuple(a[::-1]) == offsets[a_idx]

  compare_glcms(itk_glcms, p_glcms)

  fvec1 = itk_compute_itk_features(itk_glcms)

  f_vals = itk_alt_compute_glcm_features(im, ma, binCount=16)
  fvec2 = pd.Series({
    'JointEnergy': f_vals[0],
    'JointEntropy': f_vals[1],
    'Idm': f_vals[2],
    'Contrast': f_vals[3],
    'ClusterShade': f_vals[4],
    'ClusterProminence': f_vals[5],
  })


  fvec3 = compute_radiomics(im_path, ma_path)

  for k in fvec1.index:
    if k not in fvec2.index:
      continue
    print("feature %s: 1: %.4f, 2: %.4f, 3: %.4f" % (k, fvec1[k], fvec2[k], fvec3[k]))
