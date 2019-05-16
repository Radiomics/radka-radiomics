import logging

import itk
import numpy as np
import SimpleITK as sitk

import radiomics

radiomics.setVerbosity(logging.DEBUG)


def _read_itk_image_and_mask(im_path, ma_path, itk_dtype=itk.SS):
  im = itk.imread(im_path, itk_dtype)
  ma = itk.imread(ma_path, itk_dtype)

  return im, ma


def _read_simpleitk_image_and_mask(im_path, ma_path):
  s_im = sitk.ReadImage(im_path)
  s_ma = sitk.ReadImage(ma_path)

  return s_im, s_ma


def _compute_min_max(s_im, s_ma):
  im_arr = sitk.GetArrayFromImage(s_im)
  ma_arr = sitk.GetArrayFromImage(s_ma) == 1

  voxels = im_arr[ma_arr]

  max_vox = int(np.max(voxels))
  min_vox = int(np.min(voxels))

  return min_vox, max_vox


def _compute_pyrad_glcm(s_im, s_ma, binCount=16, use_alt_discretization=True):
  if use_alt_discretization:
    # Use alternative method to determine bin edges to get ITK comparable results
    def _alt_binEdges(parameterValues, **kwargs):
      binCount = kwargs.get('binCount')

      # Method 1 (Closer to PyRadiomics)
      # Copy the ROI image values and add 1 to the maximum value (Simulates ITK method)
      # Then, use numpy.histogram_bin_edges to get the bin edges (consistent with PyRadiomics behaviour)
      temp_vals = parameterValues.copy()
      temp_vals[temp_vals == max(temp_vals)] += 1
      edges = np.histogram_bin_edges(temp_vals, binCount)

      # Method 2 (Closer to ITK)
      # Compute the min and max in the image, and add 1 to the maximum
      # Then generated the edges by linearly arranging the values yields:
      # binCount + 1 values; each bin will be defined as [edge[edge_idx], edge[edge_idx + 1])
      minimum = min(parameterValues)
      maximum = max(parameterValues) + 1
      binEdges = np.linspace(minimum, maximum, binCount + 1)

      # Check that both methods yield the same edges
      assert np.allclose(edges, binEdges, rtol=1e-4)

      # Show the edges generated to the user
      print('PyRadiomics alternative bins')
      for bin_idx in range(binCount):
        print((binEdges[bin_idx], binEdges[bin_idx + 1]))

      return binEdges
    radiomics.imageoperations.getBinEdges = _alt_binEdges

  glcm_class = radiomics.getFeatureClasses()['glcm'](s_im, s_ma, binCount=binCount)

  P_glcm, angles = radiomics.cMatrices.calculate_glcm(glcm_class.imageArray,
                                                      glcm_class.maskArray,
                                                      np.array(glcm_class.settings.get('distances', [1])),
                                                      glcm_class.coefficients['Ng'],
                                                      glcm_class.settings.get('force2D', False),
                                                      glcm_class.settings.get('force2Ddimension', 0))
  P_glcm += P_glcm.copy().transpose((0, 2, 1, 3))

  return P_glcm, angles


def _compute_itk_glcm(image_to_glcm_filter, offset):

  # Calculate the GLCM
  image_to_glcm_filter.SetOffset(offset)
  image_to_glcm_filter.Update()
  itk_glcm = image_to_glcm_filter.GetOutput()

  # Convert to numpy array for easier handling
  glcm_size = tuple(itk_glcm.GetSize())
  np_glcm = np.empty(glcm_size)

  for i in range(glcm_size[0]):
    for j in range(glcm_size[1]):
      np_glcm[i, j] = itk_glcm.GetFrequency([i, j])

  return np_glcm


def main(testcase, bins=16, match_discretization=True):
  im_path, ma_path = radiomics.getTestCase(testcase)

  # Load images in itk and sitk format
  im, ma = _read_itk_image_and_mask(im_path, ma_path, itk.SS)
  s_im, s_ma = _read_simpleitk_image_and_mask(im_path, ma_path)

  # compute min and maximum in the image (needed for ITK)
  min_vox, max_vox = _compute_min_max(s_im, s_ma)

  # Get the PyRadiomics results
  P_glcm, angles = _compute_pyrad_glcm(s_im, s_ma, bins, match_discretization)

  # ITK operates in (x, y, z), PyRadiomics in (z, y, x). Therefore flip the angles
  itk_offsets = [tuple(int(da) for da in a[::-1]) for a in angles]

  # Compute the ITK GLCMs
  sitm = itk.ScalarImageToCooccurrenceMatrixFilter[im].New()
  sitm.SetInput(im)
  sitm.SetMaskImage(ma)
  sitm.SetNumberOfBinsPerAxis(bins)

  # ITK adds +1 to the maximum prior to generating the bin edges, but PyRadiomics adds this after calculation of the
  # edges. To get more comparable results when using native-pyradiomics discretization, subtract 1 here.
  # N.B. this does cause ITK to ignore the voxels with value == max_vox!
  if not match_discretization:
    sitm.SetPixelValueMinMax(min_vox, max_vox - 1)
  else:  # Alternative binning is used, matching ITK's implementation. Therefore, use default behaviour in ITK too.
    sitm.SetPixelValueMinMax(min_vox, max_vox)

  # ITK GLCM is calculated for each offset separately!
  itk_glcms = np.empty((bins, bins, len(itk_offsets)))
  for o_idx, o in enumerate(itk_offsets):
     itk_glcms[:, :, o_idx] = _compute_itk_glcm(sitm, o)

  # Get last output
  itk_glcm = sitm.GetOutput()
  print('ITK bins')
  for i in range(bins):
    bin = (itk_glcm.GetBinMin(0, i), itk_glcm.GetBinMax(0, i))
    print(bin)
  # Compute the difference between ITK and PyRadiomics GLCMs
  diff = P_glcm[0, :, :, :] - itk_glcms

  # Sum over the different offsets to show differences for i and j, across offsets
  diff_angle_sum = np.sum(np.abs(diff), axis=2)

  print('GLCM difference (summed over offsets)')
  print(diff_angle_sum)

  print('Total difference: %i' % np.sum(diff_angle_sum))


if __name__ == '__main__':
  test_case = 'brain1'
  match_discretization = True

  main(test_case, bins=16, match_discretization=match_discretization)
