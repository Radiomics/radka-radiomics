import itk
import numpy as np
import SimpleITK as sitk
import tqdm

import radiomics

from itk_functions import itk_read_image_and_mask, itk_alt_compute_glcm_features


FEATURES = [
    'JointEnergy',
    'JointEntropy',
    'Idm',
    'Contrast',
    'ClusterShade',
    'ClusterProminence']

def _alt_binEdges(parameterValues, **kwargs):
  binCount = kwargs.get('binCount')

  # Method 1 (Closer to PyRadiomics)
  # Copy the ROI image values and add 1 to the maximum value (Simulates ITK method)
  # Then, use numpy.histogram_bin_edges to get the bin edges (consistent with PyRadiomics behaviour)
  #temp_vals = parameterValues.copy()
  #temp_vals[temp_vals == max(temp_vals)] += 1
  #edges = np.histogram_bin_edges(temp_vals, binCount)

  # Method 2 (Closer to ITK)
  # Compute the min and max in the image, and add 1 to the maximum
  # Then generated the edges by linearly arranging the values yields:
  # binCount + 1 values; each bin will be defined as [edge[edge_idx], edge[edge_idx + 1])
  minimum = min(parameterValues)
  maximum = max(parameterValues) + 1
  binEdges = np.linspace(minimum, maximum, binCount + 1)

  # Check that both methods yield the same edges
  #assert np.allclose(edges, binEdges, rtol=1e-4)

  # Show the edges generated to the user
  #print('PyRadiomics alternative bins')
  #for bin_idx in range(binCount):
  #  print((binEdges[bin_idx], binEdges[bin_idx + 1]))

  return binEdges


radiomics.imageoperations.getBinEdges = _alt_binEdges


def kernel_iterator(im_path, ma_path, **kwargs):
  global FEATURES

  im, ma = itk_read_image_and_mask(im_path, ma_path)
  s_im = sitk.ReadImage(im_path)

  radius = kwargs.get('kernelRadius', 1)

  # Slow slow slow function!
  buffered_region = ma.GetBufferedRegion()
  size = [buffered_region.GetIndex(3), buffered_region.GetIndex(4), buffered_region.GetIndex(5)]

  output_array = np.zeros(size + [6], dtype='float')
  with tqdm.tqdm(range(size[0]), desc='x') as bar:
    for x in bar:
      for y in range(size[1]):
        for z in range(size[2]):
          if ma.GetPixel((x, y, z)) > 0:
            # Masked pixel! Extract region
            if x - radius < 0 or y - radius < 0 or z - radius < 0 or \
              x + radius >= size[0] or y + radius >= size[1] or z >= size[2]:
              raise Exception('index out of bounds')
            eif = itk.ExtractImageFilter[im, im].New()
            eif.SetInput(im)

            reg = itk.ImageRegion[3]()
            reg.SetSize([radius * 2 + 1] * 3)
            reg.SetIndex((x - radius, y - radius, z - radius))
            eif.SetExtractionRegion(reg)
            eif.Update()

            e_im = eif.GetOutput()
            fvec = itk_alt_compute_glcm_features(e_im, None, **kwargs)
            output_array[x, y, z] = fvec

  prefix = kwargs.get('prefix', '')
  for f_idx, f in enumerate(FEATURES):
    fmap = output_array[:, :, :, f_idx].transpose((2, 1, 0))
    fmap_im = sitk.GetImageFromArray(fmap)
    fmap_im.CopyInformation(s_im)

    sitk.WriteImage(fmap_im, '%s%s_fmap.nrrd' % (prefix, f))


if __name__ == "__main__":
  #from radiomics import getTestCase
  #im_path, ma_path = getTestCase('brain1')
  #kernel_iterator(im_path, ma_path, binCount=16, kernelRadius=1, prefix='brain1_')

  im_path = '../../cases/subject0055.nrrd'
  ma_path = '../../cases/subject0055-WholeGland.nrrd'
  kernel_iterator(im_path, ma_path, binCount=16, kernelRadius=2, prefix='2_prostateX_0055-wholegland_5x5x5_')
