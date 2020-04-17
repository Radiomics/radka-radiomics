"""
Email from Radka, Feb 3, 2020:

(1) (global) We wanted to normalize the image intensities in the ROI (in our case – the entire prostate) between 0 – 255 and calculate the texture features in 128 bins in moving box of 5x5x5.

(2) (local) Calculate the texture features in 128 bins in moving box of 5x5x5, where the intensities are normalized between 0 – 255 in the 5x5x5 box.
"""

from radiomics import firstorder, shape, glcm, logger
import os, sys, json, six
import numpy
import SimpleITK as sitk
import pandas as pd


def pprint_dict(input_dict):
  print(json.dumps(input_dict, indent=2))

# inject functions from https://github.com/Radiomics/pyradiomics/pull/337,
# adding imports

def _getKernelGenerator(self):
  from radiomics import cMatrices

  kernelRadius = self.settings.get('kernelRadius', 1)

  ROI_mask = sitk.GetArrayFromImage(self.inputMask) == self.label
  ROI_indices = numpy.array(numpy.where(ROI_mask))

  # Get the size of the input, which depends on whether it is in masked mode or not
  if self.masked:
    size = numpy.max(ROI_indices, 1) - numpy.min(ROI_indices, 1) + 1
  else:
    size = numpy.array(self.imageArray.shape)

  # Take the minimum size along each x, y and z dimension from either the size of the ROI or the kernel
  # First add the kernel radius to the size, yielding shape (2, 3), then take the minimum along axis 0, getting back
  # to shape (3,)
  self.boundingBoxSize = numpy.min(numpy.insert([size], 1, kernelRadius * 2 + 1, axis=0), axis=0)

  # Calculate the offsets, which are used to generate a list of kernel Coordinates
  kernelOffsets = cMatrices.generate_angles(self.boundingBoxSize,
                                            numpy.array(six.moves.range(1, kernelRadius + 1)),
                                            True,  # Bi-directional
                                            self.settings.get('force2D', False),
                                            self.settings.get('force2Ddimension', 0))

  # Generator loop that yields a kernel mask: a boolean array that defines the voxels included in the kernel
  kernelMask = numpy.zeros(self.imageArray.shape, dtype='bool')  # Boolean array to hold mask defining current kernel

  for idx in ROI_indices.T:  # Flip axes to get sets of 3 elements (z, y and x) for each voxel
    kernelMask[:] = False  # Reset kernel mask

    # Get coordinates for all potential voxels in this kernel
    kernelCoordinates = kernelOffsets + idx

    # Exclude voxels outside image bounds
    kernelCoordinates = numpy.delete(kernelCoordinates, numpy.where(numpy.any(kernelCoordinates < 0, axis=1)), axis=0)
    kernelCoordinates = numpy.delete(kernelCoordinates,
                                       numpy.where(numpy.any(kernelCoordinates >= self.imageArray.shape, axis=1)), axis=0)

    idx = tuple(idx)

    # Transform indices to boolean mask array
    kernelMask[tuple(kernelCoordinates.T)] = True
    kernelMask[idx] = True  # Also include center voxel

    if self.masked:
      # Exclude voxels outside ROI
      kernelMask = numpy.logical_and(kernelMask, ROI_mask)

      # check if there are enough voxels to calculate texture, skip voxel if this is not the case.
      if numpy.sum(kernelMask) <= 1:
        continue

    # Also yield the index, identifying which voxel this kernel belongs to
    yield idx, kernelMask


# redefine _calculateVoxels
import copy
def _calculateVoxels(self):
  initValue = self.settings.get('initValue', 0)
  self.kernels = self._getKernelGenerator()

  # per-kernel normalization
  if SETTING_normalizationType == "global":
    # Radka option 1: global binning + voxel feature calculation with radius 2 kernel
    self.imageArray = sitk.GetArrayFromImage(self.inputImage) # Non-discretized image values
    self.imageArray = self._applyBinning(self.imageArray)
  elif SETTING_normalizationType == "kernel":
    # Radka option 2: initialize and bin kernel within the inner loop
    input_imageArray = sitk.GetArrayFromImage(self.inputImage)
  else:
    print("Invalid normalizationType: should be \"global\" or \"kernel\", but is {}!".format(SETTING_normalizationType))
    assert(False)

  # Initialize the output with empty numpy arrays
  for feature, enabled in six.iteritems(self.enabledFeatures):
    if enabled:
      self.featureValues[feature] = numpy.full(self.imageArray.shape, initValue, dtype='float')

  # Calculate the feature values for all enabled features
  with self.progressReporter(self.kernels, 'Calculating voxels') as bar:
    for vox_idx, kernelMask in self.kernels:
      self.maskArray = kernelMask
      self.labelledVoxelCoordinates = numpy.where(self.maskArray)

      if SETTING_normalizationType == "kernel":
        # Radka option 2: binning within individual radius 2 kernels
        self.imageArray = self._applyBinning(input_imageArray)

      # Calculate the feature values for the current kernel
      for success, featureName, featureValue in self._calculateFeatures():
        if success:  # Do not store results in case of an error
          self.featureValues[featureName][vox_idx] = featureValue

  # Convert the output to simple ITK image objects
  for feature, enabled in six.iteritems(self.enabledFeatures):
    if enabled:
      self.featureValues[feature] = sitk.GetImageFromArray(self.featureValues[feature])
      self.featureValues[feature].CopyInformation(self.inputImage)


def process(args_in):
  import argparse, logging

  parser = argparse.ArgumentParser(usage="%(prog)s input_mask input_image")
  parser.add_argument("-i", "--image", required=True, help="Input image for feature extraction", dest="input_image")
  parser.add_argument("-m", "--mask", required=True, help="Input image for feature extraction", dest="input_mask")
  parser.add_argument("-n", "--normalization", choices=["kernel","global"], dest="normalization_type", default="kernel")
  parser.add_argument("-b", "--bins", type=int, dest="bin_count", default=16)
  parser.add_argument("-l", "--mask_label", type=int, dest="mask_label", default=1)
  parser.add_argument("-r", "--kernel_radius", type=int, dest="kernel_radius", default=2)
  parser.add_argument("-t", "--force2D", type=bool, dest="force_2D", default=True)
  parser.add_argument("-k", "--masked_kernel", type=bool, dest="masked_kernel", default=True)
  parser.add_argument("-f", "--feature", type=str, default="JointEntropy", dest="feature")

  args = parser.parse_args()

  from radiomics import base
  base.RadiomicsFeaturesBase._calculateVoxels = _calculateVoxels
  base.RadiomicsFeaturesBase._getKernelGenerator = _getKernelGenerator

  # based on helloVoxel.py example from pyradiomics
  parameters = \
  {
    "imageType": {
      "Original": {}
    },
    "featureClass": {
      "glcm": [args.feature]
    },
    "setting": {
      "binCount": args.bin_count,
      "force2D": args.force_2D,
      "label": args.mask_label
    },
    "voxelSetting": {
      "kernelRadius": args.kernel_radius,
      "maskedKernel": args.masked_kernel,
      # Slicer does not like nan's - changed the default to 0!
      "initValue": 0
    }
  }

  print("Normalization: %s" % args.normalization_type)

  pprint_dict(parameters)

  try:
    image = sitk.ReadImage(args.input_image)
    mask = sitk.ReadImage(args.input_mask)
  except:
    print("Fatal error: failed to read input image or mask")
    return

  import six, numpy
  from radiomics import featureextractor, getFeatureClasses

  logger.setLevel(logging.DEBUG)

  global SETTING_normalizationType
  SETTING_normalizationType = args.normalization_type # "kernel" or "global"
  extractor = featureextractor.RadiomicsFeatureExtractor(parameters)

  featureClasses = getFeatureClasses()

  for cls, features in six.iteritems(extractor.enabledFeatures):
    if features is None or len(features) == 0:
      features = [f for f, deprecated in six.iteritems(featureClasses[cls].getFeatureNames()) if not deprecated]
    for f in features:
      print(getattr(featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

  print("Calculating features")

  featureVector = extractor.execute(image, mask, voxelBased=True)

  for featureName, featureValue in six.iteritems(featureVector):
    if isinstance(featureValue, sitk.Image):
      featureFileName = "%s_%s_%s.nrrd" % (os.path.split(args.input_image)[1], featureName, args.normalization_type)

      sitk.WriteImage(featureValue, featureFileName)
      print('Computed %s, stored as "%s"' % (featureName, featureFileName))
    else:
      print('%s: %s' % (featureName, featureValue))

if __name__ == '__main__':
  process(sys.argv)
