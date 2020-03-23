"""
Email from Radka, Feb 3, 2020:

(1) (global) We wanted to normalize the image intensities in the ROI (in our case – the entire prostate) between 0 – 255 and calculate the texture features in 128 bins in moving box of 5x5x5.

(2) (local) Calculate the texture features in 128 bins in moving box of 5x5x5, where the intensities are normalized between 0 – 255 in the 5x5x5 box.
"""

from radiomics import firstorder, shape, glcm
import os, sys, json, six
import numpy as np
import SimpleITK as sitk
import pandas as pd

import tqdm
import radiomics

radiomics.progressReporter = tqdm.tqdm


def _alt_binEdges(parameterValues, **kwargs):
  binCount = kwargs.get('binCount')

  minimum = min(parameterValues)
  maximum = max(parameterValues) + 1
  binEdges = np.linspace(minimum, maximum, binCount + 1)

  return binEdges


radiomics.imageoperations.getBinEdges = _alt_binEdges


def normalizeArray(nparray, new_max):
    arr_min = np.min(nparray)
    arr_max = np.max(nparray)
    return (nparray-arr_min)*(new_max / (arr_max - arr_min))

def normalizedImage(image, mask=None):
    im_arr = sitk.GetArrayFromImage(image).astype(float)
    if mask is not None:
      mask_arr = sitk.GetArrayFromImage(mask).astype(bool)
      arr_min = np.min(im_arr[mask_arr])
      arr_max = np.max(im_arr[mask_arr])

      im_arr[~mask_arr] = arr_min  # results in non-masked voxels being 0 after normalization
    else:
      arr_min = np.min(im_arr)
      arr_max = np.max(im_arr)

    im_arr -= arr_min
    im_arr *= 256 / (arr_max - arr_min)

    im = sitk.GetImageFromArray(im_arr.astype(int))
    im.CopyInformation(image)
    return im

def getTexture(image,mask,features):
    settings = {"voxelBased": True, "kernelRadius": 2, "binCount": 128}
    gf = glcm.RadiomicsGLCM(image, mask, **settings)
    for f in features:
        gf.enableFeatureByName(f, True)
    return gf.execute()

def pprint_dict(input_dict):
  print(json.dumps(input_dict, indent=2))


imageFiles = ["subject0000.nrrd","subject0001.nrrd"]
subjects = ["0000", "0001","0015","0026","0055"]
maskTypes = ["T_PZ", "F_TZ","T_TZ","T_TZ","T_TZ"]
#maskFileSuffixes = ["subject0000 ROI-1_T_PZ.nrrd", "subject0001 ROI-1_F_TZ.nrrd"]

usfValues = []
prValues = []
subjectValues = []
featureTypes = []



import os
import SimpleITK as sitk

s = "0055"
filePrefix = "../cases"
imageFile = os.path.join(filePrefix, "subject%s.nrrd" % s)
maskFile = os.path.join(filePrefix, "subject%s ROI-1_%s.nrrd" % (s, maskTypes[subjects.index(s)]))
wholeGlandMaskFile = os.path.join(filePrefix, "subject%s-WholeGland.nrrd" % s)
image = sitk.ReadImage(imageFile)

mask = sitk.ReadImage(maskFile)
wholeGlandMask = sitk.ReadImage(wholeGlandMaskFile)



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

    #print(np.sum(kernelMask))

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

  #print("Before the loop")
  #tempImageArray = sitk.GetArrayFromImage(self.inputImage)[numpy.where(self.maskArray)]
  #print(str(tempImageArray))

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

  self.progressReporter = tqdm.tqdm

  # Calculate the feature values for all enabled features
  with self.progressReporter(self.kernels, 'Calculating voxels') as bar:
    for vox_idx, kernelMask in bar:
      self.maskArray = kernelMask
      self.labelledVoxelCoordinates = numpy.where(self.maskArray)

      if SETTING_normalizationType == "kernel":
        # Radka option 2: binning within individual radius 2 kernels
        self.imageArray = self._applyBinning(input_imageArray)

      # Calculate the feature values for the current kernel
      for success, featureName, featureValue in self._calculateFeatures():
        if success:  # Do not store results in case of an error
          self.featureValues[featureName][vox_idx] = featureValue

      ### DEBUG ->
      """
      tempImage = sitk.GetImageFromArray(self.imageArray.astype(int))
      tempImage.CopyInformation(self.inputImage)
      sitk.WriteImage(tempImage, "binned_image.nrrd")
      print("Min: {} Max: {}".format(self.imageArray[numpy.where(self.maskArray)].min(), self.imageArray[numpy.where(self.maskArray)].max()))

      sitk.WriteImage(self.inputImage, "original_image.nrrd")
      print("Min: {} Max: {}".format(sitk.GetArrayFromImage(self.inputImage)[numpy.where(self.maskArray)].min(), sitk.GetArrayFromImage(self.inputImage)[numpy.where(self.maskArray)].max()))

      tempImage = sitk.GetImageFromArray(self.maskArray.astype(int))
      tempImage.CopyInformation(self.inputImage)
      sitk.WriteImage(tempImage, "mask_image.nrrd")

      assert(False)
      """
      ### DEBUG <-

  # Convert the output to simple ITK image objects
  for feature, enabled in six.iteritems(self.enabledFeatures):
    if enabled:
      self.featureValues[feature] = sitk.GetImageFromArray(self.featureValues[feature])
      self.featureValues[feature].CopyInformation(self.inputImage)


from radiomics import base
base.RadiomicsFeaturesBase._calculateVoxels = _calculateVoxels
base.RadiomicsFeaturesBase._getKernelGenerator = _getKernelGenerator



# load parameters from file (from featureexractor.py)
# TODO: to add a static method to do this?

paramsFile = os.path.abspath(r'exampleVoxel.yaml')

# Ensure pykwalify.core has a log handler (needed when parameter validation fails)
global logger
import pykwalify, logging
from radiomics import getParameterValidationFiles

if len(pykwalify.core.log.handlers) == 0 and len(logging.getLogger().handlers) == 0:
  # No handler available for either pykwalify or root logger, provide first radiomics handler (outputs to stderr)
  pykwalify.core.log.addHandler(logging.getLogger('radiomics').handlers[0])

schemaFile, schemaFuncs = getParameterValidationFiles()
c = pykwalify.core.Core(source_file=paramsFile, source_data=None,
                        schema_files=[schemaFile], extensions=[schemaFuncs])
params = c.validate()

if "binWidth" in params["setting"].keys():
  del params["setting"]["binWidth"]
if "voxelBatch" in params["voxelSetting"].keys():
  del params["voxelSetting"]["voxelBatch"]

params["setting"]["binCount"] = 16
params["featureClass"]["glcm"] = ["Contrast"]

pprint_dict(params)


# based on helloVoxel.py example
# Slicer does not like nan's - changed the default to 0!
import six, numpy
from radiomics import featureextractor, getFeatureClasses

global SETTING_normalizationType
SETTING_normalizationType = "kernel" # "kernel" or "global"
extractor = featureextractor.RadiomicsFeatureExtractor(params)

featureClasses = getFeatureClasses()

print("Active features:")
for cls, features in six.iteritems(extractor.enabledFeatures):
  if features is None or len(features) == 0:
    features = [f for f, deprecated in six.iteritems(featureClasses[cls].getFeatureNames()) if not deprecated]
  for f in features:
    print(f)
    print(getattr(featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

print("Calculating features")

featureMask = wholeGlandMask
featureVector = extractor.execute(image, featureMask, voxelBased=True)

### save feature imageShape
for featureName, featureValue in six.iteritems(featureVector):
  if isinstance(featureValue, sitk.Image):
    if featureMask == mask:
      featureFileName = os.path.join(filePrefix, "subject%s ROI-1_%s_%s_%s.nrrd" % (s, maskTypes[0], featureName, SETTING_normalizationType))
    elif featureMask == wholeGlandMask:
      featureFileName = os.path.join(filePrefix, "subject%s ROI-WholeGland_%s_%s.nrrd" % (s, featureName, SETTING_normalizationType))

    sitk.WriteImage(featureValue, featureFileName)
    print('Computed %s, stored as "%s"' % (featureName, featureFileName))
  else:
    print('%s: %s' % (featureName, featureValue))
