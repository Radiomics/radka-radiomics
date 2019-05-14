# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 17:20:16 2019

@author: j.v.griethuysen
"""
import os

import itk


import radiomics
import logging
radiomics.setVerbosity(logging.DEBUG)
from radiomics.featureextractor import RadiomicsFeatureExtractor
#radiomics.imageoperations.logger.addHandler(radiomics.logger.handlers[0])
#radiomics.imageoperations.logger.setLevel(logging.DEBUG)
#radiomics.imageoperations.logger.handlers[0].setLevel(logging.DEBUG)

from radiomics.glcm import RadiomicsGLCM
from radiomics import cMatrices
import numpy as np
import SimpleITK as sitk


bins= 16
root = r'E:\Git-Repos\pyradiomics\data'

test_case = 'breast1'

im_path = os.path.join(root, '%s_image.nrrd' % test_case)
ma_path = os.path.join(root, '%s_SignedShort_label.nrrd' % test_case)

if not os.path.isfile(ma_path):
    sma = sitk.ReadImage(os.path.join(root, '%s_label.nrrd' % test_case))
    sitk.WriteImage(sitk.Cast(sma, sitk.sitkInt16), ma_path, True)

im = itk.imread(im_path)
ma = itk.imread(ma_path, itk.SS)
InputType = itk.Image[itk.SS, 3]

# %%
#itk.ImageFileReader.GetTypes() # (os.path.join(root, 'brain1_label.nrrd'))

#cst = itk.CastImageFilter[itk.Image[itk.SS, 3], itk.Image[itk.F, 3]].New()
#cst.SetInput(ma)
#cst.Update()

#ma_fl = cst.GetOutput()
#ma_fl

siff = itk.ScalarImageToTextureFeaturesFilter[InputType].New()
dir(siff)

siff.SetNumberOfBinsPerAxis(128)
siff.SetInput(im)
siff.SetMaskImage(ma)
#siff.SetInsidePixelValue(1)

siff.FastCalculationsOff()

siff.Update()
res = siff.GetFeatureMeans()
siff.GetInputNames()

#itk.GetArrayFromVnlVector(res)
dir(res)
res.GetElement(0)


offset = itk.Offset[3]()
idx = 0
while siff.GetOffsets().GetElementIfIndexExists(idx, offset):
    idx += 1
    print(offset)

idx = 0
fm = siff.GetFeatureMeans()
#fm = siff.GetFeatureStandardDeviations()
while fm.GetElementIfIndexExists(idx, None):
    print(fm.GetElement(idx))
    idx += 1
    


feats = ['Energy', 'Entropy', 'InverseDifferenceMoment', 'Inertia', 'ClusterShade', 'ClusterProminence']
feature_values = {}
for f_idx,  f in enumerate(feats):
    feature_values[f] = res.GetElement(f_idx)
# %%

extractor = RadiomicsFeatureExtractor(binCount=16)
extractor.disableAllFeatures()
extractor.enableFeaturesByName(glcm=['JointEnergy', 'JointEntropy', 'Idm', 'Contrast', 'ClusterShade', 'ClusterProminence'])
extractor.addProvenance(False)


pyrad_features = extractor.execute(im_path, ma_path)
print(pyrad_features)
print(feature_values)
extractor.settings

cMatrices.generate_angles(np.array([5, 5, 5]), np.array([1]), 0, 0, 0)
# %%

s_im = sitk.ReadImage(os.path.join(root, '%s_image.nrrd' % test_case))
sitk.WriteImage(sitk.Cast(s_im, sitk.sitkFloat32), os.path.join(root, '%s_float_image.nrrd' % test_case))

s_im = sitk.ReadImage(os.path.join(root, '%s_image.nrrd' % test_case))
s_ma = sitk.ReadImage(os.path.join(root, '%s_label.nrrd' % test_case))

im_arr = sitk.GetArrayFromImage(s_im)
ma_arr = sitk.GetArrayFromImage(s_ma) == 1

voxels = im_arr[ma_arr]

binedges = np.histogram(voxels, bins)[1]

max_vox = int(np.max(voxels))
min_vox = int(np.min(voxels))

# %%
sitm = itk.ScalarImageToCooccurrenceMatrixFilter[InputType].New()
sitm.SetInput(im)
sitm.SetMaskImage(ma)
sitm.SetOffset((1, 0, 0))
sitm.SetNumberOfBinsPerAxis(bins)
sitm.SetPixelValueMinMax(min_vox, max_vox-1)
#sitm.FastCalculationsOff()

sitm.Update()
out = sitm.GetOutput()


glcm = np.empty((bins, bins))
#with tqdm.tqdm(range(256), desc=j) as bar:
for i in range(bins):
    for j in range(bins):
        glcm[i, j] = out.GetFrequency([i, j])
# %%
#glcm[:, 0]


#dir(out)

itk_bins = []

for i in range(bins):
    itk_bins.append(out.GetBinMax(0, i))

# %%

itk.Offset.GetTypes()
offser = itk.Offset[3]

itk.ScalarImageToCooccurrenceMatrixFilter.GetTypes()
InputType

# %%
glcm_class = RadiomicsGLCM(s_im, s_ma, binCount=16)

P_glcm, angles = cMatrices.calculate_glcm(glcm_class.imageArray,
                                          glcm_class.maskArray,
                                          np.array(glcm_class.settings.get('distances', [1])),
                                          glcm_class.coefficients['Ng'],
                                          glcm_class.settings.get('force2D', False),
                                          glcm_class.settings.get('force2Ddimension', 0))

P_glcm += P_glcm.copy().transpose((0, 2, 1, 3))

# %%
glcm_class.enabledFeatures = {k: True for k in ['JointEnergy', 'JointEntropy', 'Idm', 'Contrast', 'Correlation', 'ClusterShade', 'ClusterProminence']}

pyrad_features = glcm_class.execute()

itk_glcm = glcm.copy()
itk_glcm /= np.sum(itk_glcm)

# %%

diff = glcm - P_glcm[0, :, :, 12]

# %%

hist = itk.Histogram[itk.F].New(MeasurementVectorSize=2)
hist_size = itk.Array.UL(2)
hist_size.Fill(bins)

hist.Initialize(hist_size, (min_vox, min_vox), (max_vox, max_vox))

# %%

#UL_glcm = P_glcm.astype('uint64')

for i in range(bins):
    for j in range(bins):
        hist.SetFrequencyOfIndex((i, j), int(P_glcm[0, i, j, 12]))

# %%

f_calculator = itk.HistogramToTextureFeaturesFilter[hist].New()
f_calculator.SetInput(hist)

f_calculator.Update()

f_vec = {
    'Energy': f_calculator.GetEnergy(),
    'Correlation': f_calculator.GetCorrelation(),
    'HaralickCorrelation': f_calculator.GetHaralickCorrelation(),
    'Entropy': f_calculator.GetEntropy(),
    'Contrast': f_calculator.GetInertia(),
    'IDM': f_calculator.GetInverseDifferenceMoment(),
    'ClusterProminence': f_calculator.GetClusterProminence(),
    'ClusterShade': f_calculator.GetClusterShade()        
}

# %%

energy = 0
t_freq = float(hist.GetTotalFrequency())
for i in range(bins**2):
    a_freq = hist.GetFrequency(i)
    r_freq = a_freq / t_freq
    energy += r_freq ** 2
    
# %%
    
n_glcm = P_glcm[0, :, :, 12].copy()
n_glcm /= np.sum(n_glcm)
print(np.sum(n_glcm))

p_ene = np.sum((n_glcm ** 2))

# %%

px = np.sum(n_glcm, axis=0)
py = np.sum(n_glcm, axis=1)

ux = np.mean(px)
uy = np.mean(py)

sx = np.std(px)
sy = np.std(py)

i = np.array(range(bins))
j = i

h_corr = np.sum(n_glcm * i[:, None] * j[None, :])
h_corr = (h_corr - ux * uy) / (sx * sy)

# %%

u = np.sum(px * i)
v = np.sum((i - u) * (i - u) * px)
corr = np.sum((n_glcm * (i[:, None] - u) * (j[None, :] - u)) / (v ** 2))

# %%

u = np.sum(n_glcm * i[:, None])
v = np.sum(n_glcm * (i[:, None] - u) * (j[:, None] - u))

corr = np.sum((n_glcm * (i[:, None] - u) * (j[None, :] - u)) / (v ** 2))

