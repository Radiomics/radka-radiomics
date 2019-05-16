#############################################################
Comparison of USF Feature calculation to PyRadiomics (v2.1.2)
#############################################################

This repository contains a review on the code used by USF for the ProstateX challenge, and how it compares to an
equivalent extraction using PyRadiomics.

Comments on the code are made on commit `fa14007 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/commit/fa14007#r33400911>`_.

Kernel Definition
=================
Voxel-based analysis is achieved by defining kernels (in `itkHaralickTextureFeaturesImageFilter.cpp <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp>`_) and then using
native ITK implementation for calculating segment-based features. Kernels are defined by instantiating a neighborhood
iterator (`L1022 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L1022>`_),
which iterates over voxels in the mask. For each segmented voxel, a region of the image (representing the kernel) is
then extracted using ITK's ExtractImageFilter (`L1048-L1050 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L1048-L1050>`_).

Prior to extraction, the region defined by the iterator is adapted, so that no invalid indices (out-of-bounds) are
included: `L1036-L1043 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L1036-L1043>`_

Finally, features are calculated `L1085-L1092 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L1085-L1092>`_.

Feature Calculation
===================
In the USF code, features are calculated using the native ITK implementations. Comparisons between PyRadiomics and ITK
feature calculation are detailed in the ``pyradiomics`` subfolder.

Gray Value discretization
-------------------------
Besides differences in feature calculation, PyRadiomics also applies discretization differently in the case of
voxel-based extraction. In PyRadiomics, discretization is applied once, based on either all values inside the ROI (when
using a masked kernel, ``maskedKernel=True``) or on all values in the cropped image sent to the feature class (when
using full kernel, ``maskedKernel=False``). ITK, or rather the ``itkHaralickTextureFeaturesImageFilter.cpp``, applies
discretization on each kernel separately. This is due to `L1074-L1082 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L1074-L1082>`_,
which re-calculated the minimum and maximum in each selected kernel. Commenting-out these lines results in behaviour
similar to PyRadiomics (discretization once, based on all intensities within the ROI). This is achieved by
`L974-L981 <https://github.com/Radiomics/UM_WIth_Olmo_Comments/blob/master/include/itkHaralickTextureFeaturesImageFilter.cpp#L974-L981>`_,
which calculate the minimum and maximum in the entire image, and which is not used now (overridden by re-calculation in
each separate kernel).
