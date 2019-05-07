##################################
Comparison to PyRadiomics (v2.1.2)
##################################

Pre-Processing
==============

Normalization
-------------


Gray value discretization
-------------------------

Method
++++++
In PyRadiomics, default behaviour for gray value discretization is to specify a fixed bin width, thereby ensuring that a
similar amount of contrast (i.e. width of the bin, difference in gray value that is considered a significant difference
and not just noise). In ITK, discretization is forced to a fixed bin number. PyRadiomics also supports a fixed bin number,
through parameter ``binCount``. Both ITK and PyRadiomics use this number to calculate the edges of the bins, which are used
to determine the discretized gray value (= bin number). However, there is a small difference is the method used.

This is due to how the binedges are treated: Each bin has a lower and upper edge, which define a half-open range of values
that will belong to the bin: :math:`[lower, upper)`. If bins are divided equally between minimum and maximum gray value in the image,
this would result in the topmost bin to have a range of :math:`[maximum - x, maximum)`, meaning that the maximum value would fall
outside of this bin. ITK prevents this by adding +1 to the maximum value BEFORE calculating the edges, PyRadiomics prevents this by
adding +1 to the topmost edge, resulting in a range of :math:`[maximum -x, maximum + 1)`. The rationale here is that PyRadiomics is
designed to also work with images of ``float`` datatype, which may have a range of gray values <<1 (e.g. ADC values). If
then the ITK method is used for edge calculation, a very skewed distribution would result, possibly with all values
assigned to the first 1 or 2 bins.

Application
+++++++++++
Besides a difference in calculation of the bin edges, PyRadiomics also applies discretization differently. In PyRadiomics,
discretization is applied once, based on either all values inside the ROI (when using a masked kernel,
``maskedKernel=True``) or on all values in the cropped image sent to the feature class (when using full kernel,
``maskedKernel=False``). ITK, or rather the ``itkHaralickTextureFeaturesImageFilter.cpp``, applies discretization on
each kernel separately.

Feature Calculation
===================
ITK feature calculation for GLCM features is defined in `itkScalarImageToTextureFeaturesFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx>`_.
However this class is more or less a wrapper that defines a `pipeline <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx#L125-L126>`_: first a glcm matrix is calculated,
which is then used to calculate feature values. The wrapper does contain the code to create the default offsets and
the code to `compute means <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx#L137-L181>`_ etc.
Below is a more detailed discussion on the separate steps of the calculation.

Matrix calculation
------------------
Defined in `itkScalarImageToCooccurrenceMatrixFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToCooccurrenceMatrixFilter.hxx>`_.
This class calculates the GLCM for an input image (and optional mask), for a certain offset (angle).
Furthermore, this also handles image discretization (see above).

Feature calculation
-------------------
Defined in `itkHistogramToTextureFeaturesFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkHistogramToTextureFeaturesFilter.hxx>`_.
This class takes the output of `itkScalarImageToCooccurrenceMatrixFilter`, an `itkHistogram`, and uses it to
compute the features. This histogram specifies 1 2D histogram, representing the GLCM matrix for a specific offset.

**Table1: Mapping ITK names to PyRadiomics' features**

========================= ==================
ITK name                  Pyradiomics Name
========================= ==================
energy                    `JointEnergy <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getJointEnergyFeatureValue>`_
entropy                   `JointEntropy <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getJointEntropyFeatureValue>`_
correlation               `Correlation <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getCorrelationFeatureValue>`_
inverseDifferenceMoment   `Idm <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getIdmFeatureValue>`_
inertia                   `Contrast <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getContrastFeatureValue>`_
clusterShade              `ClusterShape <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getClusterShadeFeatureValue>`_
clusterProminence         `ClusterProminence <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getClusterProminenceFeatureValue>`_
haralickCorrelation       `Correlation <https://pyradiomics.readthedocs.io/en/latest/features.html#radiomics.glcm.RadiomicsGLCM.getCorrelationFeatureValue>`_

========================= ==================

**N.B. correlation and haralickCorrelation are NOT part of the standard output of itkScalarImageToTextureFeaturesFilter**

Moreover, ``correlation`` and ``haralickCorrelation`` result in an identical value, differences found are caused by machine precision errors.


