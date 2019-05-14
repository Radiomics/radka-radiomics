#############################################################
Comparison of ITK Feature calculation to PyRadiomics (v2.1.2)
#############################################################

This readme details the comparison of segment-based features calculated by ITK,
versus features calculated by PyRadiomics.

GLCM
****

ITK feature calculation for GLCM features is defined in `itkScalarImageToTextureFeaturesFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx>`_.
However this class is more or less a wrapper that defines a `pipeline <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx#L125-L126>`_: first a glcm matrix is calculated,
which is then used to calculate feature values. The wrapper does contain the code to create the default offsets and
the code to `compute means <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToTextureFeaturesFilter.hxx#L137-L181>`_ etc.
Below is a more detailed discussion on the separate steps of the calculation.

Pre-Processing
==============

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

Matrix calculation
==================
Defined in `itkScalarImageToCooccurrenceMatrixFilter <https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Numerics/Statistics/include/itkScalarImageToCooccurrenceMatrixFilter.hxx>`_.
This class calculates the GLCM for an input image (and optional mask), for a certain offset (angle).
Furthermore, this also handles image discretization (see above).

Feature calculation
===================
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

Comparison of calculated values
-------------------------------

For the comparison of features calculated in ITK and in PyRadiomics, the GLCM matrix as calculated by PyRadiomics was used
for both ITK feature calculation (using ``itkHistogramToTextureFeaturesFilter``) and PyRadiomics feature calculation.

Most features yielded similar results when calculated by PyRadiomics and ITK, with exception of ``correlation`` and
``haralickCorrelation`` (Table shows results for PyRadiomics' testcase ``brain1``)

=================== ============= ===================== ==============
Feature             ITK           PyRadiomics           PyRadiomics (alternative formulas)
=================== ============= ===================== ==============
ClusterProminence     1778.634857    1778.6348565874569
ClusterShade             4.514500    4.5144995523700375
Contrast                11.976152     11.97615161066193
Correlation              0.039620   0.39126009388396177       0.039620
HaralickCorrelation  35488.036648                         35488.036648
Idm                      0.333354   0.33335376144309525
JointEnergy              0.010592  0.010592150839272757
JointEntropy             6.923483     6.923482939238772
=================== ============= ===================== ==============

``correlation``: The difference in calculation stems from the fact that in ITK, the square of pixel variances is used,
whereas in PyRadiomics, the square of standard deviations is used (:math:`sigx * sigy`, which is equal to
:math:`sigx^2` for symmetrical GLCM). Furthermore, the ITK documentation also shows the formula using
:math:`\sigma^2`, rather than :math:`var^2`.

``haralickCorrelation``: In PyRadiomics ``haralickCorrelation`` is not implemented, as the formula would just yield
comparable results to ``correlation`` (equivalent formulas). When comparing to PyRadiomics' ``correlation`` feature,
differences still exist however. These are caused by the fact that the mean and standard deviations used in ITK's
implementation do not take the gray value into account (i.e. average and standard deviation of the gray values in
the ROI), but rather just the mean and standard deviation of the probabilities. The literal interpretation of the
feature as it is proposed in the original Haralick paper reflects the implementation in ITK ("where :math:`\mu_x`,
:math:`\mu_y`, :math:`\sigma_x` and :math:`\sigma_y` are the means and standard deviations of :math:`p_x` and
:math:`p_y`" [1]_), but seems contrary to an abstract definition of correlation (see related wiki article).
Finally, Haralick et al. does not provide equations for the calculation of the means and standard deviations used.

The code used for these comparisons can be found in
`itk_glcm_feature_comparison.py <https://github.com/Radiomics/UM_WIth_Olmo_Comments/tree/master/pyradiomics/itk_glcm_feature_comparison.py>`_

Related wikipedia article: `Correlation and Dependence <https://en.wikipedia.org/wiki/Correlation_and_dependence>`_.

.. [1] Haralick, R., Shanmugan, K., Dinstein, I; Textural features for image classification;
    IEEE Transactions on Systems, Man and Cybernetics; 1973(3), p610-621
