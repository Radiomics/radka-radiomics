/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHaralickTextureFeaturesImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-04-06 00:19:17 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

  =========================================================================*/
#ifndef __itkHaralickTextureFeaturesImageFilter_cpp
#define __itkHaralickTextureFeaturesImageFilter_cpp

#include "itkHaralickTextureFeaturesImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstantBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkImageFileWriter.h"

#include "itkNumericTraits.h"
#include "itkVectorContainer.h"
#include "itkVector.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"

#include "itkMedianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"

#include <ctime>
#include <vcl_cmath.h>
#include <vcl_algorithm.h>

using namespace std;
namespace itk {

    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::HaralickTextureFeaturesImageFilter() {
        m_Radius.Fill(1);
        m_NumberOfBins = 16; //64;
        m_Mask = NULL;
        m_WindowSize = 5;

        // normalize parameters
        m_NormalizeTextures = true;
        m_ApplyMedianFiltering = true;

        m_NormalizeMax = 255.0;
        m_NormalizeMin = 0.0;

        //m_NormalizeMaxEnergy = 1.0;
        m_NormalizeMaxEnergy = 255.0;
        m_NormalizeMinEnergy = 0.0;

        //m_NormalizeMaxEntropy = 1.0;
        m_NormalizeMaxEntropy = 255.0;
        m_NormalizeMinEntropy = 0.0;

        //m_NormalizeMaxCorrelation = 2.0;
        m_NormalizeMaxCorrelation = 255.0;
        m_NormalizeMinCorrelation = 0.0;

        //m_NormalizeMaxDifferenceMoment = 1.0;
        m_NormalizeMaxDifferenceMoment = 255.0;
        m_NormalizeMinDifferenceMoment = 0.0;

        //m_NormalizeMaxInertia = 1.0;
        m_NormalizeMaxInertia = 255.0;
        m_NormalizeMinInertia = 0.0;

        m_StructureLabel = 1;
        m_AnalyzeRegion = false;

        m_DoSliceWiseAnalysis = false;

    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    void
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError) {
        std::cout << "GenerateInputRequestedRegion:" << std::endl;
        // call the superclass' implementation of this method
        Superclass::GenerateInputRequestedRegion();

        // get pointers to the input and output
        typename Superclass::InputImagePointer inputPtr = const_cast<TInputImage *>(this->GetInput());
        typename Superclass::OutputImagePointer outputPtr = this->GetOutput();

        if (!inputPtr || !outputPtr) {
            return;
        }

        // get a copy of the input requested region (should equal the output
        // requested region)
        typename TInputImage::RegionType inputRequestedRegion;
        inputRequestedRegion = inputPtr->GetRequestedRegion();
        cout << inputRequestedRegion << endl;

        std::cout << "Padding the input region by:" << m_Radius << std::endl;
        // pad the input requested region by the operator radius
        inputRequestedRegion.PadByRadius(m_Radius);

        std::cout << "Cropping the 'requested region' " << std::endl;
        // crop the input requested region at the input's largest possible region
        if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion())) {
            inputPtr->SetRequestedRegion(inputRequestedRegion);
            inputPtr->SetBufferedRegion(inputRequestedRegion);
            return;
        } else {

            std::cout << "Failed to crop! " << std::endl;
            // Couldn't crop the region (requested region is outside the largest
            // possible region).  Throw an exception.

            // store what we tried to request (prior to trying to crop)
            inputPtr->SetRequestedRegion(inputRequestedRegion);
            inputPtr->SetBufferedRegion(inputRequestedRegion);

            // build an exception
            InvalidRequestedRegionError e(__FILE__, __LINE__);
            e.SetLocation(ITK_LOCATION);
            e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
            e.SetDataObject(inputPtr);
            throw e;
        }
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetEnergyImage() {
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer energyImage = TInputImage::New();
        energyImage->CopyInformation(outputImage);
        energyImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        energyImage->Allocate();
        energyImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(energyImage, energyImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[0]);
        }
        return energyImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetEntropyImage() {
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer entropyImage = TInputImage::New();
        entropyImage->CopyInformation(outputImage);
        entropyImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        entropyImage->Allocate();
        entropyImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(entropyImage, entropyImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[1]);//HERE SEEMS TO BE THE SELECTION
        }
        return entropyImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetCorrelationImage() {
        std::cout << "GetCorrelationImage:" << std::endl;
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer correlationImage = TInputImage::New();
        correlationImage->CopyInformation(outputImage);
        correlationImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        correlationImage->Allocate();
        correlationImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(correlationImage, correlationImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[2]);//HERE SEEMS TO BE THE SELECTION
        }
        return correlationImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetDifferenceMomentImage() {
        std::cout << "GetDifferenceMomentImage:" << std::endl;
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer differenceMomentImage = TInputImage::New();
        differenceMomentImage->CopyInformation(outputImage);
        differenceMomentImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        differenceMomentImage->Allocate();
        differenceMomentImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(differenceMomentImage, differenceMomentImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[3]);//HERE SEEMS TO BE THE SELECTION
        }
        return differenceMomentImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::
    InputImagePointer HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetInertiaImage() {
        std::cout << "GetInertiaImage:" << std::endl;
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer inertiaImage = TInputImage::New();
        inertiaImage->CopyInformation(outputImage);
        inertiaImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        inertiaImage->Allocate();
        inertiaImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(inertiaImage, inertiaImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[4]);//HERE SEEMS TO BE THE SELECTION
        }
        return inertiaImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetClusterShadeImage() {
        std::cout << "GetClusterShadeImage:" << std::endl;
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer clusterShadeImage = TInputImage::New();
        clusterShadeImage->CopyInformation(outputImage);
        clusterShadeImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        clusterShadeImage->Allocate();
        clusterShadeImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(clusterShadeImage, clusterShadeImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[5]);//HERE SEEMS TO BE THE SELECTION
        }
        return clusterShadeImage;
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    const typename HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::GetClusterProminenceImage() {
        std::cout << "GetClusterProminenceImage:" << std::endl;
        typename TOutputImage::Pointer outputImage = const_cast<TOutputImage *> (this->GetOutput());
        typename TInputImage::Pointer clusterPromImage = TInputImage::New();
        clusterPromImage->CopyInformation(outputImage);
        clusterPromImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        clusterPromImage->Allocate();
        clusterPromImage->FillBuffer(0);
        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());
        ImageRegionIterator<TInputImage> en(clusterPromImage, clusterPromImage->GetBufferedRegion());
        for (out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en) {
            typename TOutputImage::PixelType vec = out.Get();
            en.Set(vec[6]);//HERE SEEMS TO BE THE SELECTION
        }
        return clusterPromImage;

    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    void
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::BeforeThreadedGenerateData() {
        std::cout << "BeforeThreadedGenerateData " << std::endl;
        // initialize the outputImage
        typename TOutputImage::Pointer output = this->GetOutput();
        if (output->GetBufferedRegion().GetNumberOfPixels() == 0) {
            output->SetBufferedRegion(this->GetInput()->GetBufferedRegion());
            output->Allocate();
            typename TOutputImage::PixelType pix;
            for (unsigned i = 0; i < pix.GetVectorDimension(); i++) {
                pix[i] = 0;
            }
            output->FillBuffer(pix);
        }

        m_Radius.Fill(m_WindowSize);
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    void
    HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::AfterThreadedGenerateData() {

        std::cout << "AfterThreadedGenerateData " << std::endl;
        typename TOutputImage::Pointer outputImage = this->GetOutput();//const_cast< TOutputImage*> (this->GetOutput());

        ImageRegionIterator<TOutputImage> out(outputImage, outputImage->GetBufferedRegion());

        // write out the texture images
        typedef double InternalPixelType;
        typedef Image<InternalPixelType, 3> InternalImageType;
        typename InternalImageType::Pointer inImage = InternalImageType::New();
        inImage->CopyInformation(outputImage);
        inImage->SetBufferedRegion(outputImage->GetBufferedRegion());
        inImage->Allocate();
        inImage->FillBuffer(0);

        typedef RescaleIntensityImageFilter<InternalImageType,
                InternalImageType> RescaleIntensityFilterType;


        // compute the texture statistics for the mask
        m_MeanEnergy = 0.0;
        m_MeanEntropy = 0.0;
        m_MeanCorrelation = 0.0;
        m_MeanDifferenceMoment = 0.0;
        m_MeanInertia = 0.0;
        m_MeanClusterShade = 0.0;
        m_MeanClusterProminence = 0.0;

        m_MedianEnergy = 0.0;
        m_MedianEntropy = 0.0;
        m_MedianCorrelation = 0.0;
        m_MedianDifferenceMoment = 0.0;
        m_MedianInertia = 0.0;
        m_MedianClusterShade = 0.0;
        m_MedianClusterProminence = 0.0;

        m_StdEnergy = 0.0;
        m_StdEntropy = 0.0;
        m_StdCorrelation = 0.0;
        m_StdDifferenceMoment = 0.0;
        m_StdInertia = 0.0;
        m_StdClusterShade = 0.0;
        m_StdClusterProminence = 0.0;

        m_KurtosisEnergy = 0.0;
        m_KurtosisEntropy = 0.0;
        m_KurtosisCorrelation = 0.0;
        m_KurtosisDifferenceMoment = 0.0;
        m_KurtosisInertia = 0.0;
        m_KurtosisClusterShade = 0.0;
        m_KurtosisClusterProminence = 0.0;

        m_SkewEnergy = 0.0;
        m_SkewEntropy = 0.0;
        m_SkewCorrelation = 0.0;
        m_SkewDifferenceMoment = 0.0;
        m_SkewInertia = 0.0;
        m_SkewClusterShade = 0.0;
        m_SkewClusterProminence = 0.0;

        for (unsigned n = 0; n < 7; ++n) {
            ImageRegionIterator<InternalImageType> in(inImage, inImage->GetBufferedRegion());
            for (in.GoToBegin(), out.GoToBegin(); !out.IsAtEnd(); ++in, ++out) {
                typename TOutputImage::PixelType vec = out.Get();
                in.Set(vec[n]);
            }

            if (m_NormalizeTextures) {
                std::cout << " Normalizing texture " << std::endl;

                // compute median filtering if that has been turned on
                if (m_ApplyMedianFiltering) {

                    typedef MedianImageFilter<InternalImageType, InternalImageType> MedianImageFilterType;
                    typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
                    medianFilter->SetInput(inImage);

                    typename TInputImage::SizeType radiusMedian;
                    radiusMedian.Fill(1);
                    radiusMedian[2] = 0;
                    std::cout << " Applying median filtering " << " Radius: " << radiusMedian << std::endl;
                    //https://itk.org/Doxygen/html/classitk_1_1MedianImageFilter.html0
                    medianFilter->SetRadius(radiusMedian);
                    medianFilter->Update();

                    // now rescale the intensity of the median filtered image
                    typename RescaleIntensityFilterType::Pointer rescaleFilter = RescaleIntensityFilterType::New();
                    rescaleFilter->SetInput(medianFilter->GetOutput());
                    std::cout << " Rescaling... " << std::endl;
                    switch (n) {
                        case 0: // Energy
                        {
                            std::cout << " Normalizing energy feature ... max " << m_NormalizeMaxEnergy << " min "
                                      << m_NormalizeMinEnergy << std::endl;
                            rescaleFilter->SetOutputMinimum(static_cast<unsigned long> (m_NormalizeMinEnergy));
                            rescaleFilter->SetOutputMaximum(static_cast<unsigned long> (m_NormalizeMaxEnergy));
                            break;
                        }
                        case 1: // Entropy
                        {
                            std::cout << " Normalizing entropy feature ... max " << m_NormalizeMaxEntropy << " min "
                                      << m_NormalizeMinEntropy << std::endl;
                            rescaleFilter->SetOutputMinimum(static_cast<unsigned long> (m_NormalizeMinEntropy));
                            rescaleFilter->SetOutputMaximum(static_cast<unsigned long> (m_NormalizeMaxEntropy));
                            break;
                        }
                        case 2: // Correlation
                        {
                            std::cout << " Normalizing correlation feature ... max " << m_NormalizeMaxCorrelation
                                      << " min " << m_NormalizeMinCorrelation << std::endl;
                            rescaleFilter->SetOutputMinimum(static_cast<unsigned long> (m_NormalizeMinCorrelation));
                            rescaleFilter->SetOutputMaximum(static_cast<unsigned long> (m_NormalizeMaxCorrelation));
                            break;
                        }
                        case 3: // Homogeneity
                        {
                            std::cout << " Normalizing contrast feature ... max " << m_NormalizeMaxDifferenceMoment
                                      << " min " << m_NormalizeMinDifferenceMoment << std::endl;
                            rescaleFilter->SetOutputMinimum(
                                    static_cast<unsigned long> (m_NormalizeMinDifferenceMoment));
                            rescaleFilter->SetOutputMaximum(
                                    static_cast<unsigned long> (m_NormalizeMaxDifferenceMoment));
                            break;
                        }
                        case 4: // Contrast
                        {
                            std::cout << " Normalizing homogeneity feature ... max " << m_NormalizeMaxInertia << " min "
                                      << m_NormalizeMinInertia << std::endl;
                            rescaleFilter->SetOutputMinimum(static_cast<unsigned long> (m_NormalizeMinInertia));
                            rescaleFilter->SetOutputMaximum(static_cast<unsigned long> (m_NormalizeMaxInertia));
                            break;
                        }
                        default: {
                            std::cout << " Normalizing features ... max " << m_NormalizeMax << " min " << m_NormalizeMin
                                      << std::endl;
                            rescaleFilter->SetOutputMinimum(static_cast<unsigned long> (m_NormalizeMin));
                            rescaleFilter->SetOutputMaximum(static_cast<unsigned long> (m_NormalizeMax));
                            break;
                        }
                    }
                    rescaleFilter->Update();

                    std::cout << " Done with rescaling " << n << std::endl;

                    // put this back into the output
                    ImageRegionIterator<InternalImageType> scalar(rescaleFilter->GetOutput(),
                                                                  rescaleFilter->GetOutput()->GetBufferedRegion());
                    ImageRegionIterator<TOutputImage> vector(outputImage,
                                                             outputImage->GetBufferedRegion());
                    for (scalar.GoToBegin(), vector.GoToBegin(); !scalar.IsAtEnd();
                         ++scalar, ++vector) {
                        typename TOutputImage::PixelType vec = vector.Get();
                        vec[n] = scalar.Get();
                        vector.Set(vec);
                    }
                } else {// No median filter
                    // simply rescale the intensity
                    typename RescaleIntensityFilterType::Pointer rescaleFilter = RescaleIntensityFilterType::New();
                    rescaleFilter->SetInput(inImage);
                    rescaleFilter->SetOutputMinimum(m_NormalizeMin);
                    rescaleFilter->SetOutputMaximum(m_NormalizeMax);
                    rescaleFilter->Update();

                    // put this back into the output
                    ImageRegionIterator<InternalImageType> scalar(rescaleFilter->GetOutput(),
                                                                  rescaleFilter->GetOutput()->GetBufferedRegion());
                    ImageRegionIterator<TOutputImage> vector(outputImage,
                                                             outputImage->GetBufferedRegion());
                    for (scalar.GoToBegin(), vector.GoToBegin(); !scalar.IsAtEnd();
                         ++scalar, ++vector) {
                        typename TOutputImage::PixelType vec = vector.Get();
                        vec[n] = scalar.Get();
                        vector.Set(vec);
                    }
                } // Else median filter
            } else {// If  Normalize and no median filter
                std::cout << "No Normalize and No MEDIAN FILTER " << std::endl;

                // put this back into the output
                ImageRegionIterator<InternalImageType> scalar(inImage, inImage->GetBufferedRegion());
                ImageRegionIterator<TOutputImage> vector(outputImage, outputImage->GetBufferedRegion());
                for (scalar.GoToBegin(), vector.GoToBegin(); !scalar.IsAtEnd();
                     ++scalar, ++vector) {
                    typename TOutputImage::PixelType vec = vector.Get();
                    vec[n] = scalar.Get();
                    vector.Set(vec);
                }
            }

            // Compute histogram features
            bool doComputation = false;
            if (doComputation) {
                out.GoToBegin();
                ImageRegionIterator<MaskImageType> mask(m_Mask, m_Mask->GetBufferedRegion());
                std::vector<InternalPixelType> vals;
                InternalPixelType meanval = 0.0;

                for (in.GoToBegin(), mask.GoToBegin(); !out.IsAtEnd(); ++in, ++out, ++mask) {

                    typename TOutputImage::PixelType vec = out.Get();
                    in.Set(vec[n]);
                    if (mask.Get() == m_StructureLabel) {
                        if (n == 2) // adjust to -1 to 1 for correlation
                        {
                            vals.push_back(in.Get() - 1.0);
                            meanval = meanval + (in.Get() - 1.0);
                        } else {
                            vals.push_back(in.Get());
                            meanval = meanval + in.Get();
                        }

                    }
                }

                std::sort(vals.begin(), vals.end());
                meanval /= (vals.size() > 0) ? (InternalPixelType) (vals.size()) : 1.0;

                InternalPixelType stdval = 0.0;
                InternalPixelType kurtosisVal = 0.0;
                InternalPixelType skewVal = 0.0;

                for (int i = 0; i < vals.size(); ++i) {
                    stdval = stdval + (vals[i] - meanval) * (vals[i] - meanval);
                    kurtosisVal = kurtosisVal + std::pow((vals[i] - meanval), 4.0);
                    skewVal = skewVal + std::pow((InternalPixelType) (vals[i] - meanval), (InternalPixelType) 3.0);
                }
                stdval /= (vals.size() > 1) ? (InternalPixelType) (vals.size() - 1) : 1.0;
                stdval = std::sqrt(stdval);

                kurtosisVal /= (std::max((InternalPixelType) 1.0, (InternalPixelType) (vals.size() - 1)) *
                                std::pow(stdval, 4.0));
                skewVal /= (std::max((InternalPixelType) 1.0, (InternalPixelType) (vals.size() - 1)) *
                            std::pow(stdval, 3.0));

                kurtosisVal -= 3;

                int mid = vals.size() / 2;
                switch (n) {
                    case 0: // Energy
                    {
                        m_MeanEnergy = meanval;
                        m_StdEnergy = stdval;
                        m_KurtosisEnergy = kurtosisVal;
                        m_SkewEnergy = skewVal;
                        m_MedianEnergy = vals[mid];
                        break;
                    }
                    case 1: // Entropy
                    {
                        m_MeanEntropy = meanval;
                        m_StdEntropy = stdval;
                        m_KurtosisEntropy = kurtosisVal;
                        m_SkewEntropy = skewVal;
                        m_MedianEntropy = vals[mid];
                        break;
                    }
                    case 2: // Correlation
                    {
                        m_MeanCorrelation = meanval;
                        m_StdCorrelation = stdval;
                        m_KurtosisCorrelation = kurtosisVal;
                        m_SkewCorrelation = skewVal;
                        m_MedianCorrelation = vals[mid];
                        break;
                    }
                    case 3: // Difference Moment
                    {
                        m_MeanDifferenceMoment = meanval;
                        m_StdDifferenceMoment = stdval;
                        m_KurtosisDifferenceMoment = kurtosisVal;
                        m_SkewDifferenceMoment = skewVal;
                        m_MedianDifferenceMoment = vals[mid];
                        break;
                    }
                    case 4: // Inertia
                    {
                        m_MeanInertia = meanval;
                        m_StdInertia = stdval;
                        m_KurtosisInertia = kurtosisVal;
                        m_SkewInertia = skewVal;
                        m_MedianInertia = vals[mid];
                        break;
                    }
                    case 5: // Cluster Shade
                    {
                        m_MeanClusterShade = meanval;
                        m_StdClusterShade = stdval;
                        m_KurtosisClusterShade = kurtosisVal;
                        m_SkewClusterShade = skewVal;
                        m_MedianClusterShade = vals[mid];
                        break;
                    }
                    case 6: {
                        // Cluster Prominence
                        m_MeanClusterProminence = meanval;
                        m_StdClusterProminence = stdval;
                        m_KurtosisClusterProminence = kurtosisVal;
                        m_SkewClusterProminence = skewVal;
                        m_MedianClusterProminence = vals[mid];
                        break;
                    }
                } // end switch

                std::cout << " Label " << m_StructureLabel << " Energy mean " << m_MeanEnergy << " kurtosis "
                          << m_KurtosisEnergy << " std " << m_StdEnergy
                          << std::endl;
                std::cout << " Entropy mean " << m_MeanEntropy << " kurtosis " << m_KurtosisEntropy << " std "
                          << m_StdEntropy << std::endl;
                std::cout << " Correlation mean " << m_MeanCorrelation << " kurtosis " << m_KurtosisCorrelation
                          << " std " << m_StdCorrelation << std::endl;
                std::cout << " Diff Moment mean " << m_MeanDifferenceMoment << " kurtosis "
                          << m_KurtosisDifferenceMoment << " std " << m_StdDifferenceMoment << std::endl;
                std::cout << " Inertia mean " << m_MeanInertia << " kurtosis " << m_KurtosisInertia << " std "
                          << m_StdInertia << std::endl;
                std::cout << " ClusterShade mean " << m_MeanClusterShade << " kurtosis " << m_KurtosisClusterShade
                          << " std " << m_StdClusterShade << std::endl;
                std::cout << " ClusterProminence mean " << m_MeanClusterProminence << " kurtosis "
                          << m_KurtosisClusterProminence << " std " <<
                          m_StdClusterProminence << std::endl;

            }
        }
    }


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    void HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>::GenerateData() {
        std::cout << " *** GenerateData ***" << std::endl;
        //defs in order to define "TextureFeaturesFilterType" to SetRequestedFeatures
        //for "itkScalarImageToTextureFeaturesFilter.h
        //in scope itk already
        typedef Statistics::ScalarImageToCooccurrenceMatrixFilter<TInputImage, Statistics::DenseFrequencyContainer2> CooccurrenceMatrixFilterType;
        typedef typename CooccurrenceMatrixFilterType::HistogramType HistType;
        typedef Statistics::HistogramToTextureFeaturesFilter<HistType> TextureFeaturesFilterType;
        typedef short TextureFeatureName;
        typedef VectorContainer<unsigned char, TextureFeatureName> FeatureNameVector;
        typedef typename FeatureNameVector::Pointer FeatureNameVectorPointer;

        //call the superclass to execute the threaded version for local texture computation
        if (!m_AnalyzeRegion) {
            Superclass::GenerateData();
            return;
        }

        // otherwise we execute the single threaded version and compute the texture for the entire mask
        typename TInputImage::Pointer inputImage = const_cast<TInputImage *>(this->GetInput());
        typedef Statistics::ScalarImageToTextureFeaturesFilter<TInputImage> TextureFilterType;

        InputPixelType minVal = NumericTraits<InputPixelType>::max(InputPixelType());
        InputPixelType maxVal = NumericTraits<InputPixelType>::min(InputPixelType()); //changed by Duc 11/11/2013
        ImageRegionIterator<TInputImage> it(inputImage, inputImage->GetBufferedRegion());
        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
            if (it.Get() < minVal) {
                minVal = it.Get();
            }
            if (it.Get() > maxVal) {
                maxVal = it.Get();
            }
        }
        typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
        typedef CastImageFilter<MaskImageType, TInputImage> CastImageFilterType;
        typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
        if (m_Mask.IsNotNull()) {
            castFilter->SetInput(m_Mask);
            castFilter->Update();

            textureFilter->SetInput(inputImage);
            textureFilter->SetPixelValueMinMax(minVal, maxVal);
            textureFilter->SetMaskImage(castFilter->GetOutput());
            textureFilter->SetNumberOfBinsPerAxis(m_NumberOfBins);
            //Duc 11/11/13: the value for the label needs to be set for labels that
            //are not 1, the default value for the inside of the mask and needs to be
            //somehow set in order to make this work
            //textureFilter->SetInsidePixelValue(6);
            //
            //this line should make it work
            textureFilter->SetInsidePixelValue(m_StructureLabel);

            //set the features we want
            FeatureNameVectorPointer requestedFeatures = FeatureNameVector::New();
            requestedFeatures->push_back(TextureFeaturesFilterType::Energy);
            requestedFeatures->push_back(TextureFeaturesFilterType::Entropy);
            requestedFeatures->push_back(TextureFeaturesFilterType::Correlation);
            requestedFeatures->push_back(TextureFeaturesFilterType::InverseDifferenceMoment);
            requestedFeatures->push_back(TextureFeaturesFilterType::Inertia);
            requestedFeatures->push_back(TextureFeaturesFilterType::ClusterShade);
            requestedFeatures->push_back(TextureFeaturesFilterType::ClusterProminence);
            textureFilter->SetRequestedFeatures(requestedFeatures);
        } else {
            //null mask
            textureFilter->SetInput(inputImage);
            textureFilter->SetPixelValueMinMax(minVal, maxVal);
            textureFilter->FastCalculationsOn();
            textureFilter->SetNumberOfBinsPerAxis(m_NumberOfBins);

            //don't forget this as well
            textureFilter->SetInsidePixelValue(m_StructureLabel);

            //don't forget to set the features we want
            FeatureNameVectorPointer requestedFeatures = FeatureNameVector::New();
            requestedFeatures->push_back(TextureFeaturesFilterType::Energy);
            requestedFeatures->push_back(TextureFeaturesFilterType::Entropy);
            requestedFeatures->push_back(TextureFeaturesFilterType::Correlation);
            requestedFeatures->push_back(TextureFeaturesFilterType::InverseDifferenceMoment);
            requestedFeatures->push_back(TextureFeaturesFilterType::Inertia);
            requestedFeatures->push_back(TextureFeaturesFilterType::ClusterShade);
            requestedFeatures->push_back(TextureFeaturesFilterType::ClusterProminence);
            textureFilter->SetRequestedFeatures(requestedFeatures);
        }

        try {
            textureFilter->Update();
        }
        catch (itk::ExceptionObject &etext) {
            std::cout << " error in texture filter " << etext << std::endl;
            return;
        }

        // now get the texture features and add them to the output
        typename TextureFilterType::FeatureValueVector::Pointer featureVector = TextureFilterType::FeatureValueVector::New();
        featureVector = textureFilter->GetFeatureMeans();
        unsigned nfeatures = featureVector->Size();
        for (unsigned j = 0; j < nfeatures; ++j) {
            switch (j) {
                case 0: {
                    m_MeanEnergy = featureVector->GetElement(j);
                    break;
                }
                case 1: {
                    m_MeanEntropy = featureVector->GetElement(j);
                    break;
                }
                case 2: {
                    m_MeanCorrelation = featureVector->GetElement(j);
                    break;
                }
                case 3: {
                    m_MeanDifferenceMoment = featureVector->GetElement(j);
                    break;
                }
                case 4: {
                    m_MeanInertia = featureVector->GetElement(j);
                    break;
                }
                case 5: {
                    m_MeanClusterShade = featureVector->GetElement(j);
                    break;
                }
                case 6: {
                    m_MeanClusterProminence = featureVector->GetElement(j);
                    break;
                }
            }
        }

        //added by Duc 01/01/2014
        //check if there are Nans in the output. If that happens, it shows that there is a problem with the textureanalyzer
        if (std::isnan(m_MeanEnergy) ||
            std::isnan(m_MeanEntropy) ||
            std::isnan(m_MeanCorrelation) ||
            std::isnan(m_MeanDifferenceMoment) ||
            std::isnan(m_MeanInertia) ||
            std::isnan(m_MeanClusterShade) ||
            std::isnan(m_MeanClusterProminence)) {
            std::cout << "THERE ARE SOME NAN VALUES" << std::endl;
            //find slices and see if that fixes the problem
            ImageRegionConstIteratorWithIndex<TInputImage> mit(castFilter->GetOutput(),
                                                               castFilter->GetOutput()->GetBufferedRegion());

            std::vector<int> xpos;
            std::vector<int> ypos;
            std::vector<int> zpos;
            //find the slices by finding the position of the label and finding the direction that remains constant
            for (mit.GoToBegin(); !mit.IsAtEnd(); ++mit) {
                bool found;
                if (mit.Get() == m_StructureLabel) {
                    found = false;
                    for (int i = 0; !found && i < xpos.size(); ++i) {
                        found = xpos[i] == mit.GetIndex()[0];
                    }
                    if (!found)
                        xpos.push_back(mit.GetIndex()[0]);

                    found = false;
                    for (int i = 0; !found && i < ypos.size(); ++i) {
                        found = ypos[i] == mit.GetIndex()[1];
                    }
                    if (!found)
                        ypos.push_back(mit.GetIndex()[1]);

                    found = false;
                    for (int i = 0; !found && i < zpos.size(); ++i) {
                        found = zpos[i] == mit.GetIndex()[2];
                    }
                    if (!found)
                        zpos.push_back(mit.GetIndex()[2]);

                }//if on m_structurelabel
            }
            m_SliceDir = 3;
            m_SliceNum;
            if (xpos.size() == 1) {
                m_SliceDir = 0;
                m_SliceNum = xpos[0];
            }
            if (ypos.size() == 1) {
                m_SliceDir = 1;
                m_SliceNum = ypos[0];

            }
            if (zpos.size() == 1) {
                m_SliceDir = 2;
                m_SliceNum = zpos[0];
            }

            //if m_SliceDir==3 then there are a couple of single slice pictures that are disconnected,
            //for now I am assuming that this does not happen. Damn fringe cases.... well I'll talk to Harini
            //about it and if she thinks this is essential I'll change it. let's get the single slice to work first
            //An initial idea for testing is to do something similar to the diff function in matlab. if diff is not 1
            //from element to element, that means there are jumps in slice number. Lets get the other stuff the single slice to work first
            if (m_SliceDir != 3) {
                std::cout << " HERE and value is: " << m_SliceDir << std::endl;
                //extract the slice image
                typedef itk::Image<float, 2> SliceType;
                typedef SliceType::Pointer SlicePointerType;

                typedef ExtractImageFilter<TInputImage, SliceType> ExFilterType;
                typename ExFilterType::Pointer exFilt = ExFilterType::New();
                exFilt->SetDirectionCollapseToIdentity();

                typename InputImageType::RegionType exRegion;
                typename InputImageType::IndexType exStart;
                typename InputImageType::SizeType exSize;

                exRegion = inputImage->GetBufferedRegion();
                exSize = exRegion.GetSize();
                //std::cout << iRegion << iSize << std::endl;
                exSize[m_SliceDir] = 0;

                exStart = exRegion.GetIndex();
                exStart[m_SliceDir] = m_SliceNum;

                //build region that extrtacts slice
                typename InputImageType::RegionType desiredRegion;
                desiredRegion.SetSize(exSize);
                desiredRegion.SetIndex(exStart);

                //extract slice
                exFilt->SetExtractionRegion(desiredRegion);
                exFilt->SetInput(inputImage);
                exFilt->Update();

                //extract the same slice from the labelmask
                typename ExFilterType::Pointer exLabel = ExFilterType::New();
                exLabel->SetDirectionCollapseToIdentity();
                exLabel->SetExtractionRegion(desiredRegion); //same as for the input image
                exLabel->SetInput(
                        castFilter->GetOutput()); //the cast makes it possible to use the same extarctfiltertype
                exLabel->Update();

                //build filter
                float minVal = NumericTraits<float>::max(float());
                float maxVal = NumericTraits<float>::min(
                        float()); //added by Duc, saying 0 takes care of the negative values....
                //float maxVal = -1 * NumericTraits< float > ::max(float());
                ImageRegionIterator<SliceType> sit(exFilt->GetOutput(), exFilt->GetOutput()->GetBufferedRegion());
                for (sit.GoToBegin(); !sit.IsAtEnd(); ++sit) {
                    if (sit.Get() < minVal) {
                        minVal = sit.Get();
                    }
                    if (sit.Get() > maxVal) {
                        maxVal = sit.Get();
                    }
                }

                typedef Statistics::ScalarImageToTextureFeaturesFilter<SliceType> TextureFilter2DType;
                typename TextureFilter2DType::Pointer textureFilter2D = TextureFilter2DType::New();
                textureFilter2D->SetInput(exFilt->GetOutput());
                textureFilter2D->SetPixelValueMinMax(minVal, maxVal);
                textureFilter2D->SetMaskImage(exLabel->GetOutput());
                textureFilter2D->SetNumberOfBinsPerAxis(m_NumberOfBins);
                textureFilter2D->SetInsidePixelValue(m_StructureLabel);

                //set the features we want
                FeatureNameVectorPointer requestedFeatures2D = FeatureNameVector::New();
                requestedFeatures2D->push_back(TextureFeaturesFilterType::Energy);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::Entropy);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::Correlation);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::InverseDifferenceMoment);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::Inertia);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::ClusterShade);
                requestedFeatures2D->push_back(TextureFeaturesFilterType::ClusterProminence);
                textureFilter2D->SetRequestedFeatures(requestedFeatures2D);

                try {
                    textureFilter2D->Update();
                }
                catch (itk::ExceptionObject &etext) {
                    std::cout << " error in 2D texture filter " << etext << std::endl;
                    return;
                }

                // now get the texture features and add them to the output
                typename TextureFilterType::FeatureValueVector::Pointer featureVector2D = TextureFilterType::FeatureValueVector::New();
                featureVector2D = textureFilter2D->GetFeatureMeans();
                unsigned nfeatures = featureVector2D->Size();
                for (unsigned j = 0; j < nfeatures; ++j) {
                    switch (j) {
                        case 0: {
                            m_MeanEnergy = featureVector2D->GetElement(j);
                            break;
                        }
                        case 1: {
                            m_MeanEntropy = featureVector2D->GetElement(j);
                            break;
                        }
                        case 2: {
                            m_MeanCorrelation = featureVector2D->GetElement(j);
                            break;
                        }
                        case 3: {
                            m_MeanDifferenceMoment = featureVector2D->GetElement(j);
                            break;
                        }
                        case 4: {
                            m_MeanInertia = featureVector2D->GetElement(j);
                            break;
                        }
                        case 5: {
                            m_MeanClusterShade = featureVector2D->GetElement(j);
                            break;
                        }
                        case 6: {
                            m_MeanClusterProminence = featureVector2D->GetElement(j);
                            break;
                        }
                    }//switch
                }//for nfeatures
            }//if slicedir ~=3
        }//isnan m_Mean....

        std::cout << " ---- out GenerateData --- " << std::endl;
    }//generatedata


    template<class TInputImage, class TOutputImage, class TMaskPixelType>
    void HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
    ::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType itkNotUsed(threadId))
    //::GenerateData()
    {

        cout << "=======  Region size: "<< outputRegionForThread.GetSize() << endl;
//        std::cout << " *** ThreadedGenerateData***" << std::endl;
        //  Initialize();
        typename TInputImage::Pointer inputImage = const_cast<TInputImage *>(this->GetInput());

        InputPixelType minVal = NumericTraits<InputPixelType>::max(InputPixelType());
        InputPixelType maxVal = NumericTraits<InputPixelType>::min(InputPixelType());
        ImageRegionIterator<TInputImage> it(inputImage, outputRegionForThread); //inputImage->GetBufferedRegion()); //outputRegionForThread);
        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
            if (it.Get() < minVal) {
                minVal = it.Get();
            }
            if (it.Get() > maxVal) {
                maxVal = it.Get();
            }
        }
        ConstantBoundaryCondition<MaskImageType> nbc;
        nbc.SetConstant(0);
        ConstNeighborhoodIterator<MaskImageType> bit;
        NeighborhoodIterator<TOutputImage> mit;
        typename OutputImageType::Pointer output = this->GetOutput();

        // Find the data-set boundary "faces"
        typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
        NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
        faceList = bC(inputImage, outputRegionForThread,
                      m_Radius); //inputImage->GetBufferedRegion(), m_Radius); //outputRegionForThread, m_Radius);
        typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
        typedef Statistics::ScalarImageToTextureFeaturesFilter<TInputImage> TextureFilterType;
        typedef Statistics::ScalarImageToTextureFeaturesFilter<Image<InputPixelType, 2> > TextureFilter2DType;

        // support progress methods/callbacks
        //ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

        typename TInputImage::SizeType imsize = inputImage->GetBufferedRegion().GetSize();
        //typename TInputImage::SizeType imsize = outputRegionForThread.GetSize();

        typename InputImageType::RegionType roiRegion;
        typename InputImageType::IndexType roiStart;
        typename InputImageType::SizeType roiSize;
        typename InputImageType::SpacingType roiSpacing = inputImage->GetSpacing();
        typename InputImageType::PointType roiOrigin = inputImage->GetOrigin();

        // Process each of the boundary faces.  These are N-d regions which border
        // the edge of the buffer.
        typedef ExtractImageFilter<InputImageType, InputImageType> ExtractFilterType;

        //added by Duc
        typedef ExtractImageFilter<InputImageType, Image<InputPixelType, 2> > ExtractFilter2DType;
//        cout << "Number of faces: " << faceList.size() << endl;

        for (fit = faceList.begin(); fit != faceList.end(); ++fit) {
            if (m_Mask.IsNotNull()) { //Mask not null
                if (m_AnalyzeRegion) {
                    continue;
                }
                bit = ConstNeighborhoodIterator<MaskImageType>(m_Radius, m_Mask, *fit);
                unsigned int neighborhoodSize = bit.Size();
//                bit.OverrideBoundaryCondition(&nbc);
                bit.GoToBegin();
                mit = NeighborhoodIterator<OutputImageType>(m_Radius, output, *fit);
                mit.GoToBegin();
//                cout << "\n Neighborhood size: " << bit.Size();
//                cout << " Iterator dimension: " << bit.Dimension << endl;
//                cout << " Radius: " << bit.GetRadius(0)<< ","<< bit.GetRadius(1)<< ","<< bit.GetRadius(2)<< endl;
//                cout << " Size (rect): " << bit.GetSize(0)<< ","<< bit.GetSize(1)<< ","<< bit.GetSize(2)<< endl;
                while (!bit.IsAtEnd()) {
                    if (bit.GetCenterPixel()) {
                        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
                        extractFilter->SetInput(inputImage);
                        typename MaskImageType::IndexType first = bit.GetIndex(0);
                        typename MaskImageType::IndexType last = bit.GetIndex(neighborhoodSize - 1);
//                        cout << "ImgSize: " << imsize << " First: " << first << "   Last: " << last << endl;
                        for (unsigned n = 0; n < inputImage->GetImageDimension(); n++) {
                            roiStart[n] = first[n] >= 0 ? first[n] : 0;
                            roiSize[n] = (last[n] < imsize[n]) ? (last[n] - roiStart[n]) : (imsize[n] - 1 - roiStart[n]);
//                            std::cout << "(" << last[n] << "<" << imsize[n] << ")? (" << last[n] << "-" << roiStart[n] << "):(" << imsize[n] << "-1-" << roiStart[n] << ")" << std::endl;
                        }
//                        std::cout << "Final Size: " << roiSize << " Roi Start: " << roiStart << " Pixel index:" << bit.GetIndex() << std::endl<< std::endl;
                        if (!m_DoSliceWiseAnalysis) //updated by Duc 02/18/2016
                        {
//							std::cout<< "@@@@@@@@@@@@ DOSliceWiseAnalysis???" << std::endl;
                            roiRegion.SetSize(roiSize);
                            roiRegion.SetIndex(roiStart);
                            extractFilter->SetExtractionRegion(roiRegion);
//                            std::cout << "Roi Size: " << roiSize << " Roi Start: " << roiStart << std::endl;
                            try {
                                extractFilter->Update();
                            }
                            catch (itk::ExceptionObject &eExtract) {
                                std::cout << " Error in extract object " << roiStart << " " << roiSize << " "
                                          << eExtract << std::endl;
                                typename TOutputImage::PixelType vPix;
                                for (unsigned j = 0; j < vPix.GetVectorDimension(); j++) {
                                    vPix[j] = 0;
                                }
                                mit.SetCenterPixel(vPix);
                                ++bit;
                                ++mit;
                                //progress.CompletedPixel();
                                continue;
                            }
                            InputPixelType minVal = NumericTraits<InputPixelType>::max(InputPixelType());
                            InputPixelType maxVal = NumericTraits<InputPixelType>::min(InputPixelType());
                            ImageRegionIterator<TInputImage> it(extractFilter->GetOutput(),
                                                                extractFilter->GetOutput()->GetBufferedRegion()); //outputRegionForThread);
//							std::cout<< "@@@@@@@@@@@@ 1" << std::endl;
                            int count = 0;
                            for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
                                if (it.Get() < minVal) {
                                    minVal = it.Get();
                                }
                                if (it.Get() > maxVal) {
                                    maxVal = it.Get();
                                }
                                count++;
                            }
//                            std::cout << "Number of iterations: " << count <<  "  should be around 5x5x5" << std::endl;
//                            std::cout << "Max: " << maxVal<<  "  Min:" << minVal << std::endl;
                            typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
                            textureFilter->SetInput(extractFilter->GetOutput());
                            textureFilter->SetPixelValueMinMax(minVal, maxVal);
                            textureFilter->SetNumberOfBinsPerAxis(m_NumberOfBins);
                            try {
//								std::cout<< "@@@@@@@@@@@@ 2" << std::endl;
                                textureFilter->Update();
                            }
                            catch (itk::ExceptionObject &etext) {
                                std::cout << " error in texture filter " << etext << std::endl;
                                return;
                            }
                            // now get the texture features and add them to the output
                            typename TextureFilterType::FeatureValueVector::Pointer featureVector = TextureFilterType::FeatureValueVector::New();
                            featureVector = textureFilter->GetFeatureMeans();
                            typename TOutputImage::PixelType vPix;
                            for (unsigned j = 0; j < vPix.GetVectorDimension(); j++) {
                                vPix[j] = 0;
                            }
                            unsigned nfeatures =
                                    vPix.GetVectorDimension() < featureVector->Size() ? vPix.GetVectorDimension() :
                                    featureVector->Size();
                            for (unsigned j = 0; j < nfeatures; j++) {
                                vPix[j] = featureVector->GetElement(j);
                                //std::cout<<" texture : "<<vPix[j]<<" ";
                            }
                            //std::cout<<""<<std::endl;
                            mit.SetCenterPixel(vPix);
                        } //if on DoSLiceWiseAnazlysis
                        else {
                            std::cout << "ERROR SHOULDN'T BE HERE!" << std::endl;
                        } //else on DoSLiceWiseAnazlysis
                    }//if on bit
                    ++bit;
                    ++mit;
                    //progress.CompletedPixel();
                } // end while
            } else//MASK IS NULL
            {
                std::cout << "ERROR SHOULDN'T BE HERE! MASK NULL" << std::endl;
            }
        }
    }

    /**
     * Standard "PrintSelf" method
     */
    template<class TInputImage, class TOutput, class TMaskPixelType>
    void HaralickTextureFeaturesImageFilter<TInputImage, TOutput, TMaskPixelType>
    ::PrintSelf(
            std::ostream &os,
            Indent indent) const {
        Superclass::PrintSelf(os, indent);
        os << indent << "Radius: " << m_Radius << std::endl;
        os << indent << "Number of Bins: " << m_NumberOfBins << std::endl;
    }

} // end namespace itk

#endif
