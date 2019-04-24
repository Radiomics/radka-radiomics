/**
 * Compute Texture Features in Fine, Medium, and Coarse settings for a batch of files
 *
 * Author: Harini Veeraraghavan
 * veerarah at mskcc dot org
 * Date: May 29nd 2013
 **/

#include "itkImage.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include <itkImageSeriesReader.h>
#include "itkImageSeriesWriter.h"
#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskImageFilter.h"

#include "itkCastImageFilter.h"

#include "itkImageSeriesReader.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include <itkImageRegionIteratorWithIndex.h>

#include "itkHaralickTextureFeaturesImageFilter.h"
#include "itkHaralickTextureFeaturesImageFilter.cpp"

#include "itkVector.h"

#include<fstream>
#include<string>

#include <stdlib.h>

const unsigned int Dimension = 3;

typedef short LabelPixelType;
typedef float InPixelType;

typedef itk::Image<InPixelType, Dimension> InImageType;
typedef itk::Image<LabelPixelType, Dimension> MaskImageType;

typedef InImageType::Pointer InImagePointerType;

typedef itk::Vector< double, 7> TexturePixelType;
typedef itk::Image< TexturePixelType, Dimension> TextureImageType;

std::string m_patientName = "";
std::string m_imageType = "";
std::string m_patientID = "";
std::string m_studyDate = "";
std::string m_studyUID = "";

bool ReadDICOMImageSeries(std::string dirname, InImagePointerType inputImage)
{
	std::cout << "Reading directory " << dirname << std::endl;

	typedef itk::ImageSeriesReader< InImageType > ReaderType;
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;

	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetUseSeriesDetails(false);
	nameGenerator->SetLoadSequences(true);
	nameGenerator->SetLoadPrivateTags(true);
	nameGenerator->SetInputDirectory(dirname.c_str());
	std::vector< std::string > filenames;

	// Get filenames
	try
	{
		const ReaderType::FileNamesContainer &fileNamesInSeries = nameGenerator->GetInputFileNames();
		if (fileNamesInSeries.size() > 1)
		{
			std::size_t nFiles = fileNamesInSeries.size();
			for (std::size_t k = 0; k < nFiles; k++)
			{
				filenames.push_back(fileNamesInSeries[k]);
			}
		}
		else
		{
			std::cout << "	error in getting the input dcm files - No files found in " << dirname.c_str() << std::endl;
			return false;
		}
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "	error in getting the input dcm files " << err.GetDescription() << std::endl;
		return false;
	}

	// Create input image
	ReaderType::Pointer reader = ReaderType::New();
	itk::GDCMImageIO::Pointer gdcmIO = itk::GDCMImageIO::New();
	reader->SetImageIO(gdcmIO);
	reader->SetFileNames(filenames);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "	error in reading the dicom series " << err.GetDescription() << std::endl;
		return false;
	}
	// Save DICOM infos
	itk::MetaDataDictionary &dict = gdcmIO->GetMetaDataDictionary();
	itk::ExposeMetaData<std::string>(dict, "0010|0010", m_patientName); // Patient name
	itk::ExposeMetaData<std::string>(dict, "0010|0020", m_patientID); // Patient id
	itk::ExposeMetaData<std::string>(dict, "0008|0020", m_studyDate); // Study date
	itk::ExposeMetaData<std::string>(dict, "0020|000d", m_studyUID); // Study instance UID

	inputImage->Graft(reader->GetOutput());

	return true;
}

int toDicom(std::string mhdFilePath, std::string outPath, std::string patientName, std::string patientID, std::string studyDate, std::string seriesDescription, std::string studyUID, std::string studyDescription)
{
	std::cout << "Inside to Dicom..." << std::endl;
	typedef signed short    PixelType;
	const unsigned int      Dimension = 3;
	typedef itk::Image< PixelType, Dimension >      ImageType;
	typedef itk::ImageFileReader< ImageType >       ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(mhdFilePath);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::GDCMImageIO                        ImageIOType;
	typedef itk::NumericSeriesFileNames             NamesGeneratorType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	itksys::SystemTools::MakeDirectory(outPath);
	typedef signed short    OutputPixelType;
	const unsigned int      OutputDimension = 2;
	typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
	typedef itk::ImageSeriesWriter<ImageType, Image2DType >  SeriesWriterType;

	//General infos
	itk::MetaDataDictionary &dict = gdcmIO->GetMetaDataDictionary();
	itk::EncapsulateMetaData<std::string>(dict, "0008|0060", "MR"); // Modality

	//Patient infos
	itk::EncapsulateMetaData<std::string>(dict, "0010|0010", patientName); // Patient name
	itk::EncapsulateMetaData<std::string>(dict, "0010|0020", patientID); // Patient id
	itk::EncapsulateMetaData<std::string>(dict, "0008|0020", studyDate); // Study date
	itk::EncapsulateMetaData<std::string>(dict, "0008|103e", seriesDescription); // Series description
	std::string seriesID = patientName + "." + patientID + "." + seriesDescription + "." + studyUID;
	itk::EncapsulateMetaData<std::string>(dict, "0020|000e", seriesID); // Series instance id
	itk::EncapsulateMetaData<std::string>(dict, "0020|0010", patientID); // Study id same as patient ID
	itk::EncapsulateMetaData<std::string>(dict, "0008|1030", studyDescription); // Study description
	itk::EncapsulateMetaData<std::string>(dict, "0020|000d", studyUID); // Study instance UID
	gdcmIO->KeepOriginalUIDOn();
	std::cout << patientName << " " << seriesID << " " << seriesDescription << " " << patientID << " " << studyDescription << " " << studyUID << std::endl;

	SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	seriesWriter->SetInput(reader->GetOutput());
	seriesWriter->SetImageIO(gdcmIO);
	ImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
	ImageType::IndexType start = region.GetIndex();
	ImageType::SizeType  size = region.GetSize();

	std::string format = outPath;
	format += "/IMG%04d.dcm";
	NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
	namesGenerator->SetSeriesFormat(format.c_str());
	namesGenerator->SetStartIndex(start[2]);
	namesGenerator->SetEndIndex(start[2] + size[2] - 1);
	namesGenerator->SetIncrementIndex(1);

	seriesWriter->SetFileNames(namesGenerator->GetFileNames());
	try
	{
		seriesWriter->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Exception thrown while writing the series " << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{

	char* values[6];
	values[1] = "/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/PX0000/ProstateX-0000_ProstateX-0000_MR_2011-07-07_114731_MR.prostaat.kanker.detectie.WDS.mc.MCAPRODETW_t2.tse.tra_n19__00000";
	values[2] = "/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/PX0000/ProstateX-0000_ProstateX-0000_MR_2011-07-07_114731_MR.prostaat.kanker.detectie.WDS.mc.MCAPRODETW_T2.Mask_n19__00000";
	values[3] = "/media/osz1/DATA/DATA/RADIOMICS_and_MRI_Normalization/TestDataForCppITK/OUTPUT/T2";
	values[4] = "2";
	values[5] = "16";
	argv = values;
	argc = 6;

	// Get arguments
	//--------------
	// Error message for missing arguments
	if (argc != 6 && argc != 9)
	{
		std::cout << "Usage <mhd File Path> <dest Folder> <patientName> <patientID> "
						"<studyDate> <seriesDescription> <studyUID> <studyDescription>" << std::endl;
		std::cout << "OR" << std::endl;
		std::cout << "Usage <source Folder> <mask Folder> <dest Folder> <window Setting> <n Bins>" << std::endl;
		std::cout << "Number of arguments: " << argc << std::endl;
		std::cout << "Arguments given:" << std::endl;
		for (int iArg = 0; iArg < argc; ++iArg)
		{
			std::cout << argv[iArg] << std::endl;
		}
		return -1;
	}
	// Dicom treatment for  mhd input
	if (argc == 9)
	{
		std::string sourceName = argv[1];
		std::string destFolder = argv[2];
		std::string patientName = argv[3];
		std::string patientID = argv[4];
		std::string studyDate = argv[5];
		std::string seriesDescription = argv[6];
		std::string studyUID = argv[7];
		std::string studyDescription = argv[8];
		return toDicom(sourceName, destFolder, patientName, patientID, studyDate, seriesDescription, studyUID, studyDescription);
	}

	std::cout << "Gooing with the 6 parameter options, this means RADIOMICS!!!!" << std::endl;
	// Treatment of arguments
	//------------------------
	// Image type
	std::string sourceName = argv[1];
	if (sourceName.find("Apparent") != std::string::npos || sourceName.find("ADC") != std::string::npos) {
		m_imageType = "ADC";
	}
	else if (sourceName.find("T2") != std::string::npos || sourceName.find("t2") != std::string::npos) {
		m_imageType = "T2";
	}
	else if (sourceName.find("k.trans") != std::string::npos || sourceName.find("k_trans") != std::string::npos) {
		m_imageType = "k_trans";
	}
	else if (sourceName.find("BVAL") != std::string::npos || sourceName.find("b-val") != std::string::npos) {
		m_imageType = "BVAL";
	}
	std::cout << "Working with: " << m_imageType<< std::endl;
	//Mask folder
	std::string maskName = argv[2];
	// Destination folder
	std::string destFolder = argv[3];
	// Window size
//	int windowSizeInput = atoi(argv[4]);
	int windowSizeInput = 2;
	if (windowSizeInput != 1 && windowSizeInput != 2 && windowSizeInput != 3)
	{
		std::cout << "Window setting : 1 -> 3x3, 2 -> 5x5, 3 -> 7x7" << std::endl;
		return -1;
	}else{
		std::cout << "Window size 1 -> 3x3, 2 -> 5x5, 3 -> 7x7 using:"<< windowSizeInput << std::endl;
	}
	// Number of bins
	std::string nBinsStr = argv[5];
	int nBins = atoi(argv[5]);
	if (nBins != 8 && nBins != 16 && nBins != 32 && nBins != 64 && nBins != 128)
	{
		std::cout << "Number of bins should be equal to : 8, 16, 32, 64, or 128" << std::endl;
		return -1;
	}else{
		std::cout << "Number of bins:"<< nBins << std::endl;
	}

	// Treatment of the input image
	//------------------------------
	// Read in the source image
	InImageType::Pointer sourceImage = InImageType::New();
	typedef itk::ImageFileReader<InImageType> ImageReaderType;
	std::cout << "Reading source images " << std::endl;
	try
	{
		if (ReadDICOMImageSeries(sourceName, sourceImage) == false)
		{
			return -1;
		}
		std::cout << "	Done reading source series " << std::endl;
	}
	catch (itk::ExceptionObject &ex1)
	{
		std::cout << "	Error reading source series " << sourceName << "  " << ex1.GetDescription() << std::endl;
		return -1;
	}
	// Read in the mask image
	MaskImageType::Pointer maskImage = MaskImageType::New();
	typedef itk::ImageFileReader<MaskImageType> MaskImageReaderType;
	if (maskName == "")
	{
		std::cout << "Creating mask image " << std::endl;
		maskImage->CopyInformation(sourceImage);
		maskImage->SetBufferedRegion(sourceImage->GetBufferedRegion());
		maskImage->Allocate();
		maskImage->FillBuffer(1);
	}
	else
	{
		std::cout << "Reading mask image " << std::endl;
		try
		{
			InImageType::Pointer maskInputImage = InImageType::New();
			typedef itk::ImageFileReader<InImageType> ImageReaderType;
			if (ReadDICOMImageSeries(maskName, maskInputImage) == false)
			{
				return -1;
			}
			std::cout << "	Done reading mask series " << std::endl;
			maskImage->CopyInformation(maskInputImage);
			maskImage->SetBufferedRegion(maskInputImage->GetBufferedRegion());
			maskImage->Allocate();
			maskImage->FillBuffer(0);
			// Populate the mask
			typedef itk::ImageRegionIteratorWithIndex<MaskImageType> IteratorType;
			IteratorType itr(maskImage, maskImage->GetBufferedRegion());
			typedef itk::ImageRegionIteratorWithIndex<InImageType> IteratorImageType;
			IteratorImageType itrImage(maskInputImage, maskInputImage->GetBufferedRegion());
			int max1 = 0;
			int max2 = 0;
			int max3 = 0;
			for (itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
			{
				MaskImageType::IndexType index = itr.GetIndex();
				if (index[0] > max1)
				{
					max1 = index[0];
				}
				if (index[1] > max2)
				{
					max2 = index[1];
				}
				if (index[2] > max3)
				{
					max3 = index[2];
				}
			}
			for (itr.GoToBegin(); !itr.IsAtEnd(); ++itr, ++itrImage)
			{
				MaskImageType::IndexType index = itr.GetIndex();
				if (index[0] > 0 && index[1] > 0 && index[2] > 0 &&
					index[0] < max1 && index[1] < max2 && index[2] < max3)
				{
					int value = itrImage.Get();
					if (value == 1)
					{
						itr.Set(1);
					}
				}
			}
		}
		catch (itk::ExceptionObject &ex1)
		{
			std::cout << "	error reading mask series " << maskName << "  " << ex1.GetDescription() << std::endl;
			return -1;
		}
	}
	// Get the number of labels in the image
	std::vector< LabelPixelType > labels;
	itk::ImageRegionIterator<MaskImageType> label(maskImage, maskImage->GetBufferedRegion());
	for (label.GoToBegin(); !label.IsAtEnd(); ++label)
	{
		if (label.Get() > 0)
		{
			bool found = false;
			for (unsigned i = 0; !found && i < labels.size(); ++i)
			{
				found = labels[i] == label.Get();
			}
			if (!found)
			{
				labels.push_back(label.Get());
				std::cout << "Label found: " << label.Get() <<std::endl;
			}
		}
	}

	// Prepare computation
	//---------------------
	typedef itk::HaralickTextureFeaturesImageFilter<InImageType, TextureImageType, short> HaralickTextureFeaturesFilterType;
	for (unsigned n = 0; n < labels.size(); ++n)
	{
		// Set up the call to the Haralick Texture Feature Filter
		HaralickTextureFeaturesFilterType::Pointer hFilter = HaralickTextureFeaturesFilterType::New();
		hFilter->SetInput(sourceImage);
		hFilter->SetMask(maskImage);
		std::string windowSize = "";
		switch (windowSizeInput) {
			case 1: {
			    cout << "Setting window size 3x3x3" << endl;
				hFilter->SetWindowSize(1); // fine Radius
				windowSize = "3x3x3";
				break;
			}
			case 2: {
				cout << "Setting window size 5x5x5" << endl;
				hFilter->SetWindowSize(3); // medium radius
				windowSize = "5x5x5";
				break;
			}
			case 3: {
				cout << "Setting window size 7x7x7" << endl;
				hFilter->SetWindowSize(5); // coarse radius
				windowSize = "7x7x7";
				break;
			}
		}

		// Set parameters
		hFilter->SetApplyMedianFiltering(false);
		hFilter->SetNormalizeTextures(false);
		hFilter->SetNumberOfBins(nBins);
		hFilter->SetStructureLabel(labels[n]);

		// Launch the texture computation
		//-------------------------------
		try
		{
			std::cout << " Executing filters.... " << std::endl;
			hFilter->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << " Error executing texture extractor " << ex.GetDescription() << std::endl;
			return -1;
		}

		std::cout << " Destination folder is " << destFolder << std::endl;
		// Store the texture images
		//--------------------------
		// Create writers
		typedef itk::ImageFileWriter< InImageType > ImageWriterType;
		ImageWriterType::Pointer writerEnergy = ImageWriterType::New();
		ImageWriterType::Pointer writerEntropy = ImageWriterType::New();
		ImageWriterType::Pointer writerCorrelation = ImageWriterType::New();
		ImageWriterType::Pointer writerDiffMoment = ImageWriterType::New();
		ImageWriterType::Pointer writerInertia = ImageWriterType::New();
		// Create files
		std::string fnameEnergy = destFolder + "//energy.mhd";
		std::string fnameEntropy = destFolder + "//entropy.mhd";
		std::string fnameCorrelation = destFolder + "//correlation.mhd";
		std::string fnameDiffMoment = destFolder + "//homogeneity.mhd";
		std::string fnameInertia = destFolder + "//contrast.mhd";
		//Write energy
		std::cout << " Writing images .... " << fnameEnergy << std::endl;
		writerEnergy->SetInput(hFilter->GetEnergyImage());
		writerEnergy->SetFileName(fnameEnergy.c_str());
		writerEnergy->Update();
		if (toDicom(fnameEnergy, destFolder + "//Energy", m_patientName, m_patientID, m_studyDate, m_imageType + " Energy", m_studyUID, "Texture analysis - " + windowSize + " - " + nBinsStr + " bins") != 0)
		{
			return EXIT_FAILURE;
		}
		//Write entropy
		std::cout << " Writing image " << fnameEntropy << std::endl;
		writerEntropy->SetInput(hFilter->GetEntropyImage());
		writerEntropy->SetFileName(fnameEntropy.c_str());
		writerEntropy->Update();
		if (toDicom(fnameEntropy, destFolder + "//Entropy", m_patientName, m_patientID, m_studyDate, m_imageType + " Entropy", m_studyUID, "Texture analysis - " + windowSize + " - " + nBinsStr + " bins") != 0)
		{
			return EXIT_FAILURE;
		}
		//Write correlation
		std::cout << " Writing image " << fnameCorrelation << std::endl;
		writerCorrelation->SetInput(hFilter->GetCorrelationImage());
		writerCorrelation->SetFileName(fnameCorrelation.c_str());
		writerCorrelation->Update();
		if (toDicom(fnameCorrelation, destFolder + "//Correlation", m_patientName, m_patientID, m_studyDate, m_imageType + " Correlation", m_studyUID, "Texture analysis - " + windowSize + " - " + nBinsStr + " bins") != 0)
		{
			return EXIT_FAILURE;
		}
		//Write difference moment
		std::cout << " Writing image " << fnameDiffMoment << std::endl;
		writerDiffMoment->SetInput(hFilter->GetDifferenceMomentImage());
		writerDiffMoment->SetFileName(fnameDiffMoment.c_str());
		writerDiffMoment->Update();
		if (toDicom(fnameDiffMoment, destFolder + "//Homogeneity", m_patientName, m_patientID, m_studyDate, m_imageType + " Homogeneity", m_studyUID, "Texture analysis - " + windowSize + " - " + nBinsStr + " bins") != 0)
		{
			return EXIT_FAILURE;
		}
		//Write inertia
		std::cout << " Writing image " << fnameInertia << std::endl;
		writerInertia->SetInput(hFilter->GetInertiaImage());
		writerInertia->SetFileName(fnameInertia.c_str());
		writerInertia->Update();
		if (toDicom(fnameInertia, destFolder + "//Contrast", m_patientName, m_patientID, m_studyDate, m_imageType + " Contrast", m_studyUID, "Texture analysis - " + windowSize + " - " + nBinsStr + " bins") != 0)
		{
			return EXIT_FAILURE;
		}
	}

	std::cout << "The end !! " << std::endl;

	return 0;
}
