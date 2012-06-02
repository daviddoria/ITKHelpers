/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "ITKHelpers.h" // make syntax parser happy
#include "ITKStatistics.h"

// STL
#include <iomanip> // for setfill()
#include <stdexcept>

// ITK
#include "itkAddImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkJoinImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

// Custom
#include "Helpers/Helpers.h"
#include "Helpers/Statistics.h"
#include "ITKTypeTraits.h"

namespace ITKHelpers
{

template<typename TImage>
bool HasNeighborWithValue(const itk::Index<2>& pixel, const TImage* const image,
                          const typename TImage::PixelType& value)
{
  std::vector<itk::Offset<2> > offsets = Get8NeighborOffsets();

  for(unsigned int neighborId = 0; neighborId < offsets.size(); ++neighborId)
    {
    if(image->GetPixel(pixel + offsets[neighborId]) == value)
      {
      return true;
      }
    }

  return false;
}

/** Copy the input to the output*/
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::Image<TInputPixel, 2>* const input, itk::Image<TOutputPixel, 2>* const output)
{
  //std::cout << "DeepCopy()" << std::endl;
  if(output->GetLargestPossibleRegion() != input->GetLargestPossibleRegion())
    {
    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();
    }
  DeepCopyInRegion(input, input->GetLargestPossibleRegion(), output);
}

// This is a specialization that ensures that the number of pixels per component also matches.
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::VectorImage<TInputPixel, 2>* const input, itk::VectorImage<TOutputPixel, 2>* const output)
{
  //std::cout << "DeepCopy<FloatVectorImageType>()" << std::endl;
  bool changed = false;
  if(input->GetNumberOfComponentsPerPixel() != output->GetNumberOfComponentsPerPixel())
    {
    output->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
    //std::cout << "Set output NumberOfComponentsPerPixel to " << input->GetNumberOfComponentsPerPixel() << std::endl;
    changed = true;
    }

  if(input->GetLargestPossibleRegion() != output->GetLargestPossibleRegion())
    {
    output->SetRegions(input->GetLargestPossibleRegion());
    changed = true;
    }
  if(changed)
    {
    output->Allocate();
    }

  DeepCopyInRegion(input, input->GetLargestPossibleRegion(), output);
}

template<typename TInputImage, typename TOutputImage>
void DeepCopyInRegion(const TInputImage* input, const itk::ImageRegion<2>& region, TOutputImage* output)
{
  // This function assumes that the size of input and output are the same.

  itk::ImageRegionConstIterator<TInputImage> inputIterator(input, region);
  itk::ImageRegionIterator<TOutputImage> outputIterator(output, region);

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

template<typename TImage>
void ReplaceValue(TImage* image, const typename TImage::PixelType& queryValue, const typename TImage::PixelType& replacementValue)
{
  // This function replaces all pixels in 'image' equal to 'queryValue' with 'replacementValue'
  itk::ImageRegionIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == queryValue)
      {
      imageIterator.Set(replacementValue);
      }
    ++imageIterator;
    }
}

template <class T>
std::vector<T> MaxValuesVectorImage(const itk::VectorImage<T, 2>* const image)
{
  typedef itk::VectorImage<T, 2> VectorImageType;
  typedef itk::Image<T, 2> ScalarImageType;

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetInput(image);

  std::vector<T> maxValues;
  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
    {
    indexSelectionFilter->SetIndex(channel);
    indexSelectionFilter->Update();
    maxValues.push_back(MaxValue<ScalarImageType>(indexSelectionFilter->GetOutput()));
    }
  return maxValues;
}

template <class TImage>
float MinValue(const TImage* const image, const itk::ImageRegion<2>& region)
{
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(region);
  extractFilter->SetInput(image);
  extractFilter->Update();

  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(extractFilter->GetOutput());
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMinimum();
}

template <class TImage>
float MaxValue(const TImage* const image, const itk::ImageRegion<2>& region)
{
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(region);
  extractFilter->SetInput(image);
  extractFilter->Update();

  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(extractFilter->GetOutput());
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMaximum();
}

template <class TImage>
float MinValue(const TImage* const image)
{
  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMinimum();
}

template <class TImage>
float MaxValue(const TImage* const image)
{
  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMaximum();
}

template <class TImage>
itk::Index<2> MaxValueLocation(const TImage* const image)
{
  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetIndexOfMaximum();
}

template <class TImage>
itk::Index<2> MinValueLocation(const TImage* const image)
{
  typedef typename itk::MinimumMaximumImageCalculator<TImage>
          ImageCalculatorFilterType;

  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetIndexOfMinimum();
}


template <class TImage>
void CopyPatchIntoImage(const TImage* patch, TImage* const image, const itk::Index<2>& centerPixel)
{
  // This function copies 'patch' into 'image' centered at 'position'.

  // The PasteFilter expects the lower left corner of the destination position, but we have passed the center pixel.
  itk::Index<2> cornerPixel;
  cornerPixel[0] = centerPixel[0] - patch->GetLargestPossibleRegion().GetSize()[0]/2;
  cornerPixel[1] = centerPixel[1] - patch->GetLargestPossibleRegion().GetSize()[1]/2;

  typedef itk::PasteImageFilter <TImage, TImage> PasteImageFilterType;

  typename PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
  pasteFilter->SetInput(0, image);
  pasteFilter->SetInput(1, patch);
  pasteFilter->SetSourceRegion(patch->GetLargestPossibleRegion());
  pasteFilter->SetDestinationIndex(cornerPixel);
  pasteFilter->InPlaceOn();
  pasteFilter->Update();

  image->Graft(pasteFilter->GetOutput());

}

template <class TImage>
void CopySelfRegion(TImage* const image, const itk::ImageRegion<2>& sourceRegion, const itk::ImageRegion<2>& targetRegion)
{
  CopyRegion(image, image, sourceRegion, targetRegion);
}

template <class TImage>
void CopyRegion(const TImage* const sourceImage, const itk::ImageRegion<2>& sourceRegion,
                TImage* const targetImage, const itk::ImageRegion<2>& targetRegion)
{
  CopyRegion(sourceImage, targetImage, sourceRegion, targetRegion);
}

template <class TImage>
void CopyRegion(const TImage* sourceImage, TImage* targetImage,
                const itk::ImageRegion<2>& sourceRegion, const itk::ImageRegion<2>& targetRegion)
{
  if(targetRegion.GetSize() != sourceRegion.GetSize())
    {
    std::cerr << "Can't copy regions that aren't the same size!" << std::endl;
    return;
    }

  itk::ImageRegionConstIterator<TImage> sourceIterator(sourceImage, sourceRegion);
  itk::ImageRegionIterator<TImage> targetIterator(targetImage, targetRegion);

  while(!sourceIterator.IsAtEnd())
    {
    targetIterator.Set(sourceIterator.Get());

    ++sourceIterator;
    ++targetIterator;
    }
}

template <class TImage>
void CopyRegion(const TImage* sourceImage, TImage* targetImage,
               const itk::Index<2>& sourcePosition, const itk::Index<2>& targetPosition, const unsigned int radius)
{
  // Copy a patch of radius 'radius' centered at 'sourcePosition' from 'sourceImage' to 'targetImage' centered at 'targetPosition'
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(GetRegionInRadiusAroundPixel(sourcePosition, radius));
  extractFilter->SetInput(sourceImage);
  extractFilter->Update();

  CopyPatchIntoImage<TImage>(extractFilter->GetOutput(), targetImage, targetPosition);
}

template <typename TInputImage, typename TOutputImage>
void ColorToGrayscale(const TInputImage* const colorImage, TOutputImage* const grayscaleImage)
{
  grayscaleImage->SetRegions(colorImage->GetLargestPossibleRegion());
  grayscaleImage->Allocate();

  typename itk::ImageRegionConstIterator<TInputImage> colorImageIterator(colorImage, colorImage->GetLargestPossibleRegion());
  typename itk::ImageRegionIterator<TOutputImage> grayscaleImageIterator(grayscaleImage, grayscaleImage->GetLargestPossibleRegion());

  typename TInputImage::PixelType largestPixel;
  largestPixel.Fill(255);

  float largestNorm = largestPixel.GetNorm();

  while(!colorImageIterator.IsAtEnd())
    {
    grayscaleImageIterator.Set(colorImageIterator.Get().GetNorm()*(255./largestNorm));

    ++colorImageIterator;
    ++grayscaleImageIterator;
    }
}


template<typename TImage>
void SetRegionToConstant(TImage* image, const itk::ImageRegion<2>& region, const typename TImage::PixelType& value)
{
  typename itk::ImageRegionIterator<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(value);

    ++imageIterator;
    }
}

template<typename TImage>
void SetImageToConstant(TImage* image, const typename TImage::PixelType& constant)
{
  SetRegionToConstant<TImage>(image, image->GetLargestPossibleRegion(), constant);
}

template<typename TImage>
unsigned int CountNonZeroPixels(const TImage* image)
{
  typename itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  unsigned int numberOfNonZeroPixels = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      numberOfNonZeroPixels++;
      }

    ++imageIterator;
    }
  return numberOfNonZeroPixels;
}

template<typename TImage>
std::vector<itk::Index<2> > GetNonZeroPixels(const TImage* image)
{
  return GetNonZeroPixels<TImage>(image, image->GetLargestPossibleRegion());
}

template<typename TImage>
std::vector<itk::Index<2> > GetNonZeroPixels(const TImage* image, const itk::ImageRegion<2>& region)
{
  std::vector<itk::Index<2> > nonZeroPixels;

  typename itk::ImageRegionConstIterator<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      nonZeroPixels.push_back(imageIterator.GetIndex());
      }

    ++imageIterator;
    }
  return nonZeroPixels;
}


template<typename TImage>
void WriteRegionAsImage(const TImage* image, const itk::ImageRegion<2>& region, const std::string& filename)
{
  // This function varies from WriteRegion() in that the Origin of the output image is (0,0).
  // Because of this, the region cannot be overlayed on the original image, but can be easily compared to other regions.
  //std::cout << "WriteRegion() " << filename << std::endl;
  //std::cout << "region " << region << std::endl;

  typedef itk::RegionOfInterestImageFilter<TImage, TImage> RegionOfInterestImageFilterType;

  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  itk::Point<float, 2> origin;
  origin.Fill(0);
  regionOfInterestImageFilter->GetOutput()->SetOrigin(origin);

  //std::cout << "regionOfInterestImageFilter " << regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

  typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();
  writer->SetFileName(filename);
  writer->SetInput(regionOfInterestImageFilter->GetOutput());
  writer->Update();
}



template<typename TImage>
void WriteRegionUnsignedChar(const TImage* image, const itk::ImageRegion<2>& region, const std::string& filename)
{
  // The file that is output has Origin = (0,0) because of how VectorImageToRGBImage copies the image.

  //std::cout << "WriteRegionUnsignedChar() " << filename << std::endl;
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> RegionOfInterestImageFilterType;

  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  //std::cout << "regionOfInterestImageFilter " << regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  VectorImageToRGBImage(regionOfInterestImageFilter->GetOutput(), rgbImage);

  //std::cout << "rgbImage " << rgbImage->GetLargestPossibleRegion() << std::endl;

  typename itk::ImageFileWriter<RGBImageType>::Pointer writer = itk::ImageFileWriter<RGBImageType>::New();
  writer->SetFileName(filename);
  writer->SetInput(rgbImage);
  writer->Update();
}

template<typename TImage>
void BlankAndOutlineRegion(TImage* image, const itk::ImageRegion<2>& region, const typename TImage::PixelType& blankValue, const typename TImage::PixelType& outlineValue)
{
  SetRegionToConstant<TImage>(image, region, blankValue);
  OutlineRegion<TImage>(image, region, outlineValue);
}

template<typename TImage>
void OutlineRegion(TImage* image, const itk::ImageRegion<2>& region, const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > boundaryPixels = GetBoundaryPixels(region);

  for(unsigned int i = 0; i < boundaryPixels.size(); ++i)
    {
    image->SetPixel(boundaryPixels[i], value);
    }
}

template<typename TImage>
void InitializeImage(TImage* image, const itk::ImageRegion<2>& region)
{
  image->SetRegions(region);
  image->Allocate();
  //image->FillBuffer(0);
  image->FillBuffer(itk::NumericTraits<typename TImage::PixelType>::Zero);
}

template<typename TImage>
void InitializeImage(const itk::VectorImage<TImage>* const image, const itk::ImageRegion<2>& region)
{
  image->SetRegions(region);
  image->Allocate();

  itk::VariableLengthVector<typename TImage::InternalPixelType> v(image->GetNumberOfComponentsPerPixel());
  image->FillBuffer(v);
}

template<typename TVectorImage>
void AnisotropicBlurAllChannels(const TVectorImage* image, TVectorImage* output, const float sigma)
{
  typedef itk::Image<typename TVectorImage::InternalPixelType, 2> ScalarImageType;

  // Disassembler
  typedef itk::VectorIndexSelectionCastImageFilter<TVectorImage, ScalarImageType> IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetInput(image);

  // Reassembler
  typedef itk::ComposeImageFilter<ScalarImageType> ImageToVectorImageFilterType;
  typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

  std::vector< typename ScalarImageType::Pointer > filteredImages;

  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->Update();

    typename ScalarImageType::Pointer imageChannel = ScalarImageType::New();
    DeepCopy<ScalarImageType>(indexSelectionFilter->GetOutput(), imageChannel);

    typedef itk::BilateralImageFilter<ScalarImageType, ScalarImageType>  BilateralFilterType;
    typename BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
    bilateralFilter->SetInput(imageChannel);
    bilateralFilter->SetDomainSigma(sigma);
    bilateralFilter->SetRangeSigma(sigma);
    bilateralFilter->Update();

    typename ScalarImageType::Pointer blurred = ScalarImageType::New();
    DeepCopy<ScalarImageType>(bilateralFilter->GetOutput(), blurred);

    filteredImages.push_back(blurred);
    imageToVectorImageFilter->SetInput(i, filteredImages[i]);
    }

  imageToVectorImageFilter->Update();

  DeepCopy<TVectorImage>(imageToVectorImageFilter->GetOutput(), output);
}


template<typename TImage>
void DilateImage(const TImage* const image, TImage* const dilatedImage, const unsigned int radius)
{
  typedef itk::BinaryBallStructuringElement<typename TImage::PixelType, 2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <TImage, TImage, StructuringElementType> BinaryDilateImageFilterType;

  typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(image);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  DeepCopy<TImage>(dilateFilter->GetOutput(), dilatedImage);

}

template<typename TImage>
void ChangeValue(TImage* const image, const typename TImage::PixelType& oldValue, const typename TImage::PixelType& newValue)
{
  itk::ImageRegionIterator<TImage> iterator(image, image->GetLargestPossibleRegion());

  while(!iterator.IsAtEnd())
    {
    if(iterator.Get() == oldValue)
      {
      iterator.Set(newValue);
      }
    ++iterator;
    }
}

template<typename TPixel>
void ExtractChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                    itk::Image<TPixel, 2>* const output)
{
  typedef itk::VectorImage<TPixel, 2> VectorImageType;
  typedef itk::Image<TPixel, 2> ScalarImageType;

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(channel);
  indexSelectionFilter->SetInput(image);
  indexSelectionFilter->Update();

  DeepCopy(indexSelectionFilter->GetOutput(), output);
}

template<typename TImage>
void ExtractChannels(const TImage* const image, const std::vector<unsigned int> channels,
                    TImage* const output)
{
  if(output->GetNumberOfComponentsPerPixel() != channels.size())
  {
    std::stringstream ss;
    ss << "ExtractChannels: 'output' must be presized. output has " << output->GetNumberOfComponentsPerPixel()
       << " components and there are " << channels.size() << " channels";
    throw std::runtime_error(ss.str());
  }

  if(channels.size() == 0)
  {
    throw std::runtime_error("Cannot extract 0 channels!");
  }

  for(unsigned int i = 0; i < channels.size(); ++i)
  {
    if(channels[i] > image->GetNumberOfComponentsPerPixel() - 1)
    {
      std::stringstream ss;
      ss << "Cannot extract channel " << channels[i] << " from an image with " << image->GetNumberOfComponentsPerPixel() << " channels!";
      throw std::runtime_error(ss.str());
    }
  }

  typedef itk::Image<float, 2> ScalarImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<TImage, ScalarImageType> IndexSelectionType;

  for(unsigned int channelIndex = 0; channelIndex < channels.size(); ++channelIndex)
  {
    std::cout << "Extracting channel " << channels[channelIndex] << " and setting channel " << channelIndex << std::endl;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(channels[channelIndex]);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    SetChannel(output, channelIndex, indexSelectionFilter->GetOutput());
  }

}

template<typename TPixel>
void ScaleChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                  const TPixel channelMax, itk::VectorImage<TPixel, 2>* const output)
{
  typedef itk::VectorImage<TPixel, 2> VectorImageType;
  typedef itk::Image<TPixel, 2> ScalarImageType;

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetIndex(channel);
  indexSelectionFilter->SetInput(image);
  indexSelectionFilter->Update();

  typedef itk::RescaleIntensityImageFilter< ScalarImageType, ScalarImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(indexSelectionFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(channelMax);
  rescaleFilter->Update();

  DeepCopy(image, output);

  ReplaceChannel(output, channel, rescaleFilter->GetOutput(), output);
}

template<typename TPixel>
void ReplaceChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                    const itk::Image<TPixel, 2>* const replacement, itk::VectorImage<TPixel, 2>* const output)
{
  if(image->GetLargestPossibleRegion() != replacement->GetLargestPossibleRegion())
    {
    throw std::runtime_error("Image and replacement channel are not the same size!");
    }

  DeepCopy(image, output);

  itk::ImageRegionConstIterator<itk::VectorImage<TPixel, 2> > iterator(image, image->GetLargestPossibleRegion());

  while(!iterator.IsAtEnd())
    {
    typename itk::VectorImage<TPixel, 2>::PixelType pixel = iterator.Get();
    pixel[channel] = replacement->GetPixel(iterator.GetIndex());
    output->SetPixel(iterator.GetIndex(), pixel);
    ++iterator;
    }
}

template<typename TImage>
void ReadImage(const std::string& fileName, TImage* const image)
{
  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();

  ITKHelpers::DeepCopy(reader->GetOutput(), image);
}


template<typename TInputImage, typename TOutputImage, typename TFilter>
void FilterImage(const TInputImage* input, TOutputImage* output)
{
  typename TFilter::Pointer filter = TFilter::New();
  filter->SetInput(input);
  filter->Update();
  DeepCopy(filter->GetOutput(), output);
}

template<typename TImage>
void NormalizeVectorImage(TImage* const image)
{
  // TImage::PixelType must have a .Normalize() function

  itk::ImageRegionIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    typename TImage::PixelType pixel = imageIterator.Get();
    pixel.Normalize();
    imageIterator.Set(pixel);
    ++imageIterator;
    }
}


template<typename T>
itk::Index<2> CreateIndex(const T& v)
{
  itk::Index<2> index = {{v[0], v[1]}};
  return index;
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerComponentType
AverageNeighborValue(const TImage* const image, const itk::Index<2>& pixel)
{
  itk::ImageRegion<2> neighborhoodRegion = GetRegionInRadiusAroundPixel(pixel, 1);
  neighborhoodRegion.Crop(image->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > neighbors = Get8NeighborsInRegion(neighborhoodRegion, pixel);

  typename TImage::PixelType pixelSum = 0.0f;

  for(unsigned int i = 0; i < neighbors.size(); ++i)
    {
    pixelSum += image->GetPixel(neighbors[i]);
    }
  return pixelSum / static_cast<float>(neighbors.size());
}

/** Get a list of the valid neighbors of a pixel.*/
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsWithValue(const itk::Index<2>& pixel, const TImage* const image,
                                                  const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > neighbors = Get8NeighborsInRegion(GetRegionInRadiusAroundPixel(pixel, 1), pixel);
  std::vector<itk::Index<2> > neighborsWithValue;
  for(unsigned int i = 0; i < neighbors.size(); ++i)
    {
    if(image->GetPixel(neighbors[i]) == value)
      {
      neighborsWithValue.push_back(neighbors[i]);
      }
    }
  return neighborsWithValue;
}

template<typename TImage>
std::vector<itk::Index<2> > GetPixelsWithValue(const TImage* const image,
                                               const typename TImage::PixelType& value)
{
  return GetPixelsWithValue(image, image->GetLargestPossibleRegion(), value);
}

template<typename TImage>
std::vector<itk::Index<2> > GetPixelsWithValue(const TImage* const image, const itk::ImageRegion<2>& region,
                                               const typename TImage::PixelType& value)
{
  std::vector<itk::Index<2> > pixelsWithValue;

  itk::ImageRegionConstIterator<TImage> regionIterator(image, region);
  while(!regionIterator.IsAtEnd())
    {
    if(regionIterator.Get() == value)
      {
      pixelsWithValue.push_back(regionIterator.GetIndex());
      }
    ++regionIterator;
    }

  return pixelsWithValue;
}

template<typename TImage>
unsigned int CountPixelsWithValue(const TImage* const image, const typename TImage::PixelType& value)
{
  return GetPixelsWithValue(image, image->GetLargestPossibleRegion(), value).size();
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageOfPixelsAtIndices(const TImage* const image, const std::vector<itk::Index<2> >& indices)
{
  std::vector<typename TImage::PixelType> pixels;
  for(unsigned int i = 0; i < indices.size(); ++i)
  {
    pixels.push_back(image->GetPixel(indices[i]));
  }
  using Statistics::Average;
  //using ITKHelpers::Average;
  return Average(pixels);
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType VarianceOfPixelsAtIndices(const TImage* const image, const std::vector<itk::Index<2> >& indices)
{
  std::vector<typename TImage::PixelType> pixels;
  for(unsigned int i = 0; i < indices.size(); ++i)
    {
    pixels.push_back(image->GetPixel(indices[i]));
    }

  return Statistics::Variance(pixels);
}

// This lets you sum either a scalar or a VariableLengthVector
template <typename T>
typename TypeTraits<T>::LargerComponentType SumOfComponents(const T& v)
{
  typename TypeTraits<T>::LargerComponentType sumOfComponents;
  SetObjectToZero(sumOfComponents);
  using Helpers::length;
  using ITKHelpers::length;
  using Helpers::index;
  using ITKHelpers::index;
  for(unsigned int i = 0; i < length(v); ++i)
    {
    sumOfComponents += index(v, i);
    }

  return sumOfComponents;
}

template <typename T>
typename TypeTraits<T>::LargerComponentType SumOfComponentMagnitudes(const T& v)
{
  typename TypeTraits<T>::LargerComponentType sumOfComponents;
  SetObjectToZero(sumOfComponents);
  using Helpers::length;
  using ITKHelpers::length;
  using Helpers::index;
  using ITKHelpers::index;
  for(unsigned int i = 0; i < length(v); ++i)
    {
    sumOfComponents += fabs(index(v, i));
    }

  return sumOfComponents;
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageOfImage(const TImage* const image)
{
  return AverageInRegion(image, image->GetLargestPossibleRegion());
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageInRegion(const TImage* const image,
                                                                            const itk::ImageRegion<2>& region)
{
  typename itk::ImageRegionConstIterator<TImage> imageIterator(image, region);
  std::vector<typename TImage::PixelType> pixels;
  while(!imageIterator.IsAtEnd())
    {
    pixels.push_back(imageIterator.Get());
    ++imageIterator;
    }

  using Statistics::Average;
  //using ITKHelpers::Average;
  return Average(pixels);
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType VarianceOfImage(const TImage* const image)
{
  return VarianceInRegion(image, image->GetLargestPossibleRegion());
}

template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType VarianceInRegion(const TImage* const image,
                                                                             const itk::ImageRegion<2>& region)
{
  typename itk::ImageRegionConstIterator<TImage> imageIterator(image, region);
  std::vector<typename TImage::PixelType> pixels;
  while(!imageIterator.IsAtEnd())
    {
    pixels.push_back(imageIterator.Get());
    ++imageIterator;
    }

  return ITKStatistics::Variance(pixels);
}


template<typename TImage, typename TDifferenceFunctor>
float AverageDifferenceInRegion(const TImage* const image, const itk::ImageRegion<2>& region1,
                                const itk::ImageRegion<2>& region2, TDifferenceFunctor differenceFunctor)
{
  return AverageDifferenceInRegion(image, region1, image, region2, differenceFunctor);
}

template<typename TImage, typename TDifferenceFunctor>
float AverageDifferenceInRegion(const TImage* const image1, const itk::ImageRegion<2>& region1,
                                const TImage* const image2, const itk::ImageRegion<2>& region2,
                                TDifferenceFunctor differenceFunctor)
{
  if(region1.GetSize() != region2.GetSize())
    {
    throw std::runtime_error("Regions must be the same size to compare!");
    }

  itk::ImageRegionConstIterator<TImage> region1Iterator(image1, region1);
  itk::ImageRegionConstIterator<TImage> region2Iterator(image2, region2);

  float totalDifference = 0.0f;
  while(!region1Iterator.IsAtEnd())
    {
    totalDifference += differenceFunctor(region1Iterator.Get(), region2Iterator.Get());

    ++region1Iterator;
    ++region2Iterator;
    }
  float averageDifference = totalDifference / static_cast<float>(region1.GetNumberOfPixels());
  return averageDifference;
}

template<typename T>
unsigned int length(const itk::RGBPixel<T>& v)
{
  return 3;
}

template<typename T, unsigned int N>
unsigned int length(const itk::CovariantVector<T, N>& v)
{
  return v.Dimension;
}

template<typename T>
unsigned int length(const itk::VariableLengthVector<T>& v)
{
  return v.GetSize();
}

template<typename T>
T index(const itk::RGBPixel<T>& v, size_t i)
{
  return v[i];
}

template<typename T>
T& index(itk::RGBPixel<T>& v, size_t i)
{
  return v[i];
}

template<typename T, unsigned int N>
T& index(itk::CovariantVector<T, N>& v, size_t i)
{
  return v[i];
}

template<typename T, unsigned int N>
T index(const itk::CovariantVector<T, N>& v, size_t i)
{
  return v[i];
}

template<typename T>
T& index(itk::VariableLengthVector<T>& v, size_t i)
{
  return v[i];
}

template<typename T>
T index(const itk::VariableLengthVector<T>& v, size_t i)
{
  return v[i];
}

template<typename T>
void SetObjectToZero(T& object)
{
  using Helpers::index;
  using ITKHelpers::index;
  using Helpers::length;
  using ITKHelpers::length;
  for(unsigned int i = 0; i < length(object); ++i)
    {
    index(object, i) = 0;
    }
}

template<typename TImage>
void SubtractRegions(const TImage* const image1, const itk::ImageRegion<2>& region1, const TImage* const image2, const itk::ImageRegion<2>& region2, TImage* const output)
{
  assert(region1.GetSize() == region2.GetSize());

  itk::ImageRegion<2> outputRegion(ZeroIndex(), region1.GetSize());
  output->SetRegions(outputRegion);
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> image1Iterator(image1, region1);
  itk::ImageRegionConstIterator<TImage> image2Iterator(image2, region2);

  while(!image1Iterator.IsAtEnd())
    {
    float difference = image1Iterator.Get() - image2Iterator.Get();
    itk::Index<2> index = Helpers::ConvertFrom<itk::Index<2>, itk::Offset<2> >(image1Iterator.GetIndex() -
                                                                               image1->GetLargestPossibleRegion().GetIndex());
    output.SetPixel(index, difference);
    ++image1Iterator;
    ++image2Iterator;
    }
}

template<typename TImage>
void ANDRegions(const TImage* const image1, const itk::ImageRegion<2>& region1, const TImage* const image2,
                const itk::ImageRegion<2>& region2, itk::Image<bool, 2>* const output)
{
  assert(region1.GetSize() == region2.GetSize());

  itk::ImageRegion<2> outputRegion(ZeroIndex(), region1.GetSize());
  output->SetRegions(outputRegion);
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> image1Iterator(image1, region1);
  itk::ImageRegionConstIterator<TImage> image2Iterator(image2, region2);

  while(!image1Iterator.IsAtEnd())
    {
    bool result = image1Iterator.Get() && image2Iterator.Get();
    itk::Index<2> index = Helpers::ConvertFrom<itk::Index<2>, itk::Offset<2> >(image1Iterator.GetIndex() -
                                                                               image1->GetLargestPossibleRegion().GetIndex());
    output->SetPixel(index, result);
    ++image1Iterator;
    ++image2Iterator;
    }
}

template<typename TImage>
void XORRegions(const TImage* const image1, const itk::ImageRegion<2>& region1, const TImage* const image2,
                const itk::ImageRegion<2>& region2, itk::Image<bool, 2>* const output)
{
  assert(region1.GetSize() == region2.GetSize());

  itk::ImageRegion<2> outputRegion(ZeroIndex(), region1.GetSize());
  output->SetRegions(outputRegion);
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> image1Iterator(image1, region1);
  itk::ImageRegionConstIterator<TImage> image2Iterator(image2, region2);

  while(!image1Iterator.IsAtEnd())
    {
    bool result = image1Iterator.Get() ^ image2Iterator.Get(); // ^ is XOR (exclusive OR)
    itk::Index<2> index = Helpers::ConvertFrom<itk::Index<2>, itk::Offset<2> >(image1Iterator.GetIndex() - region1.GetIndex());
    output->SetPixel(index, result);
    ++image1Iterator;
    ++image2Iterator;
    }
}

template<typename TImage>
void ExtractRegion(const TImage* const image, const itk::ImageRegion<2>& region,
                   TImage* const output)
{
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> ExtractFilterType;

  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetRegionOfInterest(region);
  extractFilter->SetInput(image);
  extractFilter->Update();

  DeepCopy(extractFilter->GetOutput(), output);
}

template<typename TImage>
void PrintImage(const TImage* const image)
{
  PrintRegion(image, image->GetLargestPossibleRegion());
}

template<typename TImage>
void PrintRegion(const TImage* const image, const itk::ImageRegion<2>& region)
{
  itk::ImageRegionConstIterator<TImage> imageIterator(image, region);

  unsigned int counter = 0;
  while(!imageIterator.IsAtEnd())
    {
    std::cout << static_cast<typename TypeTraits<typename TImage::PixelType>::LargerType >(imageIterator.Get()) << " ";
    counter++;
    if(counter == region.GetSize()[0])
      {
      counter = 0;
      std::cout << std::endl;
      }
    ++imageIterator;
    }
}

template<typename TImage>
std::vector<typename TImage::PixelType> GetPixelValues(const TImage* const image, const std::vector<itk::Index<2> >& indices)
{
  std::vector<typename TImage::PixelType> values;
  for(std::vector<itk::Index<2> >::const_iterator iter = indices.begin(); iter != indices.end(); ++iter)
  {
    values.push_back(image->GetPixel(*iter));
  }

  return values;
}

// This does not work for itk::VectorImage
// template<typename TImage>
// void StackImages(const TImage* const image1, const TImage* const image2, const TImage* const output)
// {
//
//   typedef itk::JoinImageFilter<TImage, TImage> JoinImageFilterType;
//   typename JoinImageFilterType::Pointer joinFilter = JoinImageFilterType::New();
//   joinFilter->SetInput1(image1);
//   joinFilter->SetInput2(image2);
//   joinFilter->Update();
//
//   DeepCopy(joinFilter->GetOutput(), output);
//}

template<typename TVectorImage, typename TScalarImage>
void SetChannel(TVectorImage* const vectorImage, const unsigned int channel, const TScalarImage* const image)
{
  assert(vectorImage->GetLargestPossibleRegion() == image->GetLargestPossibleRegion());

  if(channel > vectorImage->GetNumberOfComponentsPerPixel() - 1)
  {
    std::stringstream ss;
    ss << "SetChannel: Trying to set channel " << channel << " but vectorImage only has "
       << vectorImage->GetNumberOfComponentsPerPixel() << " channels.";
    throw std::runtime_error(ss.str());
  }
  itk::ImageRegionConstIterator<TScalarImage> inputIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TVectorImage> outputIterator(vectorImage, vectorImage->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    typename TVectorImage::PixelType pixel = outputIterator.Get();
    pixel[channel] = inputIterator.Get();
    outputIterator.Set(pixel);
    ++inputIterator;
    ++outputIterator;
    }
}

template<typename TPixel>
void ConvertTo3Channel(const itk::VectorImage<TPixel, 2>* const image,
                      typename itk::VectorImage<TPixel, 2>* const output)
{
  typedef itk::VectorImage<TPixel, 2> ImageType;

  if(image->GetNumberOfComponentsPerPixel() == 3)
  {
    DeepCopy(image,output);
    return;
  }

  if(image->GetNumberOfComponentsPerPixel() > 3)
  {
    output->SetRegions(image->GetLargestPossibleRegion());
    output->SetNumberOfComponentsPerPixel(3);
    output->Allocate();
    for(unsigned int channel = 0; channel < 3; ++channel)
      {
      typename itk::Image<TPixel, 2>::Pointer outputComponent = itk::Image<TPixel, 2>::New();
      outputComponent->SetRegions(image->GetLargestPossibleRegion());
      outputComponent->Allocate();
      ExtractChannel(image, channel, outputComponent.GetPointer());

      SetChannel(output, channel, outputComponent.GetPointer());
      }
  }

  if(image->GetNumberOfComponentsPerPixel() == 1)
  {
    output->SetRegions(image->GetLargestPossibleRegion());
    output->SetNumberOfComponentsPerPixel(3);
    output->Allocate();

    typedef itk::Image<TPixel, 2> ScalarImageType;
    typename ScalarImageType::Pointer outputComponent = ScalarImageType::New();
    outputComponent->SetRegions(image->GetLargestPossibleRegion());
    outputComponent->Allocate();
    ExtractChannel(image, 0, outputComponent.GetPointer());

    itk::ImageRegionConstIterator<ScalarImageType> inputIterator(outputComponent, outputComponent->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> outputIterator(output, output->GetLargestPossibleRegion());

    while(!inputIterator.IsAtEnd())
      {
      // Copy the first (only) channel of the input to the output
      typename ImageType::PixelType pixel = outputIterator.Get();
      pixel[0] = inputIterator.Get();
      pixel[1] = inputIterator.Get();
      pixel[2] = inputIterator.Get();
      outputIterator.Set(pixel);
      ++inputIterator;
      ++outputIterator;
      }
  }
}

template<typename TImage>
void ScaleTo255(TImage* const image)
{
  typedef itk::RescaleIntensityImageFilter<TImage, TImage> RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  DeepCopy(rescaleFilter->GetOutput(), image);
}

template<typename TImage>
void ScaleAllChannelsTo255(TImage* const image)
{
  typedef itk::Image<typename TImage::InternalPixelType, 2> ScalarImageType;
  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    typename ScalarImageType::Pointer outputComponent = ScalarImageType::New();
    outputComponent->SetRegions(image->GetLargestPossibleRegion());
    outputComponent->Allocate();
    ExtractChannel(image, channel, outputComponent.GetPointer());
    ScaleTo255(outputComponent.GetPointer());
    SetChannel(image, channel, outputComponent.GetPointer());
  }
}


template <typename TImage>
void WriteSequentialImage(const TImage* const image, const std::string& filePrefix, const unsigned int iteration)
{
  std::string fileName = Helpers::GetSequentialFileName(filePrefix, iteration, "mha");

  WriteImage(image, fileName);
}

template <typename TImage>
void WriteImageConditional(const TImage* const image, const std::string& fileName, const bool condition)
{
  if(condition)
    {
    WriteImage(image, fileName);
    }
}


template <class TImage>
void WriteScaledScalarImage(const TImage* const image, const std::string& filename)
{
//   if(T::PixelType::Dimension > 1)
//     {
//     std::cerr << "Cannot write scaled scalar image with vector image input!" << std::endl;
//     return;
//     }
  typedef itk::RescaleIntensityImageFilter<TImage, UnsignedCharScalarImageType> RescaleFilterType; // expected ';' before rescaleFilter

  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  typedef itk::ImageFileWriter<UnsignedCharScalarImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();
}


template<typename TImage>
void WriteImage(const TImage* const image, const std::string& filename)
{
  // This is a convenience function so that images can be written in 1 line instead of 4.
  typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();
}


template<typename TImage>
void WriteRGBImage(const TImage* const input, const std::string& filename)
{
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBImageType;

  RGBImageType::Pointer output = RGBImageType::New();
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<RGBImageType> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    itk::CovariantVector<unsigned char, 3> pixel;
    for(unsigned int i = 0; i < 3; ++i)
      {
      pixel[i] = inputIterator.Get()[i];
      }
    outputIterator.Set(pixel);
    ++inputIterator;
    ++outputIterator;
    }

  typename itk::ImageFileWriter<RGBImageType>::Pointer writer = itk::ImageFileWriter<RGBImageType>::New();
  writer->SetFileName(filename);
  writer->SetInput(output);
  writer->Update();
}

// template<typename TImage>
// void WritePatch(const ImagePatch<TImage>& patch, const std::string& filename)
// {
//   WriteRegion<TImage>(patch.Image, patch.Region, filename);
// }
//
// template<typename TImage>
// void WriteMaskedPatch(const Mask* mask, const ImagePatch<TImage>& patch, const std::string& filename, const typename TImage::PixelType& holeColor)
// {
//   WriteMaskedRegion<TImage>(patch.Image, mask, patch.Region, filename);
// }


template<typename TImage>
void WriteRegion(const TImage* const image, const itk::ImageRegion<2>& region, const std::string& filename)
{
  //std::cout << "WriteRegion() " << filename << std::endl;
  //std::cout << "region " << region << std::endl;
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> RegionOfInterestImageFilterType;

  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  //std::cout << "regionOfInterestImageFilter " << regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

  typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();
  writer->SetFileName(filename);
  writer->SetInput(regionOfInterestImageFilter->GetOutput());
  writer->Update();
}

// Why is RegionOfInterestImageFilter used (above) instead of ExtractFilterType ?
// template <typename TImage>
// void PatchBasedInpainting<TImage>::DebugWritePatch(const itk::ImageRegion<2>& inputRegion, const std::string& filename)
// {
//   if(!this->DebugImages)
//     {
//     return;
//     }
//
//   typedef itk::RegionOfInterestImageFilter< FloatVectorImageType,
//                                             FloatVectorImageType> ExtractFilterType;
//   itk::ImageRegion<2> region = inputRegion;
//   region.Crop(this->CurrentInpaintedImage->GetLargestPossibleRegion());
//
//   ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
//   extractFilter->SetRegionOfInterest(region);
//   extractFilter->SetInput(this->CurrentInpaintedImage);
//   extractFilter->Update();
//   /*
//   typedef itk::Image<itk::CovariantVector<unsigned char, TImage::PixelType::Dimension>, 2> OutputImageType;
//
//   typedef itk::CastImageFilter< TImage, OutputImageType > CastFilterType;
//   typename CastFilterType::Pointer castFilter = CastFilterType::New();
//   castFilter->SetInput(extractFilter->GetOutput());
//   castFilter->Update();
//
//   Helpers::WriteImage<OutputImageType>(castFilter->GetOutput(), filename);
//   */
//
// }

template <typename T>
void OutputVector(const std::vector<T>& v)
{
  for(unsigned int i = 0; i < v.size(); i++)
    {
    std::cout << v[i] << " ";
    }
  std::cout << std::endl;
}

template <typename TImage, typename TPixel>
void SetPixels(TImage* const image, const std::vector<itk::Index<2> >& pixels, const TPixel& value)
{
  for(unsigned int i = 0; i < pixels.size(); ++i)
    {
    image->SetPixel(pixels[i], value);
    }
}


template <typename TImage>
void SetPixelsInRegionToValue(TImage* const image, const itk::ImageRegion<2>& region,
                              const typename TImage::PixelType& value)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(value);
    ++imageIterator;
    }
}

template<typename TImage>
void ExtractAndNormalizeRegion(const TImage* const image, const itk::ImageRegion<2>& region, TImage* const outputImage)
{
  ITKHelpers::ExtractRegion(image, region, outputImage);

  typename TypeTraits<typename TImage::PixelType>::LargerType mean = ITKHelpers::AverageOfImage(outputImage);
  typename TypeTraits<typename TImage::PixelType>::LargerType standardDeviation =
             sqrt(ITKHelpers::VarianceOfImage(outputImage));

  typedef itk::AddImageFilter<TImage, TImage, TImage> AddImageFilterType;
  typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput(outputImage);
  addImageFilter->SetConstant2(-1.0f * mean);
  addImageFilter->Update();

  DeepCopy(addImageFilter->GetOutput(), outputImage);

  typedef itk::MultiplyImageFilter<TImage, TImage, TImage> MultiplyImageFilterType;
  typename MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(outputImage);
  multiplyImageFilter->SetConstant(1.0f/standardDeviation);
  multiplyImageFilter->Update();

  DeepCopy(multiplyImageFilter->GetOutput(), outputImage);
}

template<typename TImage>
void NormalizeImageChannels(const TImage* const image, TImage* const outputImage)
{
  DeepCopy(image, outputImage);

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    std::cout << "Normalizing channel " << channel << std::endl;
    typedef itk::Image<float, 2> ScalarImageType;
    ScalarImageType::Pointer scalarImage = ScalarImageType::New();
    ExtractChannel(image, channel, scalarImage.GetPointer());

    float mean = MeanValue(scalarImage.GetPointer());
    std::cout << "Channel " << channel << " mean is " << mean << std::endl;

//     float variance = Variance(scalarImage.GetPointer());
//     std::cout << "Channel " << channel << " variance is " << variance << std::endl;
    float standardDeviation = StandardDeviation(scalarImage.GetPointer());
    std::cout << "Channel " << channel << " standardDeviation is " << standardDeviation << std::endl;

    itk::ImageRegionIterator<ScalarImageType> imageIterator(scalarImage, scalarImage->GetLargestPossibleRegion());
    while(!imageIterator.IsAtEnd())
      {
      float newValue = (imageIterator.Get() - mean) / standardDeviation;
      imageIterator.Set(newValue);
      ++imageIterator;
      }
    ReplaceChannel(outputImage, channel, scalarImage.GetPointer(), outputImage);
  }
}

template<typename TImage>
float MeanValue(const TImage* const image)
{
  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  float sum = 0.0f;
  while(!imageIterator.IsAtEnd())
    {
    sum += imageIterator.Get();

    ++imageIterator;
    }
  return sum / static_cast<float>(image->GetLargestPossibleRegion().GetNumberOfPixels());
}


template<typename TImage>
float Variance(const TImage* const image)
{
  float average = MeanValue(image);

  float channelVarianceSummation = 0.0f;

  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    channelVarianceSummation += pow(imageIterator.Get() - average, 2);
    ++imageIterator;
    }
  // This (N-1) term in the denominator is for the "unbiased" sample variance.
  // This is what is used by Matlab, Wolfram alpha, etc.
  float variance = channelVarianceSummation /
         static_cast<float>(image->GetLargestPossibleRegion().GetNumberOfPixels() - 1);

  return variance;
}

template<typename TImage>
float StandardDeviation(const TImage* const image)
{
  return sqrt(Variance(image));
}

template<typename TImage>
std::vector<typename TImage::InternalPixelType> ComputeMinOfAllChannels(const TImage* const image)
{
  std::vector<typename TImage::InternalPixelType> mins(image->GetNumberOfComponentsPerPixel(), -1);

  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<TImage, FloatScalarImageType> IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    typedef itk::MinimumMaximumImageCalculator <FloatScalarImageType> ImageCalculatorFilterType;
    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
    imageCalculatorFilter->SetImage(indexSelectionFilter->GetOutput());
    imageCalculatorFilter->Compute();

    mins[i] = imageCalculatorFilter->GetMinimum();
    }

  return mins;
}

template<typename TImage>
std::vector<typename TImage::InternalPixelType> ComputeMaxOfAllChannels(const TImage* const image)
{
  std::vector<typename TImage::InternalPixelType> maxs(image->GetNumberOfComponentsPerPixel(), -1);

  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<TImage, FloatScalarImageType> IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    typedef itk::MinimumMaximumImageCalculator <FloatScalarImageType> ImageCalculatorFilterType;
    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
    imageCalculatorFilter->SetImage(indexSelectionFilter->GetOutput());
    imageCalculatorFilter->Compute();

    maxs[i] = imageCalculatorFilter->GetMaximum();
    }

  return maxs;
}


template <class TImage>
void CentralDifferenceDerivative(const TImage* const image, const unsigned int direction, TImage* const output)
{
  typedef itk::DerivativeImageFilter<TImage, TImage> DerivativeImageFilterType;
  typename DerivativeImageFilterType::Pointer derivativeFilter = DerivativeImageFilterType::New();
  derivativeFilter->SetDirection(direction);
  derivativeFilter->SetOrder(1);
  derivativeFilter->SetInput(image);
  derivativeFilter->Update();

  DeepCopy(derivativeFilter->GetOutput(), output);
}

template <typename TPixel>
void WriteVectorImageAsRGB(const itk::VectorImage<TPixel,2>* const image, const std::string& fileName)
{
  RGBImageType::Pointer rgbImage = RGBImageType::New();
  ITKHelpers::VectorImageToRGBImage(image, rgbImage);
  WriteImage(rgbImage.GetPointer(), fileName);
}

template <typename TPixel>
void WriteVectorImageRegionAsRGB(const itk::VectorImage<TPixel,2>* const image, const itk::ImageRegion<2>& region, const std::string& filename)
{
  typedef itk::VectorImage<TPixel,2> VectorImageType;

  //std::cout << "WriteRegion() " << filename << std::endl;
  //std::cout << "region " << region << std::endl;
  typedef itk::RegionOfInterestImageFilter<VectorImageType, VectorImageType> RegionOfInterestImageFilterType;

  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  ITKHelpers::VectorImageToRGBImage(regionOfInterestImageFilter->GetOutput(), rgbImage);
  WriteImage(rgbImage.GetPointer(), filename);
}

template <typename TPixel>
void VectorImageToRGBImage(const itk::VectorImage<TPixel,2>* const image, RGBImageType* const rgbImage)
{
  typedef itk::VectorImage<TPixel,2> VectorImageType;

  // Only the first 3 components are used (assumed to be RGB)
  rgbImage->SetRegions(image->GetLargestPossibleRegion());
  rgbImage->Allocate();

  itk::ImageRegionConstIteratorWithIndex<VectorImageType> inputIterator(image, image->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    typename VectorImageType::PixelType inputPixel = inputIterator.Get();
    RGBImageType::PixelType outputPixel;
    outputPixel.SetRed(inputPixel[0]);
    outputPixel.SetGreen(inputPixel[1]);
    outputPixel.SetBlue(inputPixel[2]);

    rgbImage->SetPixel(inputIterator.GetIndex(), outputPixel);
    ++inputIterator;
    }
}

template <typename TVector>
std::string VectorToString(const TVector& vec)
{
  std::stringstream ss;
  ss << "(";
  for(unsigned int i = 0; i < vec.GetSize(); ++i)
  {
    ss << static_cast<float>(vec[i]);
    if(i == vec.GetSize() - 1)
    {
      ss << ")";
    }
    else
    {
      ss << ", ";
    }
  }
  return ss.str();
}

template <typename TImage>
void Downsample(const TImage* const image, const float factor, TImage* const output)
{
  typename TImage::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();

  typename TImage::SizeType outputSize = image->GetLargestPossibleRegion().GetSize();
  outputSize[0] /= factor;
  outputSize[1] /= factor;

  typename TImage::SpacingType outputSpacing;
  outputSpacing[0] = image->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] = image->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));

  typedef itk::IdentityTransform<double, 2> TransformType;
  typedef itk::ResampleImageFilter<TImage, TImage> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
  resampleFilter->SetInput(image);
  resampleFilter->SetSize(outputSize);
  resampleFilter->SetOutputSpacing(outputSpacing);
  resampleFilter->SetTransform(TransformType::New());
  resampleFilter->UpdateLargestPossibleRegion();

  ITKHelpers::DeepCopy(resampleFilter->GetOutput(), output);
}

template <typename TImage>
void Upsample(const TImage* const image, const float factor, TImage* const output)
{
  Downsample(image, 1.0f/factor, output);
}

template <typename TImage>
void FillDifference(const TImage* const image, const itk::ImageRegion<2>& region,
                    TImage* const output, const typename TImage::PixelType& value)
{
  // If the image is not inside the region (i.e. not smaller than the region), there is nothing to do
  if(!region.IsInside(image->GetLargestPossibleRegion()))
  {
    return;
  }

  itk::ImageRegion<2> differenceRegion = region;
  differenceRegion.Crop(image->GetLargestPossibleRegion());

  itk::ImageRegionIterator<TImage> imageIterator(output, differenceRegion);

  while(!imageIterator.IsAtEnd())
    {
    imageIterator.Set(value);
    ++imageIterator;
    }
}

template <typename TImage>
void ScaleImage(const TImage* const image, const itk::Size<2>& destinationSize, TImage* const output)
{
  typename TImage::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();

  typename TImage::SpacingType outputSpacing;
  outputSpacing[0] = image->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(destinationSize[0]));
  outputSpacing[1] = image->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(destinationSize[1]));

  typedef itk::IdentityTransform<double, 2> TransformType;
  typedef itk::ResampleImageFilter<TImage, TImage> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
  resampleFilter->SetInput(image);
  resampleFilter->SetSize(destinationSize);
  resampleFilter->SetOutputSpacing(outputSpacing);
  resampleFilter->SetTransform(TransformType::New());
  resampleFilter->UpdateLargestPossibleRegion();

  ITKHelpers::DeepCopy(resampleFilter->GetOutput(), output);
}

template <typename TImage>
void CreateEvenSizeImage(const TImage* const image, TImage* const output)
{
  itk::Size<2> evenSize = ITKHelpers::MakeSizeEven(image->GetLargestPossibleRegion().GetSize());
  itk::Index<2> corner = {{0,0}};
  itk::ImageRegion<2> evenRegion(corner, evenSize);
  output->SetRegions(evenRegion);
  output->Allocate();

  CopyRegion(image, output, evenRegion, evenRegion);
}

}// end namespace ITKHelpers
