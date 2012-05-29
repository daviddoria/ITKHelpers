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

#include "ITKHelpers.h"

// ITK
#include "itkComposeImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkImageToVectorImageFilter.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

// Helpers submodule
#include "Helpers/Helpers.h"

// Custom
#include "itkRGBToLabColorSpacePixelAccessor.h"

namespace ITKHelpers
{


std::vector<itk::Index<2> > GetIndicesInRegion(const itk::ImageRegion<2>& region)
{
  std::vector<itk::Index<2> > indices;

  typedef itk::Image<unsigned char, 2> DummyImageType;
  DummyImageType::Pointer image = DummyImageType::New();
  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionConstIteratorWithIndex<DummyImageType> regionIterator(image, region);
  while(!regionIterator.IsAtEnd())
    {
    indices.push_back(regionIterator.GetIndex());
    ++regionIterator;
    }
  return indices;
}

itk::ImageRegion<2> GetQuadrant(const itk::ImageRegion<2>& region, const unsigned int requestedQuadrant)
{
  // Note: the four quadrants might not cover the entire 'region'.

  unsigned int quadrantSideLength = region.GetSize()[0]/2;
  itk::Size<2> size = {{quadrantSideLength, quadrantSideLength}};
  itk::Index<2> corner;
  if(requestedQuadrant == 0)
  {
    corner = region.GetIndex();
  }
  else if(requestedQuadrant == 1)
  {
    itk::Offset<2> offset = {{quadrantSideLength, 0}};
    corner = region.GetIndex() + offset;
  }
  else if(requestedQuadrant == 2)
  {
    itk::Offset<2> offset = {{0, quadrantSideLength}};
    corner = region.GetIndex() + offset;
  }
  else if(requestedQuadrant == 3)
  {
    itk::Offset<2> offset = {{quadrantSideLength, quadrantSideLength}};
    corner = region.GetIndex() + offset;
  }
  else
  {
    std::stringstream ss;
    ss << "There are only 4 quadrants (0-3). Requested " << requestedQuadrant;
    throw std::runtime_error(ss.str());
  }

  itk::ImageRegion<2> quadrant(corner, size);
  return quadrant;
}

itk::Index<2> ZeroIndex()
{
  itk::Index<2> index;
  index.Fill(0);
  return index;
}

unsigned int GetNumberOfComponentsPerPixelInFile(const std::string& filename)
{
  typedef itk::VectorImage<float, 2> TestImageType;
  typedef  itk::ImageFileReader<TestImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(filename);
  imageReader->Update();

  return imageReader->GetOutput()->GetNumberOfComponentsPerPixel();
}

std::string GetIndexString(const itk::Index<2>& index)
{
  std::stringstream ss;
  ss << "(" << index[0] << ", " << index[1] << ")";
  return ss.str();
}

std::string GetSizeString(const itk::Size<2>& size)
{
  std::stringstream ss;
  ss << "(" << size[0] << ", " << size[1] << ")";
  return ss.str();
}

FloatVector2Type AverageVectors(const std::vector<FloatVector2Type>& vectors)
{
  FloatVector2Type totalVector;
  totalVector.Fill(0);

  if(vectors.size() == 0)
    {
    std::cerr << "Cannot average 0 vectors!" << std::endl;
    return totalVector;
    }

  for(unsigned int i = 0; i < vectors.size(); ++i)
    {
    totalVector[0] += vectors[i][0];
    totalVector[1] += vectors[i][1];
    }

  FloatVector2Type averageVector;
  averageVector[0] = totalVector[0] / static_cast<float>(vectors.size());
  averageVector[1] = totalVector[1] / static_cast<float>(vectors.size());

  return averageVector;
}

float AngleBetween(const FloatVector2Type& v1, const FloatVector2Type& v2)
{
  FloatVector2Type v1normalized = v1;
  v1normalized.Normalize();

  FloatVector2Type v2normalized = v2;
  v2normalized.Normalize();

  return acos(v1normalized*v2normalized);
}

itk::Index<2> GetNextPixelAlongVector(const itk::Index<2>& pixel, const FloatVector2Type& vector)
{
  itk::Index<2> nextPixel = pixel + GetOffsetAlongVector(vector);

  return nextPixel;
}

itk::Offset<2> GetOffsetAlongVector(const FloatVector2Type& vector)
{
  FloatVector2Type normalizedVector = vector;
  normalizedVector.Normalize();

  itk::Offset<2> offset;
  offset[0] = Helpers::RoundAwayFromZero(normalizedVector[0]);
  offset[1] = Helpers::RoundAwayFromZero(normalizedVector[1]);

  return offset;
}


itk::Size<2> SizeFromRadius(const unsigned int radius)
{
  itk::Size<2> size;
  size.Fill(Helpers::SideLengthFromRadius(radius));

  return size;
}

void ITKImageToCIELabImage(const FloatVectorImageType* const image, FloatVectorImageType* const cielabImage)
{
  // Convert the first 3 channels to CIELab (this assumes the first 3 channels are RGB)
  RGBImageType::Pointer rgbImage = RGBImageType::New();
  VectorImageToRGBImage(image, rgbImage);
  RGBImageToCIELabImage(rgbImage, cielabImage);
}

void RGBImageToCIELabImage(RGBImageType* const rgbImage, FloatVectorImageType* const cielabImage)
{
  // The adaptor expects to be able to modify the image (even though we don't in this case),
  // so we cannot pass a const RGBImageType* const.
  // Convert RGB image to Lab color space
  typedef itk::Accessor::RGBToLabColorSpacePixelAccessor<unsigned char, float> RGBToLabColorSpaceAccessorType;
  typedef itk::ImageAdaptor<RGBImageType, RGBToLabColorSpaceAccessorType> RGBToLabAdaptorType;
  RGBToLabAdaptorType::Pointer rgbToLabAdaptor = RGBToLabAdaptorType::New();
  rgbToLabAdaptor->SetImage(rgbImage);

  // Disassembler
  typedef itk::VectorIndexSelectionCastImageFilter<RGBToLabAdaptorType, FloatScalarImageType> VectorIndexSelectionFilterType;
  VectorIndexSelectionFilterType::Pointer vectorIndexSelectionFilter = VectorIndexSelectionFilterType::New();
  vectorIndexSelectionFilter->SetInput(rgbToLabAdaptor);

  std::vector<FloatScalarImageType::Pointer> channels;

  // Reassembler
  typedef itk::ImageToVectorImageFilter<FloatScalarImageType> ReassemblerType;
  typename ReassemblerType::Pointer reassembler = ReassemblerType::New();

  for(unsigned int i = 0; i < 3; ++i)
    {
    channels.push_back(FloatScalarImageType::New());
    vectorIndexSelectionFilter->SetIndex(i);
    vectorIndexSelectionFilter->Update();
    DeepCopy(vectorIndexSelectionFilter->GetOutput(), channels[i].GetPointer());
    reassembler->SetNthInput(i, channels[i]);
    }

  reassembler->Update();

  // Copy to the output
  DeepCopy(reassembler->GetOutput(), cielabImage);
}

itk::ImageRegion<2> GetRegionInRadiusAroundPixel(const itk::Index<2>& pixel, const unsigned int radius)
{
  // This function returns a Region with the specified 'radius' centered at 'pixel'. By the definition of the radius of a square patch, the output region is (radius*2 + 1)x(radius*2 + 1).
  // Note: This region is not necessarily entirely inside the image!

  // The "index" is the lower left corner, so we need to subtract the radius from the center to obtain it
  itk::Index<2> lowerLeft;
  lowerLeft[0] = pixel[0] - radius;
  lowerLeft[1] = pixel[1] - radius;

  itk::ImageRegion<2> region;
  region.SetIndex(lowerLeft);
  itk::Size<2> size;
  size[0] = radius*2 + 1;
  size[1] = radius*2 + 1;
  region.SetSize(size);

  return region;
}


itk::Index<2> GetRegionCenter(const itk::ImageRegion<2>& region)
{
  // This assumes that the region is an odd size in both dimensions
  itk::Index<2> center;
  center[0] = region.GetIndex()[0] + region.GetSize()[0] / 2;
  center[1] = region.GetIndex()[1] + region.GetSize()[1] / 2;

  return center;
}

itk::Offset<2> OffsetFrom1DOffset(const itk::Offset<1>& offset1D, const unsigned int dimension)
{
  // Manually construct a 2D offset with 0 in all dimensions except the specified dimension
  itk::Offset<2> offset;
  offset.Fill(0);
  offset[dimension] = offset1D[0];

  return offset;
}

itk::ImageRegion<2> CropToRegion(const itk::ImageRegion<2>& inputRegion, const itk::ImageRegion<2>& targetRegion)
{
  // Returns the overlap of the inputRegion with the targetRegion.

  itk::ImageRegion<2> region = targetRegion;
  region.Crop(inputRegion);

  return region;
}

void OutputImageType(const itk::ImageBase<2>* const input)
{
  if(dynamic_cast<const FloatScalarImageType*>(input))
    {
    std::cout << "Image type FloatScalarImageType" << std::endl;
    }
  else if(dynamic_cast<const UnsignedCharScalarImageType*>(input))
    {
    std::cout << "Image type UnsignedCharScalarImageType" << std::endl;
    }
  else if(dynamic_cast<const FloatVectorImageType*>(input))
    {
    std::cout << "Image type FloatVectorImageType" << std::endl;
    }
  else
    {
    std::cout << "OutputImageType: Image is Invalid type!" << std::endl;
    }
}

// The return value MUST be a smart pointer
itk::ImageBase<2>::Pointer CreateImageWithSameType(const itk::ImageBase<2>* const input)
{
  itk::LightObject::Pointer objectCopyLight = input->CreateAnother();

  itk::ImageBase<2>::Pointer objectCopy = dynamic_cast<itk::ImageBase<2>*>(objectCopyLight.GetPointer());

  return objectCopy;
}

std::vector<itk::Index<2> > Get8NeighborsInRegion(const itk::ImageRegion<2>& region, const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighborsInRegion;

  std::vector<itk::Offset<2> > neighborOffsets = Get8NeighborOffsets();
  for(unsigned int i = 0; i < neighborOffsets.size(); ++i)
    {
    itk::Index<2> index = pixel + neighborOffsets[i];
    if(region.IsInside(index))
      {
      neighborsInRegion.push_back(index);
      }
    }
  return neighborsInRegion;
}

std::vector<itk::Offset<2> > Get8NeighborOffsets()
{
  std::vector<itk::Offset<2> > offsets;

  for(int i = -1; i <= 1; ++i)
    {
    for(int j = -1; j <= 1; ++j)
      {
      if(i == 0 && j == 0)
        {
        continue;
        }
      itk::Offset<2> offset;
      offset[0] = i;
      offset[1] = j;
      offsets.push_back(offset);
      }
    }
  return offsets;
}

// itk::VariableLengthVector<float> Average(const std::vector<itk::VariableLengthVector<float> >& v)
// {
//   // std::cout << "ITKHelpers::Average" << std::endl;
//   if(v.size() == 0)
//   {
//     throw std::runtime_error("Cannot average vector with size 0!");
//   }
//   itk::VariableLengthVector<float> vectorSum;
//   vectorSum.SetSize(v[0].GetSize());
//   vectorSum.Fill(0);
// 
//   for(unsigned int i = 0; i < v.size(); ++i)
//     {
//     //std::cout << "Average: Adding value " << v[i] << std::endl;
//     vectorSum += v[i];
//     //std::cout << "Average: Current vectorSum " << vectorSum << std::endl;
//     }
// 
//   itk::VariableLengthVector<float> averageVector;
//   averageVector.SetSize(v[0].GetSize());
//   averageVector = vectorSum / static_cast<float>(v.size());
// 
//   return averageVector;
// }

std::vector<itk::Index<2> > OffsetsToIndices(const std::vector<itk::Offset<2> >& offsets, const itk::Index<2>& index)
{
  std::vector<itk::Index<2> > indices;
  for(unsigned int i = 0; i < offsets.size(); ++i)
  {
    indices.push_back(index + offsets[i]);
  }
  return indices;
}

std::vector<itk::Index<2> > OffsetsToIndices(const std::vector<itk::Offset<2> >& offsets)
{
  std::vector<itk::Index<2> > indices;
  for(unsigned int i = 0; i < offsets.size(); ++i)
  {
    indices.push_back(CreateIndex(offsets[i]));
  }
  return indices;
}

std::vector<itk::Offset<2> > IndicesToOffsets(const std::vector<itk::Index<2> >& indices, const itk::Index<2>& index)
{
  std::vector<itk::Offset<2> > offsets;
  for(unsigned int i = 0; i < indices.size(); ++i)
  {
    offsets.push_back(indices[i] - index);
  }
  return offsets;
}

std::vector<itk::Index<2> > GetBoundaryPixels(const itk::ImageRegion<2>& region)
{
  std::vector<itk::Index<2> > boundaryPixels;

  for(unsigned int i = region.GetIndex()[0]; i < region.GetIndex()[0] + region.GetSize()[0]; ++i)
    {
    itk::Index<2> index;
    index[0] = i;
    index[1] = region.GetIndex()[1];
    boundaryPixels.push_back(index);

    index[0] = i;
    index[1] = region.GetIndex()[1] + region.GetSize()[1] - 1;
    boundaryPixels.push_back(index);
    }

  for(unsigned int j = region.GetIndex()[1]; j < region.GetIndex()[1] + region.GetSize()[1]; ++j)
    {
    itk::Index<2> index;
    index[0] = region.GetIndex()[0];
    index[1] = j;
    boundaryPixels.push_back(index);

    index[0] = region.GetIndex()[0] + region.GetSize()[0] - 1;
    index[1] = j;
    boundaryPixels.push_back(index);
    }

  return boundaryPixels;
}

itk::ImageRegion<2> CornerRegion(const itk::Size<2>& size)
{
  itk::ImageRegion<2> region(ZeroIndex(), size);
  return region;
}

void StackImages(const itk::VectorImage<float, 2>* const image1, const itk::VectorImage<float, 2>* const image2,
                 itk::VectorImage<float, 2>* const output)
{
  typedef itk::VectorImage<float, 2> VectorImageType;
  typedef itk::Image<float, 2> ScalarImageType;
  
  if(image1->GetLargestPossibleRegion() != image2->GetLargestPossibleRegion())
    {
    std::stringstream ss;
    ss << "StackImages: Images must be the same size!" << std::endl
       << "Image1 is " << image1->GetLargestPossibleRegion() << std::endl
       << "Image2 is " << image2->GetLargestPossibleRegion() << std::endl;
    throw std::runtime_error(ss.str());
    }

  std::cout << "StackImages: Image1 has " << image1->GetNumberOfComponentsPerPixel() << " components." << std::endl;
  std::cout << "StackImages: Image2 has " << image2->GetNumberOfComponentsPerPixel() << " components." << std::endl;

  // Create output image
  itk::ImageRegion<2> region = image1->GetLargestPossibleRegion();

  unsigned int newPixelLength = image1->GetNumberOfComponentsPerPixel() +
                                image2->GetNumberOfComponentsPerPixel();

  std::cout << "Output image has " << newPixelLength << " components." << std::endl;

  output->SetNumberOfComponentsPerPixel(newPixelLength);
  output->SetRegions(region);
  output->Allocate();

  for(unsigned int i = 0; i < image1->GetNumberOfComponentsPerPixel(); i++)
    {
    ScalarImageType::Pointer channel = ScalarImageType::New();
    channel->SetRegions(region);
    channel->Allocate();
  
    ExtractChannel(image1, i, channel.GetPointer());
    SetChannel(output, i, channel.GetPointer());
    }

  for(unsigned int i = 0; i < image2->GetNumberOfComponentsPerPixel(); i++)
    {
    ScalarImageType::Pointer channel = ScalarImageType::New();
    channel->SetRegions(region);
    channel->Allocate();

    ExtractChannel(image2, i, channel.GetPointer());
    SetChannel(output, image1->GetNumberOfComponentsPerPixel() + i, channel.GetPointer());
    }

}

std::vector<float> MinValues(const itk::VectorImage<float, 2>* const image, const itk::ImageRegion<2>& region)
{
  std::vector<float> mins(image->GetNumberOfComponentsPerPixel());

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    typedef itk::VectorImage<float, 2> VectorImageType;
    typedef itk::Image<float, 2> ScalarImageType;

    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(channel);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    mins[channel] = MinValue(indexSelectionFilter->GetOutput(), region);
  }

  return mins;
}


std::vector<float> MaxValues(const itk::VectorImage<float, 2>* const image, const itk::ImageRegion<2>& region)
{
  std::vector<float> maxs(image->GetNumberOfComponentsPerPixel());

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    typedef itk::VectorImage<float, 2> VectorImageType;
    typedef itk::Image<float, 2> ScalarImageType;

    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(channel);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    maxs[channel] = MaxValue(indexSelectionFilter->GetOutput(), region);
  }

  return maxs;
}

std::vector<float> MinValues(const itk::VectorImage<float, 2>* const image)
{
  std::vector<float> mins(image->GetNumberOfComponentsPerPixel());
  
  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    typedef itk::VectorImage<float, 2> VectorImageType;
    typedef itk::Image<float, 2> ScalarImageType;

    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(channel);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    mins[channel] = MinValue(indexSelectionFilter->GetOutput());
  }

  return mins;
}


std::vector<float> MaxValues(const itk::VectorImage<float, 2>* const image)
{
  std::vector<float> maxs(image->GetNumberOfComponentsPerPixel());

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    typedef itk::VectorImage<float, 2> VectorImageType;
    typedef itk::Image<float, 2> ScalarImageType;

    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType > IndexSelectionType;
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(channel);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->Update();

    maxs[channel] = MaxValue(indexSelectionFilter->GetOutput());
  }

  return maxs;
}


void Write2DVectorImage(const FloatVector2ImageType* const image, const std::string& filename)
{
  Write2DVectorRegion(image, image->GetLargestPossibleRegion(), filename);
}

void Write2DVectorRegion(const FloatVector2ImageType* const image, const itk::ImageRegion<2>& region, const std::string& filename)
{
  // This is a separate function than WriteRegion because Paraview requires vectors to be 3D to glyph them.

  typedef itk::RegionOfInterestImageFilter<FloatVector2ImageType, FloatVector2ImageType> RegionOfInterestImageFilterType;

  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  itk::Point<float, 2> origin;
  origin.Fill(0);
  regionOfInterestImageFilter->GetOutput()->SetOrigin(origin);

  FloatVector3ImageType::Pointer vectors3D = FloatVector3ImageType::New();
  vectors3D->SetRegions(regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion());
  vectors3D->Allocate();

  itk::ImageRegionConstIterator<FloatVector2ImageType> iterator(regionOfInterestImageFilter->GetOutput(), regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion());

  while(!iterator.IsAtEnd())
    {
    FloatVector2Type vec2d = iterator.Get();
    FloatVector3Type vec3d;
    vec3d[0] = vec2d[0];
    vec3d[1] = vec2d[1];
    vec3d[2] = 0;

    vectors3D->SetPixel(iterator.GetIndex(), vec3d);
    ++iterator;
    }

  //std::cout << "regionOfInterestImageFilter " << regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;

  itk::ImageFileWriter<FloatVector3ImageType>::Pointer writer = itk::ImageFileWriter<FloatVector3ImageType>::New();
  writer->SetFileName(filename);
  writer->SetInput(vectors3D);
  writer->Update();
}

std::vector<itk::Index<2> > DilatePixelList(const std::vector<itk::Index<2> >& pixelList,
                                            const itk::ImageRegion<2>& region, const unsigned int radius)
{
  //std::cout << "DilatePixelList: input has " << pixelList.size() << " pixels." << std::endl;
  // Construct an image of the pixels in the list
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  typedef std::vector<itk::Index<2> > PixelVectorType;

  for(PixelVectorType::const_iterator iter = pixelList.begin(); iter != pixelList.end(); ++iter)
  {
    // Note, this must be 255, not just any non-zero number, for BinaryDilateImageFilter to work properly.
    image->SetPixel(*iter, 255);
  }

  //WriteImage(image.GetPointer(), "beforeDilation.png");

  // Dilate the image
  typedef itk::BinaryBallStructuringElement<ImageType::PixelType,2> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(image);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  //WriteImage(dilateFilter->GetOutput(), "afterDilation.png");

  PixelVectorType dilatedPixelList;

  itk::ImageRegionConstIteratorWithIndex<ImageType> imageIterator(dilateFilter->GetOutput(),
                                                         dilateFilter->GetOutput()->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get())
      {
      dilatedPixelList.push_back(imageIterator.GetIndex());
      }
    ++imageIterator;
    }

  //std::cout << "DilatePixelList: output has " << dilatedPixelList.size() << " pixels." << std::endl;
  return dilatedPixelList;
}


void IndicesToBinaryImage(const std::vector<itk::Index<2> >& indices, UnsignedCharScalarImageType* const image)
{
  // The Regions of the 'image' must be set before calling this function
  //std::cout << "Setting " << indices.size() << " points to non-zero." << std::endl;

  image->Allocate();
  image->FillBuffer(0);

  // Set the pixels of indices in list to 255
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    image->SetPixel(indices[i], 255);
    }
}

itk::Size<2> Get1x1Radius()
{
  itk::Size<2> radius;
  radius.Fill(1);
  return radius;
}


std::vector<itk::Index<2> > Get4NeighborIndicesInsideRegion(const itk::Index<2>& pixel,
                                                            const itk::ImageRegion<2>& region)
{
  std::vector<itk::Index<2> > indices;

  itk::Offset<2> offset;
  offset[0] = -1;
  offset[1] = 0;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 1;
  offset[1] = 0;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 0;
  offset[1] = -1;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  offset[0] = 0;
  offset[1] = 1;

  if(region.IsInside(pixel + offset))
  {
    indices.push_back(pixel + offset);
  }

  return indices;
}


void RGBImageToVectorImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                           itk::VectorImage<float, 2>* const outputImage)
{
  outputImage->SetRegions(image->GetLargestPossibleRegion());
  outputImage->SetNumberOfComponentsPerPixel(3);
  outputImage->Allocate();

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  typedef itk::VectorImage<float, 2> VectorImageType;

  itk::ImageRegionConstIteratorWithIndex<RGBImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Can't get a reference of a pixel from a VectorImage apparently?
    //VectorImageType::PixelType& outputPixel = outputImage->GetPixel(imageIterator.GetIndex());
    VectorImageType::PixelType outputPixel = outputImage->GetPixel(imageIterator.GetIndex());
    outputPixel[0] = imageIterator.Get().GetRed();
    outputPixel[1] = imageIterator.Get().GetGreen();
    outputPixel[2] = imageIterator.Get().GetBlue();
    outputImage->SetPixel(imageIterator.GetIndex(), outputPixel);

    ++imageIterator;
    }
}


std::pair<itk::Index<2>, itk::Index<2> > IntersectLineWithHole(std::vector<itk::Index<2> > line, UnsignedCharScalarImageType::Pointer mask, bool &hasInteriorLine)
{
  // We consider the hole to be non-zero pixels of the mask. We want to find where the line enters the mask, and where it leaves the mask.
  // This function assumes that the line starts outside the mask. Nothing is assumed about where the line ends (if it ends inside the mask, then there is no interior line).
  // 'line' is an ordered vector of indices.
  // We assume the hole is convex. Nothing will break if it is not, but the line that is computed goes "through" the mask, but may
  // actually not be entirely contained within the hole if the hole is not convex.

  std::pair<itk::Index<2>, itk::Index<2> > interiorLine; // (start pixel, end pixel)

  unsigned int startPoints = 0;
  unsigned int endPoints = 0;

  // Loop over the pixels in the line. If one of them is outside the mask and its neighbor is inside the mask, this is an intersection.
  for(unsigned int i = 0; i < line.size() - 1; i++) // loop to one before the end because we use the current and current+1 in the loop
    {
    if(mask->GetPixel(line[i]) == 0 && mask->GetPixel(line[i+1]) != 0) // Found entry point
      {
      interiorLine.first = line[i]; // We want to save the outside/valid/non-hole point. This is the first point (i) in the 'exit' case.
      startPoints++;
      }

    if(mask->GetPixel(line[i]) != 0 && mask->GetPixel(line[i+1]) == 0) // Found exit point
      {
      interiorLine.second = line[i+1]; // We want to save the outside/valid/non-hole point. This is the second point (i+1) in the 'exit' case.
      endPoints++;
      }
    }

  // If there is exactly one entry and exactly 1 exit point, the interior line is well defined
  if(startPoints == 1 && endPoints == 1)
    {
    hasInteriorLine = true;
    }
  else
    {
    hasInteriorLine = false;
    }

  // This is only valid if hasInteriorLine is true
  return interiorLine;
}

itk::ImageRegion<2> GetInternalRegion(const itk::ImageRegion<2>& wholeRegion, const unsigned int patchRadius)
{
  unsigned int width = wholeRegion.GetSize()[0];
  unsigned int height = wholeRegion.GetSize()[1];

  itk::Index<2> regionCorner = {{patchRadius, patchRadius}};
  itk::Size<2> regionSize = {{width - 2*patchRadius, height - 2*patchRadius}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  return region;
}

std::vector<itk::ImageRegion<2> > GetAllPatchesContainingPixel(const itk::Index<2>& pixel,
                                                               const unsigned int patchRadius,
                                                               const itk::ImageRegion<2>& imageRegion)
{
  typedef itk::Image<float,2> DummyImageType;
  DummyImageType::Pointer dummyImage = DummyImageType::New();
  dummyImage->SetRegions(imageRegion);
  //dummyImage->Allocate(); // Do we actually need this to iterate over the image?

  // This region includes all patch centers
  itk::Index<2> possibleRegionCorner = {{pixel[0] - patchRadius, pixel[1] - patchRadius}};
  itk::Size<2> possibleRegionSize = {{patchRadius*2 + 1, patchRadius*2 + 1}};
  itk::ImageRegion<2> possibleRegion(possibleRegionCorner, possibleRegionSize);

  // Don't include patch centers that fall outside the image
  possibleRegion.Crop(imageRegion);

  itk::ImageRegionIteratorWithIndex<DummyImageType> imageIterator(dummyImage, possibleRegion);

  std::vector<itk::ImageRegion<2> > regions;

  // Each pixel in this loop is a potential patch center
  while(!imageIterator.IsAtEnd())
  {
    itk::ImageRegion<2> region = GetRegionInRadiusAroundPixel(imageIterator.GetIndex(), patchRadius);
    if(imageRegion.IsInside(region))
    {
      regions.push_back(region);
    }
    ++imageIterator;
  }

  return regions;
}


unsigned int ClosestPoint(const std::vector<itk::CovariantVector<float, 3> >& vec, const itk::CovariantVector<float, 3>& value)
{
  typedef itk::CovariantVector<float, 3> PointType;

  std::vector<float> distances(vec.size());
  for(size_t i = 0; i < vec.size(); ++i)
  {
    distances[i] = (vec[i] - value).GetSquaredNorm();
  }

  return Helpers::argmin(distances);
}

itk::Size<2> MakeSizeEven(const itk::Size<2>& inputSize)
{
  itk::Size<2> outputSize = inputSize;
  
  if(Helpers::IsOdd(inputSize[0]))
  {
    outputSize[0]--;
  }

  if(Helpers::IsOdd(inputSize[1]))
  {
    outputSize[1]--;
  }

  return outputSize;
}

} // end namespace
