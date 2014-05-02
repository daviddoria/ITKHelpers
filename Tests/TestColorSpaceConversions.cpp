/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
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

// ITKHelpers
#include "itkRGBToHSVColorSpacePixelAccessor.h"
#include "ITKHelpers.h"

// ITK
#include "itkImage.h"
#include "itkImageAdaptor.h"
#include "itkVectorImage.h"

static bool TestAdaptor();
static bool TestConversion_Range();
static bool TestConversion_Instantiations();

int main()
{
  bool allPass = true;

  allPass &= TestAdaptor();

  allPass &= TestConversion_Range();

  allPass &= TestConversion_Instantiations();

  if(allPass)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}

bool TestAdaptor()
{
  typedef itk::CovariantVector<unsigned char, 3> RGBPixelType;

  typedef itk::Image<RGBPixelType, 2> ImageType;

  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};

  itk::ImageRegion<2> region(corner,size);

  image->SetRegions(region);
  image->Allocate();

  RGBPixelType red;
  red[0] = 255;
  red[1] = 0;
  red[2] = 0;
  image->SetPixel(corner, red);

  typedef itk::Accessor::RGBToHSVColorSpacePixelAccessor<unsigned char, float, RGBPixelType>
      RGBToHSVColorSpaceAccessorType;
  typedef itk::ImageAdaptor<ImageType, RGBToHSVColorSpaceAccessorType> RGBToHSVImageAdaptorType;
  RGBToHSVImageAdaptorType::Pointer rgbToHSVAdaptor = RGBToHSVImageAdaptorType::New();
  rgbToHSVAdaptor->SetImage(image);

  auto convertedPixel = rgbToHSVAdaptor->GetPixel(corner);
  decltype(convertedPixel) correctConversion;
  correctConversion[0] = 0;
  correctConversion[1] = 1;
  correctConversion[2] = 1;
  if(!Helpers::FuzzyCompare(convertedPixel, correctConversion))
  {
    return false;
  }

  return true;
}

bool TestConversion_Instantiations()
{
  std::cout << "TestConversion_Instantiations()" << std::endl;

  // From Image<CovariantVector> to VectorImage
  {
  typedef itk::CovariantVector<unsigned char, 3> RGBPixelType;

  typedef itk::Image<RGBPixelType, 2> ImageType;

  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};

  itk::ImageRegion<2> region(corner,size);

  image->SetRegions(region);
  image->Allocate();

  RGBPixelType red;
  red[0] = 255;
  red[1] = 0;
  red[2] = 0;
  image->SetPixel(corner, red);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;
  //typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;
  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());

  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  // From Image<RGBPixel> to Image<CovariantVector>
  {
  typedef itk::RGBPixel<unsigned char> RGBPixelType;

  typedef itk::Image<RGBPixelType, 2> ImageType;

  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};

  itk::ImageRegion<2> region(corner,size);

  image->SetRegions(region);
  image->Allocate();

  RGBPixelType red;
  red[0] = 255;
  red[1] = 0;
  red[2] = 0;
  image->SetPixel(corner, red);

  typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;
  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());

  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  return true;
}


bool TestConversion_Range()
{
  std::cout << "TestConversion_Range()" << std::endl;

  typedef itk::CovariantVector<unsigned char, 3> RGBPixelType;

  typedef itk::Image<RGBPixelType, 2> ImageType;

  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};

  itk::ImageRegion<2> region(corner,size);

  image->SetRegions(region);
  image->Allocate();

  RGBPixelType color;

  // Red
  {
  std::cout << "Red: ";
  color[0] = 255;
  color[1] = 0;
  color[2] = 0;
  image->SetPixel(corner, color);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;

  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());
  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  // Green
  {
  std::cout << "Green: ";
  color[0] = 0;
  color[1] = 255;
  color[2] = 0;
  image->SetPixel(corner, color);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;

  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());
  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  // Blue
  {
  std::cout << "Blue: ";
  color[0] = 0;
  color[1] = 0;
  color[2] = 255;
  image->SetPixel(corner, color);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;

  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());
  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  // Black
  {
  std::cout << "Black: ";
  color[0] = 0;
  color[1] = 0;
  color[2] = 0;
  image->SetPixel(corner, color);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;

  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());
  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  // White
  {
  std::cout << "White: ";
  color[0] = 255;
  color[1] = 255;
  color[2] = 255;
  image->SetPixel(corner, color);

  typedef itk::VectorImage<float, 2> FloatVectorImageType;

  FloatVectorImageType::Pointer hsvImage = FloatVectorImageType::New();
  ITKHelpers::ITKImageToHSVImage(image.GetPointer(), hsvImage.GetPointer());
  std::cout << hsvImage->GetPixel(corner) << std::endl;
  }

  return true;
}
