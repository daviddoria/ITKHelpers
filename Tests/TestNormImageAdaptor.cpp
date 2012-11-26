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

// Custom
#include "itkNormImageAdaptor.h"

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"

static void TestGetPixel();
static void TestAdaptorSize();
static void TestChangeImage();

int main()
{

  TestGetPixel();
  TestAdaptorSize();
  TestChangeImage();
  return 0;
}

void TestGetPixel()
{

  typedef itk::Image<itk::CovariantVector<float, 3> > ImageType;

  itk::Index<2> corner = {{0,0}};

  itk::Size<2> size = {{10,10}};

  itk::ImageRegion<2> region(corner, size);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  ImageType::PixelType pixel = 2; // (2,2,2)

  image->SetPixel(corner, pixel);

  typedef itk::NormImageAdaptor<ImageType> ImageAdaptorType;

  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(image);

  std::cout << adaptor->GetPixel(corner) << std::endl;

}

void TestAdaptorSize()
{

  typedef itk::Image<itk::CovariantVector<float, 3> > ImageType;

  itk::Index<2> corner = {{0,0}};

  itk::Size<2> size = {{10,10}};

  itk::ImageRegion<2> region(corner, size);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  typedef itk::NormImageAdaptor<ImageType> ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(image);

  std::cout << "image region: " << image->GetLargestPossibleRegion() << std::endl;
  std::cout << "adaptor region: " << adaptor->GetLargestPossibleRegion() << std::endl;

}

void TestChangeImage()
{
  std::cout << "TestChangeImage()" << std::endl;
  typedef itk::Image<itk::CovariantVector<float, 3> > ImageType;

  itk::Index<2> corner = {{0,0}};

  itk::Size<2> originalSize = {{0,0}};

  itk::ImageRegion<2> region(corner, originalSize);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  typedef itk::NormImageAdaptor<ImageType> ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(image);

  std::cout << "image region: " << image->GetLargestPossibleRegion() << std::endl;
  std::cout << "adaptor region: " << adaptor->GetLargestPossibleRegion() << std::endl;

  itk::Size<2> newSize = {{10,10}};
  region.SetSize(newSize);
  image->SetRegions(region);

  std::cout << "image region: " << image->GetLargestPossibleRegion() << std::endl;
  std::cout << "adaptor region: " << adaptor->GetLargestPossibleRegion() << std::endl;

}
