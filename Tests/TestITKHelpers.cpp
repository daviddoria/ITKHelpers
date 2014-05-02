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

#include "ITKHelpers.h"

static bool TestGetClosedContourOrdering();

static bool TestGetOpenContourOrdering();

static bool TestDrawRectangle();

static bool TestIsClosedLoop();

static bool TestRandomImage();

static bool TestBlurAllChannelsScalar();

static bool TestBlurAllChannelsVector();

static bool TestBilateralFilterAllChannels();

static bool TestHistogramOfGradients();

static bool TestExtractChannel();
static bool TestExtractChannels();

static bool TestSumOfComponentMagnitudes();

static bool TestClosestValueIndex();

bool TestDeepCopyFloatScalar();
bool TestDeepCopyUnsignedCharScalar();
bool TestDeepCopyFloatVector();
bool TestDeepCopyUnsignedCharVector();

static bool TestUpsample();
static bool TestDownsample();

static bool TestGetAllPatchesContainingPixel();

static bool TestBreadthFirstOrderingNonZeroPixels();

static bool TestGetBoundaryPixels();

static bool TestDivideRegion();

static bool TestMinOfIndex();

static bool TestMinOfAllIndices();

static bool TestComputeGradientsInRegion();

static bool TestCreateLuminanceImage();

static bool TestHasBracketOperator();

template<typename TImage>
static bool TestHasBracketOperator_ConstTemplate(const TImage* const image);

template<typename TImage>
static bool TestHasBracketOperator_Template(TImage* image);

template<typename TImage>
static bool TestHasBracketOperator_Template2();

int main(int, char **)
{
  bool allPass = true;

  allPass &= TestRandomImage();

  allPass &= TestGetClosedContourOrdering();

  allPass &= TestGetOpenContourOrdering();

  allPass &= TestDrawRectangle();

  allPass &= TestIsClosedLoop();

  allPass &= TestClosestValueIndex();

  allPass &= TestBreadthFirstOrderingNonZeroPixels();

  allPass &= TestBlurAllChannelsVector();
  allPass &= TestBlurAllChannelsScalar();
  allPass &= TestBilateralFilterAllChannels();

  allPass &= TestHistogramOfGradients();

  allPass &= TestExtractChannel();
  allPass &= TestExtractChannels();

  allPass &= TestSumOfComponentMagnitudes();

  allPass &= TestGetAllPatchesContainingPixel();

  allPass &= TestDownsample();
  allPass &= TestUpsample();

  allPass &= TestDeepCopyFloatScalar();
  allPass &= TestDeepCopyUnsignedCharScalar();
  allPass &= TestDeepCopyFloatVector();
  allPass &= TestDeepCopyUnsignedCharVector();

  allPass &= TestGetBoundaryPixels();

  allPass &= TestDivideRegion();

  allPass &= TestMinOfIndex();

  allPass &= TestMinOfAllIndices();

  allPass &= TestComputeGradientsInRegion();

  allPass &= TestCreateLuminanceImage();

  allPass &= TestHasBracketOperator();

  if(allPass)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}


bool TestBlurAllChannelsScalar()
{
  // Create an image
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  ITKHelpers::RandomImage(image.GetPointer());
  ITKHelpers::WriteImage(image.GetPointer(), "image.png");

  // Some blur
  {
  ImageType::Pointer blurred = ImageType::New();

  float sigma = 2.0f;
  std::cout << "Blurring with sigma = " << sigma << std::endl;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);

  ITKHelpers::WriteImage(blurred.GetPointer(), "blurred_2.png");
  }

  // No blur
  {
  ImageType::Pointer blurred = ImageType::New();

  //float sigma = 0.0f; // Don't do this! Resulting image is black
  float sigma = 1.0f;
  std::cout << "Blurring with sigma = " << sigma << std::endl;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
  ITKHelpers::WriteImage(blurred.GetPointer(), "blurred_1.png");
  }

  return true;
}

bool TestBlurAllChannelsVector()
{
  // Create an image
  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(3);
  image->Allocate();

  ITKHelpers::RandomImage(image.GetPointer());
  ITKHelpers::WriteRGBImage(image.GetPointer(), "image.png");

  // Some blur
  {
  ImageType::Pointer blurred = ImageType::New();

  float sigma = 2.0f;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);

  ITKHelpers::WriteRGBImage(blurred.GetPointer(), "blurred_2.png");
  }

  // No blur
  {
  ImageType::Pointer blurred = ImageType::New();

  // float sigma = 0.0f; // Don't do this! Resulting image is black
  float sigma = 1.0f;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
  ITKHelpers::WriteRGBImage(blurred.GetPointer(), "blurred_0.png");
  }

  return true;
}

bool TestBilateralFilterAllChannels()
{
  {
  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer blurred = ImageType::New();

  float domainSigma = 2.0f;
  float rangeSigma = 2.0f;
  ITKHelpers::BilateralFilterAllChannels(image.GetPointer(), blurred.GetPointer(), domainSigma, rangeSigma);
  }

  // This does not work
//   {
//   typedef itk::Image<itk::CovariantVector<float, 3>, 2> ImageType;
//   ImageType::Pointer image = ImageType::New();
//   ImageType::Pointer blurred = ImageType::New();
//
//   float sigma = 2.0f;
//   ITKHelpers::AnisotropicBlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
//   }

  return true;
}


bool TestGetAllPatchesContainingPixel()
{
  typedef itk::Image<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> imageRegion(corner,size);

  image->SetRegions(imageRegion);
  image->Allocate();
  image->FillBuffer(0.0f);

  itk::Index<2> queryPixel = {{5,5}};

  std::vector<itk::ImageRegion<2> > allPatches =
        ITKHelpers::GetAllPatchesContainingPixel(queryPixel, 1, imageRegion);

  std::cout << "Number of patches " << allPatches.size() << std::endl;

  return true;
}

bool TestDeepCopyFloatScalar()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> region(corner, size);

  typedef itk::Image<float, 2> ImageType;
  ImageType::Pointer imageIn = ImageType::New();
  imageIn->SetRegions(region);
  imageIn->Allocate();

  ImageType::Pointer imageOut = ImageType::New();
  ITKHelpers::DeepCopy(imageIn.GetPointer(), imageOut.GetPointer());

  return true;
}

bool TestDeepCopyUnsignedCharScalar()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> region(corner, size);

  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer imageIn = ImageType::New();
  imageIn->SetRegions(region);
  imageIn->Allocate();

  ImageType::Pointer imageOut = ImageType::New();
  ITKHelpers::DeepCopy(imageIn.GetPointer(), imageOut.GetPointer());

  return true;
}

bool TestDeepCopyFloatVector()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> region(corner, size);

  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer imageIn = ImageType::New();
  imageIn->SetRegions(region);
  imageIn->SetNumberOfComponentsPerPixel(3);
  imageIn->Allocate();

  ImageType::Pointer imageOut = ImageType::New();
  ITKHelpers::DeepCopy(imageIn.GetPointer(), imageOut.GetPointer());

  return true;
}

bool TestDeepCopyUnsignedCharVector()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> region(corner, size);

  typedef itk::VectorImage<unsigned char, 2> ImageType;
  ImageType::Pointer imageIn = ImageType::New();
  imageIn->SetRegions(region);
  imageIn->SetNumberOfComponentsPerPixel(3);
  imageIn->Allocate();

  ImageType::Pointer imageOut = ImageType::New();
  ITKHelpers::DeepCopy(imageIn.GetPointer(), imageOut.GetPointer());

  return true;
}

bool TestDownsample()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  ImageType::Pointer downsampled = ImageType::New();

  ITKHelpers::Downsample(image.GetPointer(), 2, downsampled.GetPointer());

  std::cout << "TestDownsample() output size: "
            << downsampled->GetLargestPossibleRegion().GetSize() << std::endl;

  return true;
}

bool TestUpsample()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  ImageType::Pointer upsampled = ImageType::New();

  ITKHelpers::Upsample(image.GetPointer(), 2, upsampled.GetPointer());

  std::cout << "TestUpsample() output size: "
            << upsampled->GetLargestPossibleRegion().GetSize() << std::endl;

  return true;
}

bool TestClosestValueIndex()
{
  typedef itk::CovariantVector<float, 3> PointType;
  std::vector<PointType> vec;

  PointType a;
  a.Fill(1.1);
  vec.push_back(a);

  PointType b;
  b.Fill(2.1);
  vec.push_back(b);

  PointType c;
  c.Fill(3.1);
  vec.push_back(c);

  PointType query;
  query.Fill(1.2);

  unsigned int closestId = ITKHelpers::ClosestValueIndex(vec, query);

  std::cout << "Closest point to " << query << " is " << vec[closestId] << std::endl;

  return true;
}

bool TestSumOfComponentMagnitudes()
{
  itk::VariableLengthVector<float> a(2);
  a.Fill(3.0);

  float sum = ITKHelpers::SumOfComponentMagnitudes(a);

  std::cout << "Sum of " << a << " is " << sum << std::endl;

  return true;
}

bool TestExtractChannel()
{
  // VectorImage
  {
  typedef itk::VectorImage<float, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(2);
  image->Allocate();

  typedef itk::Image<float, 2> FloatScalarImageType;
  FloatScalarImageType::Pointer floatScalarImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, floatScalarImage.GetPointer());

  typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;
  UnsignedCharScalarImageType::Pointer unsignedCharScalarImage = UnsignedCharScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, unsignedCharScalarImage.GetPointer());
  }

  // VectorImage different output type
  {
  typedef itk::VectorImage<float, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(2);
  image->Allocate();

  typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;
  UnsignedCharScalarImageType::Pointer unsignedCharScalarImage = UnsignedCharScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, unsignedCharScalarImage.GetPointer());
  }

  // Scalar Image
  {
  typedef itk::Image<float, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  typedef itk::Image<float, 2> FloatScalarImageType;
  FloatScalarImageType::Pointer floatScalarImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, floatScalarImage.GetPointer());

  typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;
  UnsignedCharScalarImageType::Pointer unsignedCharScalarImage = UnsignedCharScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, unsignedCharScalarImage.GetPointer());
  }

  // Image<CovariantVector>
  {
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  typedef itk::Image<float, 2> FloatScalarImageType;
  FloatScalarImageType::Pointer floatScalarImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, floatScalarImage.GetPointer());
  }

  // Image<Vector>
  {
  typedef itk::Image<itk::Vector<float, 3>, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  typedef itk::Image<float, 2> FloatScalarImageType;
  FloatScalarImageType::Pointer floatScalarImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannel(image.GetPointer(), 0, floatScalarImage.GetPointer());
  }

  return true;
}

bool TestExtractChannels()
{
  typedef itk::VectorImage<float, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(3);
  image->Allocate();

  // Extract the first two channels
  std::vector<unsigned int> channels;
  channels.push_back(0);
  channels.push_back(1);

  typedef itk::VectorImage<float, 2> FloatScalarImageType;
  FloatScalarImageType::Pointer floatScalarImage = FloatScalarImageType::New();
  ITKHelpers::ExtractChannels(image.GetPointer(), channels, floatScalarImage.GetPointer());

  typedef itk::VectorImage<unsigned char, 2> UnsignedCharScalarImageType;
  UnsignedCharScalarImageType::Pointer unsignedCharScalarImage = UnsignedCharScalarImageType::New();
  ITKHelpers::ExtractChannels(image.GetPointer(), channels, unsignedCharScalarImage.GetPointer());

  return true;
}

bool TestHistogramOfGradients()
{
  typedef itk::Image<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  std::vector<float> histogram =
         ITKHelpers::HistogramOfGradients(image.GetPointer(),
                                          image->GetLargestPossibleRegion(), 10);

  for(unsigned int i = 0; i < histogram.size(); ++i)
  {
    std::cout << histogram[i] << " ";
  }
  std::cout << std::endl;

  return true;
}

bool TestBreadthFirstOrderingNonZeroPixels()
{
  std::cout << "TestBreadthFirstOrderingNonZeroPixels()" << std::endl;

  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  for(itk::Index<2>::IndexValueType i = 20; i < 30; ++i)
  {
    itk::Index<2> pixel = {{i, 50}};
    image->SetPixel(pixel, 255);
  }

  itk::Index<2> start = {{25, 50}};
  std::vector<itk::Index<2> > breadthFirstOrdering =
      ITKHelpers::BreadthFirstOrderingNonZeroPixels(image.GetPointer(), start);

  for(unsigned int i = 0; i < breadthFirstOrdering.size(); ++i)
  {
    //std::cout << breadthFirstOrdering[i] << " ";
    std::cout << breadthFirstOrdering[i] << std::endl;
  }

  std::cout << std::endl;

  return true;
}

bool TestIsClosedLoop()
{
  std::cout << "TestIsClosedLoop()" << std::endl;

  // Open loop
  {
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  for(itk::Index<2>::IndexValueType i = 20; i < 30; ++i)
  {
    itk::Index<2> pixel = {{i, 50}};
    image->SetPixel(pixel, 255);
  }

  itk::Index<2> start = {{25, 50}};

  bool isClosedLoop = ITKHelpers::IsClosedLoop(image.GetPointer(), start);

  std::cout << "Is closed loop? " << isClosedLoop << std::endl;
  }

  // Closed loop
  {
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  itk::Index<2> corner0 = {{10,10}};
  itk::Index<2> corner1 = {{30,30}};
  ITKHelpers::DrawRectangle(image.GetPointer(), 255,
                            corner0, corner1);

  itk::Index<2> start = {{10, 10}};

  bool isClosedLoop = ITKHelpers::IsClosedLoop(image.GetPointer(), start);

  std::cout << "Is closed loop? " << isClosedLoop << std::endl;
  }

  return true;
}

bool TestDrawRectangle()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  itk::Index<2> corner0 = {{10,10}};
  itk::Index<2> corner1 = {{30,30}};
  ITKHelpers::DrawRectangle(image.GetPointer(), 255,
                            corner0, corner1);

  ITKHelpers::WriteImage(image.GetPointer(), "rectangle.png");

  return true;
}

bool TestGetOpenContourOrdering()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  for(itk::Index<2>::IndexValueType i = 20; i < 30; ++i)
  {
    itk::Index<2> pixel = {{i, 50}};
    image->SetPixel(pixel, 255);
  }

  itk::Index<2> start = {{25, 50}};

  std::vector<itk::Index<2> > contourOrdering = ITKHelpers::GetOpenContourOrdering(image.GetPointer(), start);

  for(unsigned int i = 0; i < contourOrdering.size(); ++i)
  {
    //std::cout << breadthFirstOrdering[i] << " ";
    std::cout << contourOrdering[i] << std::endl;
  }

  std::cout << std::endl;

  return true;
}

bool TestGetClosedContourOrdering()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  itk::Index<2> corner0 = {{10,10}};
  itk::Index<2> corner1 = {{30,30}};
  ITKHelpers::DrawRectangle(image.GetPointer(), 255,
                            corner0, corner1);

  itk::Index<2> start = {{10, 10}};

  std::vector<itk::Index<2> > ordering = ITKHelpers::GetClosedContourOrdering(image.GetPointer(), start);

  for(unsigned int i = 0; i < ordering.size(); ++i)
  {
    std::cout << ordering[i] << std::endl;
  }

  std::cout << std::endl;

  return true;
}

bool TestRandomImage()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  ITKHelpers::RandomImage(image.GetPointer());
  ITKHelpers::WriteImage(image.GetPointer(), "random.png");

  return true;
}

bool TestGetBoundaryPixels()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{3,3}};
  itk::ImageRegion<2> region(corner,size);

  unsigned int thickness = 1;
  
  std::vector<itk::Index<2> > boundaryPixels = ITKHelpers::GetBoundaryPixels(region, thickness);

  for(size_t i = 0; i < boundaryPixels.size(); ++i)
  {
    std::cout << boundaryPixels[i] << std::endl;
  }

  return true;
}

bool TestDivideRegion()
{
  // Test an evenly divisible region
  {
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{10,10}};
  itk::ImageRegion<2> region(corner, size);
  const unsigned int divisionsPerDimension = 2;
  std::vector<itk::ImageRegion<2> > subregions = ITKHelpers::DivideRegion(region, divisionsPerDimension);

  std::cout << "There are " << subregions.size() << " subregions" << std::endl;
  for(size_t i = 0; i < subregions.size(); ++i)
  {
    std::cout << "Subregion " << i << " : " << subregions[i] << std::endl;
  }
  }

  // Test a non-evenly divisible region
  {
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{11,11}};
  itk::ImageRegion<2> region(corner, size);
  const unsigned int divisionsPerDimension = 2;
  std::vector<itk::ImageRegion<2> > subregions = ITKHelpers::DivideRegion(region, divisionsPerDimension);

  std::cout << "There are " << subregions.size() << " subregions" << std::endl;
  for(size_t i = 0; i < subregions.size(); ++i)
  {
    std::cout << "Subregion " << i << " : " << subregions[i] << std::endl;
  }
  }

  return true;
}

bool TestMinOfIndex()
{
  typedef itk::CovariantVector<int, 3> CovariantVectorType;
  std::vector<CovariantVectorType> vectorOfVectors;
  for(unsigned int i = 4; i < 10; ++i)
  {
    CovariantVectorType v;
    v[0] = i;
    v[1] = i;
    v[2] = i;
    vectorOfVectors.push_back(v);
  }

  int minComponent2 = Helpers::MinOfIndex(vectorOfVectors, 2);

  std::cout << "minComponent2: " << minComponent2 << std::endl;

  return true;
}

bool TestMinOfAllIndices()
{
  typedef itk::CovariantVector<int, 3> VectorType;
  std::vector<VectorType> vectorOfVectors;
  for(unsigned int i = 4; i < 10; ++i)
  {
    VectorType v;
    v[0] = i;
    v[1] = i;
    v[2] = i;
    vectorOfVectors.push_back(v);
  }

  VectorType minComponents;

  Helpers::MinOfAllIndices(vectorOfVectors, minComponents);

  std::cout << "minComponents: ";
  for(size_t i = 0; i < minComponents.GetNumberOfComponents(); ++i)
  {
    std::cout << minComponents[i] << std::endl;
  }

  return true;
}

bool TestComputeGradientsInRegion()
{
  itk::Index<2> imageCorner = {{0,0}};
  itk::Size<2> imageSize = {{100,100}};
  itk::ImageRegion<2> imageRegion(imageCorner, imageSize);

  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(imageRegion);
  image->Allocate();

  itk::ImageRegionIterator<ImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
  {
    imageIterator.Set(rand() % 255);

    ++imageIterator;
  }

  ITKHelpers::WriteImage(image.GetPointer(), "Image.mha");

  itk::Index<2> regionCorner = {{50,50}};
  itk::Size<2> regionSize = {{10,10}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  typedef itk::Image<itk::CovariantVector<float, 2>, 2> GradientImageType;
  GradientImageType::Pointer gradientImage = GradientImageType::New();

  std::cout << "Computing gradients..." << std::endl;
  ITKHelpers::ComputeGradientsInRegion(image.GetPointer(), region, gradientImage.GetPointer());

  std::cout << "Writing..." << std::endl;
  ITKHelpers::WriteImage(gradientImage.GetPointer(), "GradientImage.mha");

  return true;
}

bool TestCreateLuminanceImage()
{
  // From RGB image
  {
  itk::Index<2> imageCorner = {{0,0}};
  itk::Size<2> imageSize = {{100,100}};
  itk::ImageRegion<2> imageRegion(imageCorner, imageSize);

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  RGBImageType::Pointer rgbImage = RGBImageType::New();
  rgbImage->SetRegions(imageRegion);
  rgbImage->Allocate();

  typedef itk::Image<float, 2> LuminanceImageType;
  LuminanceImageType::Pointer luminanceImage = LuminanceImageType::New();

  ITKHelpers::CreateLuminanceImage(rgbImage.GetPointer(), luminanceImage.GetPointer());
  }

  // From Vector image
  {
  itk::Index<2> imageCorner = {{0,0}};
  itk::Size<2> imageSize = {{100,100}};
  itk::ImageRegion<2> imageRegion(imageCorner, imageSize);

  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> VectorImageType;
  VectorImageType::Pointer vectorImage = VectorImageType::New();
  vectorImage->SetRegions(imageRegion);
  vectorImage->Allocate();

  typedef itk::Image<float, 2> LuminanceImageType;
  LuminanceImageType::Pointer luminanceImage = LuminanceImageType::New();

  ITKHelpers::CreateLuminanceImage(vectorImage.GetPointer(), luminanceImage.GetPointer());
  }

  return true;
}

bool TestHasBracketOperator()
{
  {
  typedef itk::CovariantVector<int, 3> IntVectorType;

  static_assert(Helpers::HasBracketOperator<IntVectorType>::value,
              "TestHasBracketOperator for CovariantVector<int, 3> failed!");
  }

  {
  typedef itk::CovariantVector<unsigned char, 3> UCharVectorType;

  static_assert(Helpers::HasBracketOperator<UCharVectorType>::value,
              "TestHasBracketOperator for CovariantVector<unsigned char, 3> failed!");
  }

  {
  typedef itk::CovariantVector<float, 3> FloatVectorType;

  static_assert(Helpers::HasBracketOperator<FloatVectorType>::value,
              "TestHasBracketOperator for CovariantVector<float, 3> failed!");
  }

  {
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;

  static_assert(Helpers::HasBracketOperator<FloatVectorImageType::PixelType>::value,
              "TestHasBracketOperator for CovariantVector<float, 3> failed!");
  }

  {
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;
  FloatVectorImageType::Pointer image = FloatVectorImageType::New();

  TestHasBracketOperator_Template(image.GetPointer());
  }

  {
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;
  itk::SmartPointer<FloatVectorImageType> image = FloatVectorImageType::New();

  TestHasBracketOperator_Template(image.GetPointer());
  TestHasBracketOperator_ConstTemplate(image.GetPointer());
  }

  {
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> FloatVectorImageType;

  TestHasBracketOperator_Template2<FloatVectorImageType>();
  }

  {
  typedef itk::VariableLengthVector<int> VectorType;

  static_assert(Helpers::HasBracketOperator<VectorType>::value,
              "TestHasBracketOperator for VariableLengthVector failed!");
  }

  {
  typedef std::vector<int> VectorType;

  static_assert(Helpers::HasBracketOperator<VectorType>::value,
              "TestHasBracketOperator for std::vector failed!");
  }

  // This (intentionally) fails
//  {
//  static_assert(Helpers::HasBracketOperator<float>::value,
//              "TestHasBracketOperator for float failed!");
//  }

  return true;
}

template<typename TImage>
bool TestHasBracketOperator_ConstTemplate(const TImage* const )
{
  static_assert(Helpers::HasBracketOperator<typename TImage::PixelType>::value,
              "TestHasBracketOperator for TImage failed!");
  return true;
}

template<typename TImage>
bool TestHasBracketOperator_Template(TImage* )
{
  static_assert(Helpers::HasBracketOperator<typename TImage::PixelType>::value,
              "TestHasBracketOperator for TImage failed!");
  return true;
}


template<typename TImage>
bool TestHasBracketOperator_Template2()
{
  static_assert(Helpers::HasBracketOperator<typename TImage::PixelType>::value,
              "TestHasBracketOperator for TImage failed!");
  return true;
}
