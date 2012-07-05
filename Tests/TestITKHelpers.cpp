#include "ITKHelpers.h"

void TestGetClosedContourOrdering();

void TestGetOpenContourOrdering();

void TestDrawRectangle();

void TestIsClosedLoop();

void TestBlurAllChannels();

void TestAnisotropicBlurAllChannels();

void TestHistogramOfGradients();

void TestExtractChannel();
void TestExtractChannels();

void TestSumOfComponentMagnitudes();

void TestGetAllPatchesContainingPixel();

void TestClosestPoint();

void TestDownsample();
void TestUpsample();

void TestDeepCopyFloatScalar();
void TestDeepCopyUnsignedCharScalar();
void TestDeepCopyFloatVector();
void TestDeepCopyUnsignedCharVector();

void TestBreadthFirstOrderingNonZeroPixels();

int main( int argc, char ** argv )
{
  TestGetClosedContourOrdering();

  TestGetOpenContourOrdering();

  TestDrawRectangle();

  TestIsClosedLoop();

  TestBreadthFirstOrderingNonZeroPixels();

  TestBlurAllChannels();

  TestAnisotropicBlurAllChannels();

  TestHistogramOfGradients();

  TestExtractChannel();
  TestExtractChannels();

  TestSumOfComponentMagnitudes();

  TestClosestPoint();

  TestGetAllPatchesContainingPixel();

  TestDownsample();
  TestUpsample();

  TestDeepCopyFloatScalar();
  TestDeepCopyUnsignedCharScalar();
  TestDeepCopyFloatVector();
  TestDeepCopyUnsignedCharVector();

  return 0;
}

void TestBlurAllChannels()
{

  {
  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer blurred = ImageType::New();

  float sigma = 2.0f;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
  }

  {
  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer blurred = ImageType::New();

  float sigma = 2.0f;
  ITKHelpers::BlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
  }

}

void TestAnisotropicBlurAllChannels()
{
  {
  typedef itk::VectorImage<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer blurred = ImageType::New();

  float sigma = 2.0f;
  ITKHelpers::AnisotropicBlurAllChannels(image.GetPointer(), blurred.GetPointer(), sigma);
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
}


void TestGetAllPatchesContainingPixel()
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

}

void TestDeepCopyFloatScalar()
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
}

void TestDeepCopyUnsignedCharScalar()
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
}

void TestDeepCopyFloatVector()
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
}

void TestDeepCopyUnsignedCharVector()
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
}

void TestDownsample()
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
}

void TestUpsample()
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
}

void TestClosestPoint()
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

  unsigned int closestId = ITKHelpers::ClosestPoint(vec, query);

  std::cout << "Closest point to " << query << " is " << vec[closestId] << std::endl;
}

void TestSumOfComponentMagnitudes()
{
  itk::VariableLengthVector<float> a(2);
  a.Fill(3.0);

  float sum = ITKHelpers::SumOfComponentMagnitudes(a);

  std::cout << "Sum of " << a << " is " << sum << std::endl;
}

void TestExtractChannel()
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

void TestExtractChannels()
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
}

void TestHistogramOfGradients()
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
}

void TestBreadthFirstOrderingNonZeroPixels()
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

  for(unsigned int i = 20; i < 30; ++i)
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
}

void TestIsClosedLoop()
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

  for(unsigned int i = 20; i < 30; ++i)
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
}

void TestDrawRectangle()
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
}

void TestGetOpenContourOrdering()
{
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();

  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> region(corner, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  for(unsigned int i = 20; i < 30; ++i)
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
}

void TestGetClosedContourOrdering()
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
}
