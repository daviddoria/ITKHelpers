#include "ITKHelpers.h"

void TestGetAllPatchesContainingPixel();

void TestDownsample();

void TestDeepCopyFloatScalar();
void TestDeepCopyUnsignedCharScalar();
void TestDeepCopyFloatVector();
void TestDeepCopyUnsignedCharVector();

int main( int argc, char ** argv )
{
  TestGetAllPatchesContainingPixel();

  TestDownsample();
  
  TestDeepCopyFloatScalar();
  TestDeepCopyUnsignedCharScalar();
  TestDeepCopyFloatVector();
  TestDeepCopyUnsignedCharVector();

  return 0;
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
