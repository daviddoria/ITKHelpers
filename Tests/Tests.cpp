#include "ITKHelpers.h"

void TestDeepCopyFloatScalar();
void TestDeepCopyUnsignedCharScalar();
void TestDeepCopyFloatVector();
void TestDeepCopyUnsignedCharVector();

int main( int argc, char ** argv )
{
  TestDeepCopyFloatScalar();
  TestDeepCopyUnsignedCharScalar();
  TestDeepCopyFloatVector();
  TestDeepCopyUnsignedCharVector();

  return 0;	
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
