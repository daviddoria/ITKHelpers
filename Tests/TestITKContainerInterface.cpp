#include "ITKContainerInterface.h"

int main()
{
  itk::RGBPixel<unsigned char> rgbPixel;
  std::cout << Helpers::length(rgbPixel) << std::endl;
  std::cout << Helpers::index(rgbPixel, 2) << std::endl;

  itk::VariableLengthVector<float> varVec(4);
  std::cout << Helpers::length(varVec) << std::endl;
  std::cout << Helpers::index(varVec, 1) << std::endl;

  itk::CovariantVector<float, 5> covVec;
  std::cout << Helpers::length(covVec) << std::endl;
  std::cout << Helpers::index(covVec, 1) << std::endl;

  return 0;
}
