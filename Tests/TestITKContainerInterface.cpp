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
