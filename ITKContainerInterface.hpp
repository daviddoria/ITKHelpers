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

#ifndef ITKContainerInterface_HPP
#define ITKContainerInterface_HPP

#include "ITKContainerInterface.h"

namespace Helpers
{

template<typename T>
unsigned int length(const itk::RGBPixel<T>& )
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

} // end namespace

#endif
