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

#ifndef ITKContainerInterface_H
#define ITKContainerInterface_H

#include "itkCovariantVector.h"
#include "itkRGBPixel.h"
#include "itkVariableLengthVector.h"

namespace Helpers
{

/** Return the length of the vector through the same interface that we have defined for
 * std::vector and scalars in Helpers.*/
template<typename T>
unsigned int length(const itk::RGBPixel<T>& v);

/** Return the length of the vector through the same interface that we have defined for
 * std::vector and scalars in Helpers.*/
template<typename T>
unsigned int length(const itk::VariableLengthVector<T>& v);

/** Return the length of the vector through the same interface that we have defined for
 * std::vector and scalars in Helpers.*/
template<typename T, unsigned int N>
unsigned int length(const itk::CovariantVector<T, N>& v);

/** Return a reference to the specified component of the vector using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T>
T& index(itk::VariableLengthVector<T>& v, size_t i);

/** Return the specified component of the vector using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T>
T index(const itk::VariableLengthVector<T>& v, size_t i);

/** Return a reference to the specified component of the vector using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T, unsigned int N>
T& index(itk::CovariantVector<T, N>& v, size_t i);

/** Return the specified component of the vector using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T, unsigned int N>
T index(const itk::CovariantVector<T, N>& v, size_t i);

/** Return the specified component of the object using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T>
T index(const itk::RGBPixel<T>& v, size_t i);

/** Return a reference to the specified component of the object using the same interface that we have
 * defined for std::vector and scalars in Helpers.*/
template<typename T>
T& index(itk::RGBPixel<T>& v, size_t i);

} // end namespace

#include "ITKContainerInterface.hpp"

#endif
