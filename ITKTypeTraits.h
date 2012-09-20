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

#ifndef ITKTypeTraits_H
#define ITKTypeTraits_H

// Submodules
#include <Helpers/TypeTraits.h>

// ITK
#include "itkCovariantVector.h"
#include "itkVariableLengthVector.h"
#include "itkRGBPixel.h"

/** For generic itk::VariableLengthVector, use the same type as the LargerType.
  * This is a partial specialization. */
template <typename T>
struct TypeTraits<itk::VariableLengthVector<T> >
{
  typedef itk::VariableLengthVector<T> LargerType;
  typedef typename TypeTraits<T>::LargerType LargerComponentType;
  typedef T ComponentType;
};

/** For itk::VariableLengthVector<unsigned char>, use itk::VariableLengthVector<float> as the LargerType.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::VariableLengthVector<unsigned char> >
{
  typedef itk::VariableLengthVector<float> LargerType;
  typedef float LargerComponentType;
  typedef unsigned char ComponentType;
};

/** For itk::RGBPixel<T>.
  * This is a partial specialization. */
template <typename T>
struct TypeTraits<itk::RGBPixel<T> >
{
  typedef itk::RGBPixel<T> LargerType;
  typedef T LargerComponentType;
  typedef T ComponentType;
};

/** For itk::CovariantVector<T, N>.
  * This is a partial specialization. */
template <typename T, unsigned int N>
struct TypeTraits<itk::CovariantVector<T, N> >
{
  typedef itk::CovariantVector<T, N> LargerType;
  typedef T LargerComponentType;
  typedef T ComponentType;
};

/** For itk::CovariantVector<unsigned char, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<unsigned char, N> >
{
  typedef itk::CovariantVector<float, N> LargerType;
  typedef float LargerComponentType;
  typedef unsigned char ComponentType;
};

/** For itk::CovariantVector<int, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<int, N> >
{
  typedef itk::CovariantVector<float, N> LargerType;
  typedef float LargerComponentType;
  typedef unsigned char ComponentType;
};

#endif
