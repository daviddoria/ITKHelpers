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

////////////// VariableLengthVector ///////////////

/** For generic itk::VariableLengthVector, use the same type as the LargerType.
  * This is a partial specialization. */
template <typename T>
struct TypeTraits<itk::VariableLengthVector<T> >
{
  typedef typename TypeTraits<T>::LargerType LargerComponentType;
  typedef itk::VariableLengthVector<LargerComponentType> LargerType;

  typedef T ComponentType;
};

/** For itk::VariableLengthVector<unsigned char>, use itk::VariableLengthVector<float> as the LargerType.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::VariableLengthVector<unsigned char> >
{
  typedef itk::VariableLengthVector<unsigned char> SelfType;

  typedef float LargerComponentType;
  typedef itk::VariableLengthVector<LargerComponentType> LargerType;

  typedef SelfType::ValueType ComponentType;
};

/** For itk::VariableLengthVector<int>, use itk::VariableLengthVector<float> as the LargerType.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::VariableLengthVector<int> >
{
  typedef itk::VariableLengthVector<int> SelfType;

  typedef float LargerComponentType;
  typedef itk::VariableLengthVector<LargerComponentType> LargerType;

  typedef typename SelfType::ValueType ComponentType;
};

/** For itk::VariableLengthVector<unsigned int>, use itk::VariableLengthVector<float> as the LargerType.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::VariableLengthVector<unsigned int> >
{
  typedef itk::VariableLengthVector<unsigned int> SelfType;

  typedef float LargerComponentType;
  typedef itk::VariableLengthVector<LargerComponentType> LargerType;

  typedef typename SelfType::ValueType ComponentType;
};

////////////// CovariantVector ///////////////

/** For itk::CovariantVector<T, N>.
  * This is a partial specialization. */
template <typename T, unsigned int N>
struct TypeTraits<itk::CovariantVector<T, N> >
{
  typedef T LargerComponentType;
  typedef itk::CovariantVector<LargerComponentType, N> LargerType;

  typedef T ComponentType;
};

/** For itk::CovariantVector<unsigned char, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<unsigned char, N> >
{
  typedef itk::CovariantVector<unsigned char, N> SelfType;

  typedef float LargerComponentType;
  typedef itk::CovariantVector<LargerComponentType, N> LargerType;

  typedef typename SelfType::ValueType ComponentType;
};

/** For itk::CovariantVector<int, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<int, N> >
{
  typedef itk::CovariantVector<int, N> SelfType;

  typedef float LargerComponentType;
  typedef itk::CovariantVector<LargerComponentType, N> LargerType;

  typedef typename SelfType::ValueType ComponentType;
};

/** For itk::CovariantVector<unsigned int, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<unsigned int, N> >
{
  typedef itk::CovariantVector<unsigned int, N> SelfType;

  typedef float LargerComponentType;
  typedef itk::CovariantVector<LargerComponentType, N> LargerType;

  typedef typename SelfType::ValueType ComponentType;
};

/////////////////////// RGBPixel ////////////////////
/** For itk::RGBPixel<T>.
  * This is a partial specialization. */
template <typename T>
struct TypeTraits<itk::RGBPixel<T> >
{
  typedef itk::RGBPixel<T> SelfType;

  typedef T LargerComponentType;

  typedef itk::RGBPixel<LargerComponentType> LargerType;

  typedef typename SelfType::ComponentType ComponentType;
};

/** For itk::RGBPixel<unsigned char>.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::RGBPixel<unsigned char> >
{
  typedef itk::RGBPixel<unsigned char> SelfType;

  typedef float LargerComponentType;
  typedef itk::RGBPixel<LargerComponentType> LargerType;

  typedef typename SelfType::ComponentType ComponentType;
};

/** For itk::RGBPixel<unsigned int>.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::RGBPixel<unsigned int> >
{
  typedef itk::RGBPixel<unsigned int> SelfType;

  typedef float LargerComponentType;
  typedef itk::RGBPixel<LargerComponentType> LargerType;

  typedef typename SelfType::ComponentType ComponentType;
};

/** For itk::RGBPixel<unsigned int>.
  * This is a partial specialization. */
template <>
struct TypeTraits<itk::RGBPixel<int> >
{
  typedef itk::RGBPixel<int> SelfType;

  typedef float LargerComponentType;
  typedef itk::RGBPixel<LargerComponentType> LargerType;

  typedef typename SelfType::ComponentType ComponentType;
};

#endif
