#ifndef ITKTypeTraits_H
#define ITKTypeTraits_H

#include "Helpers/TypeTraits.h"

#include "itkVariableLengthVector.h"

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

/** For itk::CovariantVector<unsigned char, N>, use itk::CovariantVector<float, N> as the LargerType.
  * This is a partial specialization. */
template <unsigned int N>
struct TypeTraits<itk::CovariantVector<unsigned char, N> >
{
  typedef itk::CovariantVector<float, N> LargerType;
  typedef float LargerComponentType;
  typedef unsigned char ComponentType;
};

#endif
