#ifndef ITKTypeTraits_H
#define ITKTypeTraits_H

/** For generic itk::VariableLengthVector, use the same type as the LargerType */
template <typename T>
struct TypeTraits<itk::VariableLengthVector<T> >
{
  typedef itk::VariableLengthVector<T> LargerType;
  typedef typename TypeTraits<T>::LargerType LargerComponentType;
  typedef T ComponentType;
};

/** For itk::VariableLengthVector<unsigned char>, use itk::VariableLengthVector<float> as the LargerType */
template <>
struct TypeTraits<itk::VariableLengthVector<unsigned char> >
{
  typedef itk::VariableLengthVector<float> LargerType;
  typedef float LargerComponentType;
  typedef unsigned char ComponentType;
};

#endif
