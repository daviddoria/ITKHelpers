/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkNormImageAdaptor_h
#define __itkNormImageAdaptor_h

#include "itkImageAdaptor.h"

namespace itk
{
namespace Accessor
{
/** \class NormPixelAccessor
 * \brief Give access to the GetNorm() function of a pixel
 *
 * NormPixelAccessor is templated over an internal type and an
 * external type representation. This class casts the input,
 * applies the function to it, then casts the result according
 * to the types defined as template parameters.
 *
 * \ingroup ImageAdaptors
 * \ingroup ITKImageAdaptors
 */
template< class TInternalType, class TExternalType >
class ITK_EXPORT NormPixelAccessor
{
public:
  /** External typedef. It defines the external aspect
   * that this class will exhibit. */
  typedef TExternalType ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data. */
  typedef TInternalType InternalType;

  static inline TExternalType Get(const TInternalType & input)
  { return (TExternalType)( (double)input.GetNorm() ); }
};
} // end namespace Accessor

/** \class NormImageAdaptor
 * \brief Presents an image as being composed of the GetNorm() of its pixels
 *
 * Additional casting is performed according to the input and output image
 * types following C++ default casting rules.
 *
 * \ingroup ImageAdaptors
 * \ingroup ITKImageAdaptors
 */
template< class TImage, class TOutputPixelType = typename TImage::PixelType::RealValueType >
class ITK_EXPORT NormImageAdaptor:public
  ImageAdaptor< TImage, Accessor::NormPixelAccessor<
                  typename TImage::PixelType,
                  TOutputPixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef NormImageAdaptor Self;
  typedef ImageAdaptor< TImage, Accessor::NormPixelAccessor<
                          typename TImage::PixelType,
                          TOutputPixelType > >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormImageAdaptor, ImageAdaptor);
protected:
  NormImageAdaptor() {}
  virtual ~NormImageAdaptor() {}
private:
  NormImageAdaptor(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
};
} // end namespace itk

#endif
