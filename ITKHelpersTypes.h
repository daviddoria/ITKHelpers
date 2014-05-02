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

#ifndef ITKHelpersTypes_H
#define ITKHelpersTypes_H

// ITK
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkVectorImage.h"

namespace ITKHelpersTypes
{
  /** Scalar types. */
  typedef itk::Image<float, 2> FloatScalarImageType;
  typedef itk::Image<unsigned char, 2> UnsignedCharScalarImageType;

  /** Fixed length vector types. */
  typedef itk::CovariantVector<float, 2> FloatVector2Type;
  typedef itk::Image<FloatVector2Type , 2> FloatVector2ImageType;

  typedef itk::CovariantVector<float, 3> FloatVector3Type;
  typedef itk::Image<FloatVector3Type , 2> FloatVector3ImageType;

  /** RGB types. */
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;

  /** Variable length vector types. */
  typedef itk::VectorImage<float, 2> FloatVectorImageType;
  typedef itk::VectorImage<unsigned char, 2> UnsignedCharVectorImageType;
}

#endif
