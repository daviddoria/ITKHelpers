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

#include "ITKTypeTraits.h"

#include "itkCovariantVector.h"

static bool TestVariableLengthVector();
static bool TestCovariantVector();

int main()
{
  bool allPass = true;

  allPass &= TestVariableLengthVector();

  allPass &= TestCovariantVector();

  if(allPass)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}

bool TestCovariantVector()
{
  bool allPass = true;

  // Same larger type
  {
    typedef itk::CovariantVector<double, 3> VectorType;
    typedef itk::CovariantVector<double, 3> CorrectLargerType;
    if(!std::is_same<TypeTraits<VectorType>::LargerType, CorrectLargerType>::value)
    {
      std::cerr << "TestCovariantVector: Same larger type failed!" << std::endl;
      allPass = false;
    }
  }

  // Larger type
  {
    typedef itk::CovariantVector<int, 3> VectorType;
    typedef itk::CovariantVector<float, 3> CorrectLargerType;
    if(!std::is_same<TypeTraits<VectorType>::LargerType, CorrectLargerType>::value)
    {
      std::cerr << "TestCovariantVector: Larger type failed!" << std::endl;
      allPass = false;
    }
  }

  // Component type
  {
    typedef itk::CovariantVector<int, 3> VectorType;
    typedef int CorrectComponentType;
    if(!std::is_same<TypeTraits<VectorType>::ComponentType, CorrectComponentType>::value)
    {
      std::cerr << "TestCovariantVector: Same component type failed!" << std::endl;
      allPass = false;
    }
  }

  // Larger component type
  {
    typedef itk::CovariantVector<int, 3> VectorType;
    typedef float CorrectLargerComponentType;
    if(!std::is_same<TypeTraits<VectorType>::LargerComponentType, CorrectLargerComponentType>::value)
    {
      std::cerr << "TestCovariantVector: Larger component type failed!" << std::endl;
      allPass = false;
    }
  }

  return allPass;
}


bool TestVariableLengthVector()
{
  bool allPass = true;

  // Same larger type
  {
    typedef itk::VariableLengthVector<double> VectorType;
    typedef itk::VariableLengthVector<double> CorrectLargerType;
    if(!std::is_same<TypeTraits<VectorType>::LargerType, CorrectLargerType>::value)
    {
      std::cerr << "TestVariableLengthVector: Same larger type failed!" << std::endl;
      allPass = false;
    }
  }

  // Larger type
  {
    typedef itk::VariableLengthVector<int> VectorType;
    typedef itk::VariableLengthVector<float> CorrectLargerType;
    if(!std::is_same<TypeTraits<VectorType>::LargerType, CorrectLargerType>::value)
    {
      std::cerr << "TestVariableLengthVector: Larger type failed!" << std::endl;
      allPass = false;
    }
  }

  // Component type
  {
    typedef itk::VariableLengthVector<int> VectorType;
    typedef int CorrectComponentType;
    if(!std::is_same<TypeTraits<VectorType>::ComponentType, CorrectComponentType>::value)
    {
      std::cerr << "TestVariableLengthVector: Same component type failed!" << std::endl;
      allPass = false;
    }
  }

  // Larger component type
  {
    typedef itk::VariableLengthVector<int> VectorType;
    typedef float CorrectLargerComponentType;
    if(!std::is_same<TypeTraits<VectorType>::LargerComponentType, CorrectLargerComponentType>::value)
    {
      std::cerr << "TestVariableLengthVector: Larger component type failed!" << std::endl;
      allPass = false;
    }
  }

  return allPass;
}
