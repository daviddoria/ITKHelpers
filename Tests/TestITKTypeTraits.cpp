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

static void TestGenericCovariantVector();
static void TestUnsignedCharCovariantVector();

int main()
{
  TestGenericCovariantVector();
  TestUnsignedCharCovariantVector();


  return 0;
}

void TestGenericCovariantVector()
{
  itk::CovariantVector<int, 3> v;
  v.Fill(0);
  std::cout << v << std::endl;

  TypeTraits<itk::CovariantVector<int, 3> >::ComponentType c = v[0];
  std::cout << c << std::endl;
}

void TestUnsignedCharCovariantVector()
{
  itk::CovariantVector<unsigned char, 3> v;
  v.Fill(0);
  //std::cout << v << std::endl;

  TypeTraits<itk::CovariantVector<unsigned char, 3> >::ComponentType c = 2.3f;
  std::cout << static_cast<float>(c) << std::endl;

  TypeTraits<itk::CovariantVector<unsigned char, 3> >::LargerComponentType larger = 2.3f;
  std::cout << larger << std::endl;
}
