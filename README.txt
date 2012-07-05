Obtaining the Source Code
-------------------------
Once you clone this repository, you'll have to also:

git submodule update --init

to get the code in Helpers/ that ITKHelpers relies on.

Building
--------
This repository does not depend on any external libraries. The only caveat is that it depends
on c++0x/11 parts of the c++ language used in the Helpers submodule.
For Linux, this means it must be built with the flag
gnu++0x. For Windows (Visual Studio 2010), nothing special must be done.

Usage
-----
This repository follows the structure of our Helpers repository (https://github.com/daviddoria/Helpers), but
contains functions that operate on ITK data structures. That is, ITKHelpers.h contains a namespace called
"ITKHelpers", and ITKStatistics.h contains a namespace called "ITKStatistics". ITKTypeTraits.h contains
overloads to perform statistical operations on ITK data structures and allow the expected return types to be
produced.

Implementation details
----------------------
A common pattern is to want to call a function on any of itk::VectorImage<T, 2>, itk::Image<T, 2>, as well as
itk::Image<itk::CovariantVector<T, N>, 2>. To accomplish this, we have a vector version of the function,

template <typename T>
FunctionName_Vector(T*)

that provides the vector functionality. We then have a:

template <typename T>
FunctionName(itk::Image<T, 2>*)
{
  // Do scalar image specific things
}

that provides the scalar functionality, and a

template <typename T, unsigned int N>
FunctionName(itk::Image<itk::CovariantVector<T, N>, 2>* image)
{
 FunctionName_Vector(image)
}

that takes the special case of an itk::Image with vector pixels and passes it along to the vector version.

We also provide a

template <typename T>
FunctionName(T* image)
{
 FunctionName_Vector(image)
}

so that the user never has to call the _Vector function directly.
