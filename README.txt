Obtaining the Source Code
-------------------------
This repository depends on submodules, and sub-submodules. You must clone it with:

git clone --recursive

OR

If you have already cloned the repository with 'git clone', you must also:

git submodule update --init

to get the code in Helpers/ that ITKHelpers relies on.

Dependencies
------------
ITK >= 4

NOTE: If you get errors like:
vxl/core/vnl/vnl_numeric_traits.h:366:29: error: ‘constexpr’ needed for in-class initialization of static data member ‘zero’ of non-integral type

It means that you have not built ITK with c++11 eneabled. To do this, you must add -std=gnu++11 to CMAKE_CXX_FLAGS when you configure ITK:

ccmake ~/src/ITK -DCMAKE_CXX_FLAGS=-std=gnu++11

NOTE: you cannot configure (ccmake) and THEN set CMAKE_CXX_FLAGS - you MUST include the gnu++11 in the ccmake command the very first time it is run.

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
