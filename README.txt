Obtaining the Source Code
-------------------------
Once you clone this repository, you'll have to also:

git submodule update --init

to get the code in Helpers/ that ITKHelpers relies on.

Building
--------
This repository does not depend on any external libraries. The only caveat is that it depends
on c++0x/11 parts of the c++ language. For Linux, this means it must be built with the flag
gnu++0x. For Windows, we are working on finding the comparable solution/flag.

Usage
-----
This repository follows the structure of our Helpers repository (https://github.com/daviddoria/Helpers), but
contains functions that operate on ITK data structures. That is, ITKHelpers.h contains a namespace called
"ITKHelpers", and ITKStatistics.h contains a namespace called "ITKStatistics". ITKTypeTraits.h contains
overloads to perform statistical operations on ITK data structures and allow the expected return types to be
produced.