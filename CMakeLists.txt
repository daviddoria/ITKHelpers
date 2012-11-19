cmake_minimum_required(VERSION 2.8)

project(ITKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ENABLE_TESTING() # This has to be in the root directory if we want to run 'ctest' from the root build directory.

add_subdirectory(CMakeHelpers)

UseSubmodule(Helpers ITKHelpers)

# Enable c++11
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology)
  INCLUDE(${ITK_USE_FILE})
endif()
set(ITKHelpers_libraries ${ITKHelpers_libraries} ${ITK_LIBRARIES})

# Give the compiler all of the required include directories
include_directories(${ITKHelpers_include_dirs})

# Create the library
add_library(ITKHelpers ITKHelpers.cpp)
target_link_libraries(ITKHelpers ${ITKHelpers_libraries})
set(ITKHelpers_libraries ${ITKHelpers_libraries} ITKHelpers)

# Add all of the non-compiled files to the project
add_custom_target(ITKHelpersSources SOURCES
ITKHelpers.hpp ITKTypeTraits.h ITKContainerInterface.h ITKContainerInterface.hpp
ITKHelpersTypes.h itkRGBToHSVColorSpacePixelAccessor.h itkRGBToLabColorSpacePixelAccessor.h)

# Allow this project to be detected and used as a submodule
CreateSubmodule(ITKHelpers)

# Build the tests if requested
option(ITKHelpers_BuildTests "Build ITKHelpers tests?" OFF)
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()
