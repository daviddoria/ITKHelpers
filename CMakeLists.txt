cmake_minimum_required(VERSION 2.8)

project(ITKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Only include Helpers if it is not already available
get_property(FoundHelpers GLOBAL PROPERTY HelpersDir SET)
if(NOT FoundHelpers)
  add_subdirectory(Helpers)
endif()

get_property(HelpersIncludeDirs GLOBAL PROPERTY HelpersIncludeDirs)
set(ITKHelpers_include_dirs ${ITKHelpers_include_dirs} ${HelpersIncludeDirs})
get_property(HelpersLibs GLOBAL PROPERTY HelpersLibs)
set(ITKHelpers_libraries ${ITKHelpers_libraries} ${HelpersLibs})

# Required because of the c++0x used in Helpers
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
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
set(ITKHelpers_include_dirs ${ITKHelpers_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY ITKHelpersIncludeDirs ${ITKHelpers_include_dirs})
set_property(GLOBAL PROPERTY ITKHelpersLibs ${ITKHelpers_libraries})

# Build the tests if requested
SET(ITKHelpers_BuildTests OFF CACHE BOOL "ITKHelpers_BuildTests")
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()
