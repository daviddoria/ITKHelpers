project(ITKHelpers)

cmake_minimum_required(VERSION 2.8)

# Only include Helpers if it is not already available
if(NOT TARGET Helpers)
  add_subdirectory(Helpers)
  include_directories(${Helpers_includes})
endif()

# Add this directory and all submodule directories to a variable that can
# be used in by a parent project to set the include path
if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(ITKHelpers_includes ${Helpers_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endif()

# Required because of the c++0x used in Helpers
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG
             ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology)
INCLUDE(${ITK_USE_FILE})

add_library(ITKHelpers ITKHelpers.cpp)
target_link_libraries(ITKHelpers Helpers ${ITK_LIBRARIES})

SET(ITKHelpers_BuildTests OFF CACHE BOOL "ITKHelpers_BuildTests")
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()
