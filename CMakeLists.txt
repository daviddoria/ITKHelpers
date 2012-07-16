project(ITKHelpers)

cmake_minimum_required(VERSION 2.8)

# Only include Helpers if it is not already available
if(NOT TARGET Helpers)
  add_subdirectory(Helpers)
endif()

set(ITKHelpers_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

# Required because of the c++0x used in Helpers
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase
ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology)
INCLUDE(${ITK_USE_FILE})

add_library(ITKHelpers ITKHelpers.cpp)
target_link_libraries(ITKHelpers Helpers ${ITK_LIBRARIES})

SET(ITKHelpers_BuildTests OFF CACHE BOOL "ITKHelpers_BuildTests")
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()
