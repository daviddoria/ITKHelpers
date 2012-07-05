cmake_minimum_required(VERSION 2.8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Only include Helpers if it is not already available
if(NOT TARGET Helpers)
  add_subdirectory(Helpers)
endif()

# This is needed for shared_ptr and the trick using enable_if and if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKImageIntensity
ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology
ITKDeprecated)
INCLUDE(${ITK_USE_FILE})

add_library(ITKHelpers ITKHelpers.cpp)
#target_link_libraries(ITKHelpers LINK_INTERFACE_LIBRARIES Helpers)
target_link_libraries(ITKHelpers Helpers ${ITK_LIBRARIES})

SET(ITKHelpers_BuildTests ON CACHE BOOL "ITKHelpers_BuildTests" FORCE)
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()
