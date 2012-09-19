project(ITKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.8)

# Only include Helpers if it is not already available
get_property(FoundHelpers GLOBAL PROPERTY HelpersDir SET)
if(NOT FoundHelpers)
  add_subdirectory(Helpers)
  message("ITKHelpers submodule in ${CMAKE_CURRENT_SOURCE_DIR} used Helpers submodule (with includes ${Helpers_includes} ).")
endif()

get_property(HelpersDir GLOBAL PROPERTY HelpersDir)
set(ITKHelpers_includes ${ITKHelpers_includes} ${HelpersDir})
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

include_directories(${ITKHelpers_includes})

add_library(ITKHelpers ITKHelpers.cpp)
target_link_libraries(ITKHelpers ${ITKHelpers_libraries})
set(ITKHelpers_libraries ${ITKHelpers_libraries} ITKHelpers)

add_custom_target(ITKHelpersSources SOURCES
ITKHelpers.hpp ITKTypeTraits.h ITKContainerInterface.h ITKContainerInterface.hpp
ITKHelpersTypes.h itkRGBToHSVColorSpacePixelAccessor.h itkRGBToLabColorSpacePixelAccessor.h)

SET(ITKHelpers_BuildTests OFF CACHE BOOL "ITKHelpers_BuildTests")
if(ITKHelpers_BuildTests)
  add_subdirectory(Tests)
endif()

set_property(GLOBAL PROPERTY ITKHelpersDir ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY ITKHelpersLibs ${ITKHelpers_libraries})
