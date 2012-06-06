# Allow headers in tests to be included like
# #include "ITKHelpers.h" rather than needing
# #include "ITKHelpers/ITKHelpers.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestITKHelpers TestITKHelpers.cpp)
target_link_libraries(TestITKHelpers ITKHelpers)

add_executable(TestITKContainerInterface TestITKContainerInterface.cpp)
target_link_libraries(TestITKContainerInterface ITKHelpers)
