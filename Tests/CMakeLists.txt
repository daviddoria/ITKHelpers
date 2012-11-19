# Allow headers in tests to be included like
# #include "ITKHelpers.h" rather than needing
# #include "ITKHelpers/ITKHelpers.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

GET_DIRECTORY_PROPERTY(ITKHelpersTestsIncludes INCLUDE_DIRECTORIES)
message("ITKHelpersTests include directories: ${ITKHelpersTestsIncludes}")

message("ITKHelpersTests ITKHelpers_includes: ${ITKHelpers_includes}")

add_executable(TestITKHelpers TestITKHelpers.cpp)
target_link_libraries(TestITKHelpers ${ITKHelpers_libraries})

add_executable(TestITKContainerInterface TestITKContainerInterface.cpp)
target_link_libraries(TestITKContainerInterface ${ITKHelpers_libraries})

add_executable(TestITKStatistics TestITKStatistics.cpp)
target_link_libraries(TestITKStatistics ${ITKHelpers_libraries})

add_executable(TestITKTypeTraits TestITKTypeTraits.cpp)
target_link_libraries(TestITKTypeTraits ${ITKHelpers_libraries})

add_executable(TestColorSpaceConversions TestColorSpaceConversions.cpp)
target_link_libraries(TestColorSpaceConversions ${ITKHelpers_libraries})

add_executable(TestNormImageAdaptor TestNormImageAdaptor.cpp)
target_link_libraries(TestNormImageAdaptor ${ITKHelpers_libraries})
