# Allow headers in tests to be included like
# #include "ITKHelpers.h" rather than needing
# #include "ITKHelpers/ITKHelpers.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestITKHelpers TestITKHelpers.cpp)
target_link_libraries(TestITKHelpers ${ITKHelpers_libraries})
add_test(TestITKHelpers TestITKHelpers)

add_executable(TestITKContainerInterface TestITKContainerInterface.cpp)
target_link_libraries(TestITKContainerInterface ${ITKHelpers_libraries})
add_test(TestITKContainerInterface TestITKContainerInterface)

add_executable(TestITKStatistics TestITKStatistics.cpp)
target_link_libraries(TestITKStatistics ${ITKHelpers_libraries})
add_test(TestITKStatistics TestITKStatistics)

add_executable(TestITKTypeTraits TestITKTypeTraits.cpp)
target_link_libraries(TestITKTypeTraits ${ITKHelpers_libraries})
add_test(TestITKTypeTraits TestITKTypeTraits)

add_executable(TestColorSpaceConversions TestColorSpaceConversions.cpp)
target_link_libraries(TestColorSpaceConversions ${ITKHelpers_libraries})
add_test(TestColorSpaceConversions TestColorSpaceConversions)

add_executable(TestNormImageAdaptor TestNormImageAdaptor.cpp)
target_link_libraries(TestNormImageAdaptor ${ITKHelpers_libraries})
add_test(TestNormImageAdaptor TestNormImageAdaptor)
