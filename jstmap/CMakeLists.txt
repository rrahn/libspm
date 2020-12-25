# Add the main build target.

set (target_name "${CMAKE_PROJECT_NAME}")
add_executable ("${CMAKE_PROJECT_NAME}" jstmap/main.cpp)
target_include_directories ("${CMAKE_PROJECT_NAME}" PUBLIC seqan3::seqan3 ${SEQAN_LIBRARIES} jstmap::index)
target_link_libraries ("${CMAKE_PROJECT_NAME}" PUBLIC seqan3::seqan3 ${SEQAN_LIBRARIES} jstmap::index)
