cmake_minimum_required (VERSION 3.14)

# Set directories for test output files, input data and binaries.
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Set the seqan3 specific external project cmake args here as well, since we use the seqan3 external procject to
# fetch goolge test and others.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS}")

# Add more cmake tooling from this project and seqan3.
find_path (JSTMAP_TEST_CMAKE_MODULE_DIR NAMES app_datasources.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "${JSTMAP_TEST_CMAKE_MODULE_DIR}")

find_path (SEQAN3_TEST_CMAKE_MODULE_DIR NAMES seqan3_test_component.cmake
                                        HINTS "${CMAKE_SOURCE_DIR}/lib/seqan3/test/cmake/")
list(APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")

# Build tests just before their execution, because they have not been built with "all" target.
# The trick is here to provide a cmake file as a directory property that executes the build command.
file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target api_test)\n"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target cli_test)\n"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target benchmark_test)")
set_directory_properties (PROPERTIES TEST_INCLUDE_FILE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake")

# Define the test targets. All depending targets are built just before the test execution.
add_custom_target (api_test)
add_custom_target (cli_test)
add_custom_target (benchmark_test)

# Test executables and libraries should not mix with the application files.
unset (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset (CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Define some helper interface libraries for the test
find_path (SEQAN3_TEST_INCLUDE_DIR NAMES seqan3/test/expect_range_eq.hpp
                                   HINTS "${CMAKE_SOURCE_DIR}/lib/seqan3/test/include/")

#--------------------------------------------------------------------
# Cmake interface targets
#--------------------------------------------------------------------

add_library (jstmap_test INTERFACE)
target_include_directories (jstmap_test INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")
target_compile_options (jstmap_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
target_compile_features (jstmap_test INTERFACE cxx_std_20)
target_link_libraries (jstmap_test INTERFACE "pthread" "seqan3::seqan3")
add_library (jstmap::test ALIAS jstmap_test)

add_library (jstmap_test_unit INTERFACE)
target_include_directories (jstmap_test_unit INTERFACE "gtest" "gtest_main" "jstmap::test")
target_link_libraries (jstmap_test_unit INTERFACE "gtest" "gtest_main" "jstmap::test")
add_library (jstmap::test::unit ALIAS jstmap_test_unit)

add_library (jstmap_test_performance INTERFACE)
target_include_directories (jstmap_test_performance INTERFACE "gbenchmark" "jstmap::test")
target_link_libraries (jstmap_test_performance INTERFACE "gbenchmark" "jstmap::test" )
add_library (jstmap::test::performance ALIAS jstmap_test_performance)

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules.
# ----------------------------------------------------------------------------

include (app_datasources)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

include (seqan3_require_benchmark)
include (seqan3_require_ccache)
include (seqan3_require_test)
include (add_subdirectories)
include (seqan3_test_component)

message (STATUS "${FontBold}You can run `make test` to build and run tests.${FontReset}")
