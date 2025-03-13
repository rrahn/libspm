cmake_minimum_required (VERSION 3.20)

CPMGetPackage (seqan3)

if (NOT DATA_ROOT_DIR)
    set(DATA_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# Set directories for test output files, input data and binaries.
set (DATA_DIR "${DATA_ROOT_DIR}/data/")
message(STATUS "Using DATA_ROOT_DIR: ${DATA_ROOT_DIR}")

file (MAKE_DIRECTORY ${DATA_DIR})
file (MAKE_DIRECTORY ${DATA_ROOT_DIR}/output)

add_definitions (-DOUTPUTDIR=\"${DATA_ROOT_DIR}/output/\")
add_definitions (-DDATADIR=\"${DATA_DIR}\")
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
                                        HINTS "${seqan3_SOURCE_DIR}/test/cmake/")
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
                                   HINTS "${seqan3_SOURCE_DIR}/test/include/")

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
target_include_directories (jstmap_test_unit INTERFACE "GTest::gtest" "GTest::gtest_main" "jstmap::test")
target_link_libraries (jstmap_test_unit INTERFACE "GTest::gtest" "GTest::gtest_main" "jstmap::test")
add_library (jstmap::test::unit ALIAS jstmap_test_unit)

# Add dedicated unit test interface for libjst using Catch2 as test framework.
add_library (libjst_test_catch2 INTERFACE)
target_link_libraries (libjst_test_catch2 INTERFACE "jstmap::test" "libjst::libjst" "Catch2::Catch2WithMain")
add_library (libjst::test::catch2 ALIAS libjst_test_catch2)

add_library (jstmap_test_performance INTERFACE)
target_include_directories (jstmap_test_performance INTERFACE "benchmark::benchmark" "jstmap::test")
target_link_libraries (jstmap_test_performance INTERFACE "benchmark::benchmark" "jstmap::test" )
add_library (jstmap::test::performance ALIAS jstmap_test_performance)

add_library (jstmap_test_asan INTERFACE)
target_compile_options(jstmap_test_asan INTERFACE "${JST_SANITIZER_FLAGS}")
target_link_options(jstmap_test_asan INTERFACE "${JST_SANITIZER_FLAGS}")
target_link_libraries (jstmap_test_asan INTERFACE "jstmap::test::unit")
add_library (jstmap::test::asan ALIAS jstmap_test_asan)

add_library (jstmap_test_tsan INTERFACE)
target_compile_options(jstmap_test_tsan INTERFACE "${JST_SANITIZER_FLAGS}")
target_link_options(jstmap_test_tsan INTERFACE "${JST_SANITIZER_FLAGS}")
target_link_libraries (jstmap_test_tsan INTERFACE "jstmap::test::unit")
add_library (jstmap::test::tsan ALIAS jstmap_test_tsan)

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules.
# ----------------------------------------------------------------------------

include (app_datasources)
include (read_simulation)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

include (add_subdirectories)
include (get_test_component)
enable_testing ()

message (STATUS "${FontBold}You can run `make test` to build and run tests.${FontReset}")
