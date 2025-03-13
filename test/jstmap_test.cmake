cmake_minimum_required (VERSION 3.20)

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


# Add more cmake tooling from this project and seqan3.
find_path (JSTMAP_TEST_CMAKE_MODULE_DIR NAMES app_datasources.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "${JSTMAP_TEST_CMAKE_MODULE_DIR}")

# Build tests just before their execution, because they have not been built with "all" target.
# The trick is here to provide a cmake file as a directory property that executes the build command.
file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target api_test)\n"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target benchmark_test)")
set_directory_properties (PROPERTIES TEST_INCLUDE_FILE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake")

# Define the test targets. All depending targets are built just before the test execution.
add_custom_target (api_test)
add_custom_target (benchmark_test)

# Test executables and libraries should not mix with the application files.
unset (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset (CMAKE_RUNTIME_OUTPUT_DIRECTORY)

#--------------------------------------------------------------------
# Cmake interface targets
#--------------------------------------------------------------------

add_library (jstmap_test INTERFACE)
target_compile_options (jstmap_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
target_compile_features (jstmap_test INTERFACE cxx_std_20)
target_link_libraries (jstmap_test INTERFACE "pthread" "libspm::libspm")
add_library (jstmap::test ALIAS jstmap_test)

add_library (jstmap_test_unit INTERFACE)
target_link_libraries (jstmap_test_unit INTERFACE "GTest::gtest" "GTest::gtest_main" "jstmap::test")
add_library (jstmap::test::unit ALIAS jstmap_test_unit)

# Add dedicated unit test interface for libjst using Catch2 as test framework.
add_library (libjst_test_catch2 INTERFACE)
target_link_libraries (libjst_test_catch2 INTERFACE "jstmap::test" "Catch2::Catch2WithMain")
add_library (libjst::test::catch2 ALIAS libjst_test_catch2)

add_library (jstmap_test_performance INTERFACE)
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
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

include (add_subdirectories)
include (get_test_component)
enable_testing ()

message (STATUS "${FontBold}You can run `make test` to build and run tests.${FontReset}")
