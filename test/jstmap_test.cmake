cmake_minimum_required (VERSION 3.8)

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

# Download and build Googletest module. The interface target 'gtest_all' contains all libs and header paths.
find_path (JSTMAP_TEST_CMAKE_MODULE_DIR NAMES app_datasources.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "${JSTMAP_TEST_CMAKE_MODULE_DIR}")

message (STATUS "Configuring tests. Googletest will be downloaded on demand only.")
include (build_googletest)

# Build tests just before their execution, because they have not been built with "all" target.
# The trick is here to provide a cmake file as a directory property that executes the build command.
file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target api_test)\n"
            "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target cli_test)")
set_directory_properties (PROPERTIES TEST_INCLUDE_FILE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake")

# Define the test targets. All depending targets are built just before the test execution.
add_custom_target (api_test)
add_custom_target (cli_test)

# Test executables and libraries should not mix with the application files.
unset (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset (CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Define some helper interface libraries for the test
find_path (SEQAN3_TEST_INCLUDE_DIR NAMES seqan3/test/expect_range_eq.hpp
                                   HINTS "${CMAKE_SOURCE_DIR}/lib/seqan3/test/include/")

add_library (jstmap_test INTERFACE)
target_include_directories (jstmap_test INTERFACE "gtest_all" "${SEQAN3_TEST_INCLUDE_DIR}")
target_compile_options (jstmap_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
target_link_libraries (jstmap_test INTERFACE "gtest_all" "pthread")
add_library (jstmap::test ALIAS jstmap_test)

# A macro that adds an api or cli test.
macro (add_app_test test_filename test_alternative target_dependencies)
    enable_testing ()
    # Extract the test target name.
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_filename_component (target "${source_file}" NAME_WE)

    include_directories (${SEQAN_INCLUDE_DIRS})
    # Create the test target.
    add_executable (${target} ${test_filename})
    target_link_libraries (${target} "jstmap::test" "${target_dependencies}")

    # Add the test to its general target (cli or api).
    if (${test_alternative} STREQUAL "CLI_TEST")
        add_dependencies (${target} "${PROJECT_NAME}") # cli test needs the application executable
        target_include_directories(${target} PUBLIC "${SEQAN3_CLONE_DIR}/test/include")
        add_dependencies (cli_test ${target})
    elseif (${test_alternative} STREQUAL "API_TEST")
        add_dependencies (api_test ${target})
    endif ()

    # Generate and set the test name.
    get_filename_component (target_relative_path "${source_file}" DIRECTORY)
    if (target_relative_path)
        set (test_name "${target_relative_path}/${target}")
    else ()
        set (test_name "${target}")
    endif ()
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (source_file)
    unset (target)
    unset (test_name)
endmacro ()

# A macro that adds an api test.
macro (add_api_test test_filename target_dependencies)
    add_app_test (${test_filename} API_TEST ${target_dependencies})
endmacro ()

# A macro that adds a cli test.
macro (add_cli_test test_filename target_dependencies)
    add_app_test (${test_filename} CLI_TEST ${target_dependencies})
endmacro ()

# Fetch data and add the tests.
include (app_datasources)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

message (STATUS "${FontBold}You can run `make test` to build and run tests.${FontReset}")
