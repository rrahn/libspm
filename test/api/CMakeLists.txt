cmake_minimum_required (VERSION 3.8)

# A macro that adds an api or cli test.
macro (add_app_test test_filename test_alternative target_dependencies)
    enable_testing ()
    # Extract the test target name.
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_filename_component (target "${source_file}" NAME_WE)

    include_directories (${SEQAN_INCLUDE_DIRS})
    # Create the test target.
    add_executable (${target} ${test_filename})
    target_link_libraries (${target} "jstmap::test::unit" "${target_dependencies}")

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

seqan3_require_test ()

# Load the tests from the subdirectories
add_subdirectories ()
