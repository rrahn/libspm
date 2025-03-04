cmake_minimum_required (VERSION 3.20)

# A macro that adds an api or cli test.
macro (add_app_test test_filename test_alternative target_dependencies)
    enable_testing ()
    # Extract the test target name.
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_test_component (target "${source_file}" TARGET_NAME)
    get_test_component (test_name "${source_file}" TEST_NAME)

    include_directories (${SEQAN_INCLUDE_DIRS})

    # Create the test target.
    set (base_lib_dependency "")
    if (${CMAKE_BUILD_TYPE} STREQUAL "Asan")
        set (base_lib_dependency "jstmap::test::asan")
    elseif (${CMAKE_BUILD_TYPE} STREQUAL "Tsan")
        set (base_lib_dependency "jstmap::test::tsan")
    else ()
        set (base_lib_dependency "jstmap::test::unit")
    endif ()

    add_executable (${target} ${test_filename})
    target_link_libraries (${target} "${base_lib_dependency}" "${target_dependencies}")

    # Add the test to its general target (cli or api).
    if (${test_alternative} STREQUAL "CLI_TEST")
        add_dependencies (${target} "${PROJECT_NAME}") # cli test needs the application executable
        target_include_directories(${target} PUBLIC "${SEQAN3_CLONE_DIR}/test/include")
        add_dependencies (cli_test ${target})
    elseif (${test_alternative} STREQUAL "API_TEST")
        add_dependencies (api_test ${target})
    endif ()

    # Set the target
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

# Special test macro to add a test that uses Catch2.
macro (add_catch2_test test_filename)
    enable_testing ()

    # Extract the test target name.
    file (RELATIVE_PATH source_file "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${test_filename}")
    get_test_component (target "${source_file}" TARGET_NAME)
    get_test_component (test_name "${source_file}" TEST_NAME)

    add_executable (${target} ${test_filename})
    target_link_libraries (${target} PRIVATE "libjst::test::catch2")

    add_dependencies (api_test ${target})

    # Set the target
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (source_file)
    unset (target)
    unset (test_name)
endmacro()

seqan3_require_ccache ()

# Load the tests from the subdirectories
add_subdirectories ()
