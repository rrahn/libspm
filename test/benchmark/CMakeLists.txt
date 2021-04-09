cmake_minimum_required (VERSION 3.14)

set (JSTMAP_BENCHMARK_MIN_TIME "1" CACHE STRING "Set --benchmark_min_time= for each benchmark. Timings are unreliable in CI.")

macro (jstmap_benchmark benchmark_cpp link_libraries)
    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${benchmark_cpp}")
    seqan3_test_component (target "${benchmark}" TARGET_NAME)
    seqan3_test_component (test_name "${benchmark}" TEST_NAME)

    add_executable (${target} ${benchmark_cpp})
    target_link_libraries (${target} jstmap::test::performance ${link_libraries})
    add_test (NAME "${test_name}" COMMAND ${target} "--benchmark_min_time=${JSTMAP_BENCHMARK_MIN_TIME}")

    add_dependencies (benchmark_test ${target})

    unset (benchmark)
    unset (target)
    unset (test_name)
endmacro ()

# Download and build Googlebenchmark module.
seqan3_require_benchmark ()

add_subdirectories ()
