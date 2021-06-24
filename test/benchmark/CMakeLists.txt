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

# Adds a macro benchmark:
function (jstmap_macro_benchmark macro_benchmark_cpp)
    cmake_parse_arguments(
        "MACRO_BENCHMARK"
        ""
        "REPETITIONS"
        "META_DATA"
        ${ARGN}
    )

    if (NOT MACRO_BENCHMARK_REPETITIONS)
        set (MACRO_BENCHMARK_REPETITIONS "3")
    endif ()

    file (RELATIVE_PATH macro_benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${macro_benchmark_cpp}")

    include (seqan3_test_component)
    seqan3_test_component (target "${macro_benchmark}" TARGET_NAME)
    seqan3_test_component (test_name "${macro_benchmark}" TEST_NAME)

    add_executable (${target} ${macro_benchmark_cpp})
    target_link_libraries (${target} jstmap::test::performance jstmap::create jstmap::search)
    add_test (NAME "${test_name}"
              COMMAND ${target}
                      --benchmark_repetitions=${MACRO_BENCHMARK_REPETITIONS}
                      --benchmark_min_time=0
                      --benchmark_out=$<TARGET_FILE_NAME:${target}>.json
                      --benchmark_out_format=json)

    add_dependencies (benchmark_test ${target})

endfunction ()

# Download and build Googlebenchmark module.
seqan3_require_ccache ()
seqan3_require_benchmark ()

add_subdirectories ()
