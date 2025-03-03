cmake_minimum_required (VERSION 3.20)

set (JSTMAP_BENCHMARK_MIN_TIME "1" CACHE STRING "Set --benchmark_min_time= for each benchmark.")
set (JSTMAP_BENCHMARK_TIME_UNIT "ns" CACHE STRING "Set --benchmark_time_unit= for each benchmark.")

macro (jstmap_benchmark)
    set(options "")
    set(oneValueArgs SOURCE REPETITIONS)
    set(multiValueArgs DEPENDS)
    cmake_parse_arguments(MACRO_BENCHMARK "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT MACRO_BENCHMARK_REPETITIONS)
        set (MACRO_BENCHMARK_REPETITIONS "3")
    endif ()

    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${MACRO_BENCHMARK_SOURCE}")

    include (seqan3_test_component)
    seqan3_test_component (target "${benchmark}" TARGET_NAME)
    seqan3_test_component (test_name "${benchmark}" TEST_NAME)

    add_executable (${target} ${MACRO_BENCHMARK_SOURCE})
    target_link_libraries (${target} jstmap::test::performance ${MACRO_BENCHMARK_DEPENDS})
    add_test (NAME "${test_name}"
              COMMAND   ${target}
                        --benchmark_repetitions=${MACRO_BENCHMARK_REPETITIONS}
                        --benchmark_min_time=0
                        --benchmark_time_unit=${JSTMAP_BENCHMARK_TIME_UNIT}
                        --benchmark_out=$<TARGET_FILE_NAME:${target}>.json
                        --benchmark_out_format=json)

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
    jstmap_benchmark(SOURCE ${macro_benchmark_cpp}
                     REPETITIONS ${MACRO_BENCHMARK_REPETITIONS}
                     DEPENDS jstmap::create jstmap::search)
endfunction ()

# Download and build Googlebenchmark module.
seqan3_require_ccache ()
seqan3_require_benchmark ()

add_subdirectories ()
