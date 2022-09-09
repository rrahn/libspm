cmake_minimum_required (VERSION 3.14)

set (LIBIO_TEST_DEPENDENCIES "seqan3::seqan3" "libio::libio")

macro (add_libio_test test_filename)
    add_api_test(${test_filename} "${LIBIO_TEST_DEPENDENCIES}")
endmacro()

add_subdirectories ()

unset (LIBIO_TEST_DEPENDENCIES)
