cmake_minimum_required (VERSION 3.14)

set (LIBJST_TEST_DEPENDENCIES "libcontrib::libcontrib")

macro (add_libcontrib_test test_filename)
    add_api_test(${test_filename} "${LIBJST_TEST_DEPENDENCIES}")
endmacro()

add_subdirectories ()

add_libcontrib_test(copyable_box_test.cpp)

unset (LIBJST_TEST_DEPENDENCIES)
