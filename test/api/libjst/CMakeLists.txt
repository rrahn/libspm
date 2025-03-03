cmake_minimum_required (VERSION 3.20)

set (LIBJST_TEST_DEPENDENCIES "jstmap::simulate" "jstmap::search" "libjst::libjst")

macro (add_libjst2_test test_filename)
    add_api_test(${test_filename} "libjst::libjst")
endmacro()

macro (add_libjst_test test_filename)
    add_api_test(${test_filename} "${LIBJST_TEST_DEPENDENCIES}")
endmacro()

add_subdirectories ()

unset (LIBJST_TEST_DEPENDENCIES)
