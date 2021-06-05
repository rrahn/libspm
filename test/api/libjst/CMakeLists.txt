cmake_minimum_required (VERSION 3.14)

set (LIBJST_TEST_DEPENDENCIES "jstmap::simulate" "jstmap::search" "libjst::libjst")

macro (add_libjst_test test_filename)
    add_api_test(${test_filename} "${LIBJST_TEST_DEPENDENCIES}")
endmacro()

add_libjst2_test (journal_test.cpp)
add_libjst2_test (journal_sequence_iterator_test.cpp)

add_subdirectories ()

unset (LIBJST_TEST_DEPENDENCIES)
