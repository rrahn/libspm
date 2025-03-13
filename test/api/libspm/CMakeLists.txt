cmake_minimum_required (VERSION 3.20)

macro (add_libspm_test test_filename)
    add_api_test(${test_filename})
endmacro()

add_subdirectories ()

add_libspm_test(copyable_box_test.cpp)
