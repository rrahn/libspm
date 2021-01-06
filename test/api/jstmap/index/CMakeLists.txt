cmake_minimum_required (VERSION 3.8)

macro (add_jstmap_index_test test_filename)
    add_api_test(${test_filename} "jstmap::index")
endmacro ()

add_jstmap_index_test (load_sequence_test.cpp)
