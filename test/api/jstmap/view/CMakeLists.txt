cmake_minimum_required (VERSION 3.14)

macro (add_jstmap_view_test test_filename)
    add_api_test(${test_filename} "jstmap::view")
endmacro ()

# add_jstmap_view_test (view_format_fasta_test.cpp)
# target_use_datasources (view_format_fasta_test FILES sim_refx5.jst)
