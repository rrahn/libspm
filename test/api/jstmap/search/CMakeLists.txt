cmake_minimum_required (VERSION 3.14)

macro (add_jstmap_test test_filename target_dependencies)
    add_api_test(${test_filename} "${target_dependencies}")
endmacro ()

add_jstmap_test (load_jst_test.cpp "jstmap::search")
target_use_datasources (load_jst_test FILES sim_refx5_p0.jst)

add_jstmap_test (load_queries_test.cpp "jstmap::search")
target_use_datasources (load_queries_test FILES sim_reads_ref1x10.fa)

add_jstmap_test (search_queries_test.cpp "jstmap::search")
