cmake_minimum_required (VERSION 3.8)

macro (add_jstmap_search_test test_filename)
    add_api_test(${test_filename} "jstmap::search")
endmacro ()

add_jstmap_search_test (load_jst_test.cpp)
target_use_datasources (load_jst_test FILES sim_refx5.jst)

add_jstmap_search_test (load_queries_test.cpp)
target_use_datasources (load_queries_test FILES sim_reads_ref1x10.fa)

add_jstmap_search_test (search_queries_test.cpp)
target_use_datasources (search_queries_test FILES sim_refx5.jst sim_reads_ref1x10.fa)
