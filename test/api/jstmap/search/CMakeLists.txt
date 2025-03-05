cmake_minimum_required (VERSION 3.20)

macro (add_jstmap_test test_filename target_dependencies)
    add_api_test(${test_filename} "${target_dependencies}")
endmacro ()

# add_jstmap_test (load_jst_test.cpp "jstmap::search")
# target_use_datasources (load_jst_test FILES sim_refx5_p0.jst)

# add_jstmap_test (load_queries_test.cpp "jstmap::search")
# target_use_datasources (load_queries_test FILES sim_reads_ref1x10.fa)

# add_jstmap_test (search_queries_test.cpp "jstmap::search")

# add_jstmap_test (seed_prefix_seek_position_test.cpp "jstmap::search")
# add_jstmap_test (seed_prefix_node_cargo_test.cpp "jstmap::search")

# add_jstmap_test (pigeonhole_filter_test.cpp "jstmap::search")
# target_use_datasources (pigeonhole_filter_test FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst)

# add_jstmap_test (bucket_searcher_test.cpp "jstmap::search")
# target_use_datasources (bucket_searcher_test FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst)
