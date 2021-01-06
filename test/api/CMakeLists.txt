cmake_minimum_required (VERSION 3.8)

include (../jstmap_test.cmake)

# TODO: infra improvements
# include (diagnostics/list_unused_unit_tests) TODO: Enable me!

# jstmap_require_ccache ()
# jstmap_require_test ()

#list_unused_unit_tests ()

### add subdirectories:

add_subdirectory (jstmap)
add_subdirectory (libjst)

# add_api_test (convert_fastq_test.cpp)
# target_use_datasources (convert_fastq_test FILES in.fastq)
