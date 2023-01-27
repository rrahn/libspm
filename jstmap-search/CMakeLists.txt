
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/search/*.hpp")

### Base interface target for the index subcommand.
add_library (jstmap_search_base INTERFACE)
target_include_directories (jstmap_search_base INTERFACE ../jstmap-search seqan3::seqan3 seqan2::seqan2)
target_link_libraries (jstmap_search_base INTERFACE OpenMP::OpenMP_CXX seqan3::seqan3 seqan2::seqan2 jstmap::global)
add_library (jstmap::search::base ALIAS jstmap_search_base)

### Create object library for better build times
add_library(jstmap_search_input_queries OBJECT jstmap/search/load_queries.cpp jstmap/search/load_queries.hpp)
target_link_libraries (jstmap_search_input_queries PUBLIC jstmap::search::base)

### Building the journaled sequence tree
add_library(jstmap_search_search OBJECT jstmap/search/search_queries2.cpp jstmap/search/search_queries.hpp)
target_link_libraries (jstmap_search_search PUBLIC jstmap::search::base libjst::libjst)

### Filtering the queries using the additional index
add_library(jstmap_search_filter OBJECT jstmap/search/filter_queries.cpp jstmap/search/filter_queries.hpp)
target_link_libraries (jstmap_search_filter PUBLIC jstmap::search::base libjst::libjst)

### Building the journaled sequence tree
add_library(jstmap_search_output OBJECT jstmap/search/write_results.cpp jstmap/search/write_results.hpp)
target_link_libraries (jstmap_search_output PUBLIC jstmap::search::base libjst::libjst)

### Filtering the queries using the additional index
add_library(jstmap_search_match_aligner OBJECT jstmap/search/match_aligner.cpp jstmap/search/match_aligner.hpp)
target_link_libraries (jstmap_search_match_aligner PUBLIC jstmap::search::base libjst::libjst)

### Create static library for index subcommand
add_library (jstmap_search STATIC jstmap/search/search_main.cpp)
target_link_libraries (jstmap_search PUBLIC jstmap_search_search
                                            jstmap_search_input_queries
                                            jstmap::search::base
                                            jstmap_search_filter
                                            jstmap_search_match_aligner
                                            jstmap::global::bam_writer
                                            # jstmap_search_output
                                            )
add_library (jstmap::search ALIAS jstmap_search)

unset (HEADER_LIST)
