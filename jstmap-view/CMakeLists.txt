
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/view/*.hpp")

### Base interface target for the index subcommand.
add_library (jstmap_view_base INTERFACE)
target_include_directories (jstmap_view_base INTERFACE ../jstmap-view seqan3::seqan3)
target_link_libraries (jstmap_view_base INTERFACE seqan3::seqan3 libjst::libjst)
add_library (jstmap::view::base ALIAS jstmap_view_base)

### Create object library for better build times
add_library(jstmap_view_load_jst OBJECT jstmap/view/load_jst.cpp
                                        jstmap/view/load_jst.hpp
                                        jstmap/view/global_types.hpp)
target_include_directories (jstmap_view_load_jst PUBLIC jstmap::view::base)
target_link_libraries (jstmap_view_load_jst PUBLIC jstmap::view::base)

add_library(jstmap_view_format_fasta OBJECT jstmap/view/view_format_fasta.cpp
                                            jstmap/view/view_format_fasta.hpp
                                            jstmap/view/global_types.hpp)
target_include_directories (jstmap_view_format_fasta PUBLIC jstmap::view::base)
target_link_libraries (jstmap_view_format_fasta PUBLIC jstmap::view::base)

### Create static library for index subcommand
add_library (jstmap_view STATIC jstmap/view/view_main.cpp)
target_include_directories (jstmap_view PUBLIC jstmap_view_format_fasta jstmap_view_load_jst)
target_link_libraries (jstmap_view PUBLIC jstmap_view_format_fasta jstmap_view_load_jst)
add_library (jstmap::view ALIAS jstmap_view)

unset (HEADER_LIST)
