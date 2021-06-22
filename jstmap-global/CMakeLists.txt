### Base interface target for the global jstmap library.
add_library (jstmap_global_base INTERFACE)
target_include_directories (jstmap_global_base INTERFACE ../jstmap-global
                                                         libjst::libjst
                                                         seqan3::seqan3
                                                         ${SEQAN_INCLUDE_DIRS})
target_compile_features (jstmap_global_base INTERFACE cxx_std_17)
target_link_libraries (jstmap_global_base INTERFACE libjst::libjst seqan3::seqan3 ${SEQAN_LIBRARIES})
add_library (jstmap::global ALIAS jstmap_global_base)
