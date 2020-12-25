
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/index/*.hpp")

add_library (jstmap_index STATIC jstmap/index/index_main.cpp ${HEADER_LIST})
target_include_directories (jstmap_index PUBLIC ../jstmap-index)
target_compile_features (jstmap_index PRIVATE cxx_std_17)
target_link_libraries (jstmap_index PUBLIC seqan3::seqan3 ${SEQAN_LIBRARIES})
add_library (jstmap::index ALIAS jstmap_index)

unset (HEADER_LIST)
