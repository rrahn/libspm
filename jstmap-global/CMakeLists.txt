### Base interface target for the global jstmap library.
add_library (jstmap_global_base INTERFACE)
target_include_directories (jstmap_global_base INTERFACE ../jstmap-global
                                                         libjst::libjst
                                                         seqan3::seqan3
                                                         ${SEQAN_INCLUDE_DIRS})
target_compile_features (jstmap_global_base INTERFACE cxx_std_20)
target_link_libraries (jstmap_global_base INTERFACE libjst::libjst seqan3::seqan3 ${SEQAN_LIBRARIES})
add_library (jstmap::global::base ALIAS jstmap_global_base)

### Create object library for the logging mechanism.
add_library(jstmap_global_logging OBJECT jstmap/global/application_logger.cpp jstmap/global/application_logger.hpp)
target_link_libraries (jstmap_global_logging PUBLIC jstmap::global::base)

### Create object library for loading a jst.
add_library(jstmap_global_load_jst STATIC jstmap/global/load_jst.cpp jstmap/global/load_jst.hpp)
target_link_libraries (jstmap_global_load_jst PUBLIC jstmap::global::base jstmap_global_logging)

### Expose global library name.
add_library (jstmap::global ALIAS jstmap_global_load_jst)
