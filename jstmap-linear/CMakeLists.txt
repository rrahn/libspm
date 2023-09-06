
### Base interface target for the index subcommand.
add_library (jstmap_linear_base INTERFACE)
target_include_directories (jstmap_linear_base INTERFACE ../jstmap-linear seqan3::seqan3)
target_link_libraries (jstmap_linear_base INTERFACE seqan3::seqan3 libjst::libjst jstmap::global)
add_library (jstmap::linear::base ALIAS jstmap_linear_base)

### Create object library for better build times

### Create static library for linear subcommand
add_library (jstmap_linear STATIC jstmap/linear/linear_main.cpp)
target_link_libraries (jstmap_linear PUBLIC jstmap_linear_base)
add_library (jstmap::linear ALIAS jstmap_linear)
