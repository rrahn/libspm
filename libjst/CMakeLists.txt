cmake_minimum_required (VERSION 3.14)

### Adds an interface library for the jst lib.
add_library (libjst_libjst INTERFACE)
target_include_directories (libjst_libjst INTERFACE ../libjst seqan3::seqan3)
target_compile_features (libjst_libjst INTERFACE cxx_std_17)
target_link_libraries (libjst_libjst INTERFACE seqan3::seqan3)
add_library (libjst::libjst ALIAS libjst_libjst)
