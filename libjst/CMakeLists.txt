cmake_minimum_required (VERSION 3.18)

### Adds an interface library for the jst lib.
add_library (libjst_libjst INTERFACE)
target_include_directories (libjst_libjst INTERFACE ../libjst)
target_compile_features (libjst_libjst INTERFACE cxx_std_20)
target_link_libraries (libjst_libjst INTERFACE libcontrib::libcontrib)
add_library (libjst::libjst ALIAS libjst_libjst)
