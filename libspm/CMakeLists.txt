cmake_minimum_required (VERSION 3.20)

### Adds an interface library for contribs to libjst.
add_library (libspm_libspm INTERFACE)
target_include_directories (libspm_libspm INTERFACE ../libspm)
target_compile_features (libspm_libspm INTERFACE cxx_std_20)
target_link_libraries (libspm_libspm INTERFACE seqan3::seqan3 seqan::seqan2)
add_library (libspm::libspm ALIAS libspm_libspm)
