cmake_minimum_required (VERSION 3.20)

### Adds an interface library for contribs to libjst.
add_library (libcontrib_libcontrib INTERFACE)
target_include_directories (libcontrib_libcontrib INTERFACE ../libcontrib)
target_compile_features (libcontrib_libcontrib INTERFACE cxx_std_20)
target_link_libraries (libcontrib_libcontrib INTERFACE seqan3::seqan3 seqan::seqan2)
add_library (libcontrib::libcontrib ALIAS libcontrib_libcontrib)

get_target_property(seqan2_compile_definitions seqan::seqan2 INTERFACE_COMPILE_DEFINITIONS)
message(STATUS "Compile definitions for seqan::seqan2: ${seqan2_compile_definitions}")

get_target_property(libspm_compile_definitions libcontrib::libcontrib INTERFACE_COMPILE_DEFINITIONS)
message(STATUS "Compile definitions for libspm::libspm: ${libspm_compile_definitions}")
