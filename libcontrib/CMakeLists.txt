cmake_minimum_required (VERSION 3.23)

### Adds an interface library for contribs to libjst.
add_library (libcontrib_libcontrib INTERFACE)
target_sources(libcontrib_libcontrib INTERFACE FILE_SET HEADERS BASE_DIRS ../libcontrib)
target_compile_features (libcontrib_libcontrib INTERFACE cxx_std_20)
target_link_libraries (libcontrib_libcontrib INTERFACE seqan3::seqan3 seqan2::seqan2)
add_library (libcontrib::libcontrib ALIAS libcontrib_libcontrib)
