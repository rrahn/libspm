cmake_minimum_required (VERSION 3.14)

### Adds an interface library for the vcf lib.
add_library (libio_lib INTERFACE)
target_include_directories (libio_lib INTERFACE ../libio seqan3::seqan3)
target_compile_features (libio_lib INTERFACE cxx_std_17)
target_link_libraries (libio_lib INTERFACE seqan3::seqan3)
add_library (libio::libio ALIAS libio_lib)
