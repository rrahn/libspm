add_libio_test (mini_test.cpp)
add_libio_test (formatted_stream_test.cpp)
add_libio_test (formatted_file_test.cpp)
target_use_datasources(formatted_file_test FILES in.fasta)
