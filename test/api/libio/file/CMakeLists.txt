add_libio_test (mini_test.cpp)
add_libio_test (formatted_stream_test.cpp)
add_libio_test (formatted_file_test.cpp)
target_use_datasources(formatted_file_test FILES in.fasta)

add_libio_test (vcf_file_test.cpp)
target_use_datasources(vcf_file_test FILES test_file.vcf 1KGP.chr22.vcf 1KGP.chr22.test.vcf)

add_libio_test (pivot_tokenizer_test.cpp)
