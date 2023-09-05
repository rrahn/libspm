cmake_minimum_required (VERSION 3.14)

jstmap_benchmark (SOURCE bit_vector_benchmark.cpp DEPENDS libjst::libjst)
jstmap_benchmark (SOURCE journaled_sequence_benchmark.cpp DEPENDS libjst::libjst)
jstmap_benchmark (SOURCE journaled_sequence_search_benchmark.cpp DEPENDS libjst::libjst)
jstmap_benchmark (SOURCE sorted_container_benchmark.cpp DEPENDS libjst::libjst)
