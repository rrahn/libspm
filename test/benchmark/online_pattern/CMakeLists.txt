cmake_minimum_required (VERSION 3.14)

# #################################################
# Declare benchmarks
# #################################################

jstmap_macro_benchmark(naive_pattern_linear_benchmark.cpp)
target_use_datasources(naive_pattern_linear_benchmark FILES
                       ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                       needle32.fa
                       needle64.fa
                       needle128.fa
                       needle256.fa)
