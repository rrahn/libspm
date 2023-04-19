cmake_minimum_required (VERSION 3.14)

# #################################################
# Declare benchmarks
# #################################################

jstmap_macro_benchmark(naive_pattern_plain_benchmark.cpp)
target_use_datasources(naive_pattern_plain_benchmark FILES
                       ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                       needle32.fa
                       needle64.fa
                       needle128.fa
                       needle256.fa)

jstmap_macro_benchmark(horspool_pattern_oblivious_jst_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_jst_benchmark FILES
                                              ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                                              needle32.fa
                                              needle64.fa
                                              needle128.fa
                                              needle256.fa)

jstmap_macro_benchmark(shiftor_pattern_resumable_jst_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_jst_benchmark FILES
                                              ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                                              needle32.fa
                                              needle64.fa
                                              needle128.fa
                                              needle256.fa)

jstmap_macro_benchmark(shiftor_pattern_oblivious_jst_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_jst_benchmark FILES
                                              ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                                              needle32.fa
                                              needle64.fa
                                              needle128.fa
                                              needle256.fa)
