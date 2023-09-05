cmake_minimum_required (VERSION 3.14)

# #################################################
# Custom target wrapper
# #################################################

macro (onlin_pattern_ibf_benchmark benchmark_cpp)
    set (JSTMAP_BENCHMARK_TIME_UNIT "s")
    jstmap_benchmark(SOURCE ${benchmark_cpp} DEPENDS jstmap::global jstmap::search)
endmacro ()

# #################################################
# Declare benchmarks
# #################################################

# onlin_pattern_ibf_benchmark(naive_pattern_plain_ibf_benchmark.cpp)
# target_use_datasources(naive_pattern_plain_ibf_benchmark FILES
#                        ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
#                        ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf
#                        needle32.fa
#                        needle64.fa
#                        needle128.fa
#                        needle256.fa)

# Shiftor resumable pattern benchmarks

onlin_pattern_ibf_benchmark(shiftor_pattern_resumable_chr22_needle32_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle32_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_resumable_chr22_needle64_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle64_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_resumable_chr22_needle128_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle128_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_resumable_chr22_needle256_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle256_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

# Shiftor oblivious pattern benchmarks

onlin_pattern_ibf_benchmark(shiftor_pattern_oblivious_chr22_needle32_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle32_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_oblivious_chr22_needle64_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle64_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_oblivious_chr22_needle128_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle128_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(shiftor_pattern_oblivious_chr22_needle256_ibf16384_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle256_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

# Horspool oblivious pattern benchmarks

onlin_pattern_ibf_benchmark(horspool_pattern_oblivious_chr22_needle32_ibf16384_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle32_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(horspool_pattern_oblivious_chr22_needle64_ibf16384_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle64_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(horspool_pattern_oblivious_chr22_needle128_ibf16384_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle128_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(horspool_pattern_oblivious_chr22_needle256_ibf16384_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle256_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

# Myers resumable pattern benchmarks

onlin_pattern_ibf_benchmark(myers_pattern_resumable_chr22_needle32_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle32_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_resumable_chr22_needle64_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle64_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_resumable_chr22_needle128_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle128_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_resumable_chr22_needle256_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle256_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

# Myers oblivious pattern benchmarks

onlin_pattern_ibf_benchmark(myers_pattern_oblivious_chr22_needle32_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle32_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_oblivious_chr22_needle64_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle64_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_oblivious_chr22_needle128_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle128_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

onlin_pattern_ibf_benchmark(myers_pattern_oblivious_chr22_needle256_ibf16384_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle256_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)

# Pigeonhole oblivious pattern benchmarks

onlin_pattern_benchmark(pigeonhole_pattern_oblivious_chr22_sim100x100Ke3_ibf16384_benchmark.cpp)
target_use_datasources(pigeonhole_pattern_oblivious_chr22_sim100x100Ke3_ibf16384_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst sim_reads_chr22_s100_c100K_e3.fa ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf)
