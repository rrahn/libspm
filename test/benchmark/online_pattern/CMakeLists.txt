cmake_minimum_required (VERSION 3.14)

# #################################################
# Custom target wrapper
# #################################################

macro (onlin_pattern_benchmark benchmark_cpp)
    set (JSTMAP_BENCHMARK_TIME_UNIT "s")
    jstmap_benchmark(SOURCE ${benchmark_cpp} DEPENDS jstmap::global jstmap::search)
endmacro ()

# #################################################
# Declare benchmarks
# #################################################

onlin_pattern_benchmark(naive_pattern_plain_benchmark.cpp)
target_use_datasources(naive_pattern_plain_benchmark FILES
                       ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst
                       needle32.fa
                       needle64.fa
                       needle128.fa
                       needle256.fa)

# Shiftor resumable pattern benchmarks

onlin_pattern_benchmark(shiftor_pattern_resumable_chr22_needle32_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle32_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa)

onlin_pattern_benchmark(shiftor_pattern_resumable_chr22_needle64_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle64_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa)

onlin_pattern_benchmark(shiftor_pattern_resumable_chr22_needle128_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle128_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa)

onlin_pattern_benchmark(shiftor_pattern_resumable_chr22_needle256_benchmark.cpp)
target_use_datasources(shiftor_pattern_resumable_chr22_needle256_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa)

# Shiftor oblivious pattern benchmarks

onlin_pattern_benchmark(shiftor_pattern_oblivious_chr22_needle32_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle32_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa)

onlin_pattern_benchmark(shiftor_pattern_oblivious_chr22_needle64_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle64_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa)

onlin_pattern_benchmark(shiftor_pattern_oblivious_chr22_needle128_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle128_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa)

onlin_pattern_benchmark(shiftor_pattern_oblivious_chr22_needle256_benchmark.cpp)
target_use_datasources(shiftor_pattern_oblivious_chr22_needle256_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa)

# Horspool oblivious pattern benchmarks

onlin_pattern_benchmark(horspool_pattern_oblivious_chr22_needle32_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle32_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa)

onlin_pattern_benchmark(horspool_pattern_oblivious_chr22_needle64_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle64_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa)

onlin_pattern_benchmark(horspool_pattern_oblivious_chr22_needle128_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle128_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa)

onlin_pattern_benchmark(horspool_pattern_oblivious_chr22_needle256_benchmark.cpp)
target_use_datasources(horspool_pattern_oblivious_chr22_needle256_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa)

# Myers oblivious pattern benchmarks

onlin_pattern_benchmark(myers_pattern_oblivious_chr22_needle32_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle32_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa)

onlin_pattern_benchmark(myers_pattern_oblivious_chr22_needle64_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle64_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa)

onlin_pattern_benchmark(myers_pattern_oblivious_chr22_needle128_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle128_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa)

onlin_pattern_benchmark(myers_pattern_oblivious_chr22_needle256_benchmark.cpp)
target_use_datasources(myers_pattern_oblivious_chr22_needle256_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa)

# Myers resumable pattern benchmarks

onlin_pattern_benchmark(myers_pattern_resumable_chr22_needle32_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle32_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle32.fa)

onlin_pattern_benchmark(myers_pattern_resumable_chr22_needle64_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle64_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle64.fa)

onlin_pattern_benchmark(myers_pattern_resumable_chr22_needle128_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle128_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle128.fa)

onlin_pattern_benchmark(myers_pattern_resumable_chr22_needle256_benchmark.cpp)
target_use_datasources(myers_pattern_resumable_chr22_needle256_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst needle256.fa)

# Pigeonhole oblivious pattern benchmarks

onlin_pattern_benchmark(pigeonhole_pattern_oblivious_chr22_sim100x100Ke3_benchmark.cpp)
target_use_datasources(pigeonhole_pattern_oblivious_chr22_sim100x100Ke3_benchmark FILES ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst sim_reads_chr22_s100_c100K_e3.fa)
