include(ProcessorCount)

# #################################################
# Declare simulation targets
# #################################################

ProcessorCount(num_cores)
math(EXPR threads "${num_cores} / 2" OUTPUT_FORMAT DECIMAL)

# simulate(SIMULATION_NAME simulate_vg
#          REFERENCE_FILE Ash1_v2.2.fa.gz
#          READ_COUNT 40000000
#          READ_LENGTH 150
#          SNP_PROB 0.01
#          INDEL_PROB 0.002
#          THREAD_COUNT "${threads}")

# # #################################################
# # Declare benchmarks
# # #################################################

# jstmap_macro_benchmark(
#     search_benchmark.cpp
#     REPETITIONS 10 # optional, default: 3
# )

# target_use_datasources (search_benchmark FILES ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst)
# # target_use_datasources (search_benchmark FILES simulate_vg)

# jstmap_macro_benchmark(pigeonhole_filter_benchmark.cpp)
# target_use_datasources(pigeonhole_filter_benchmark FILES sim_ref_10Kb_SNP_INDELs.jst)

# jstmap_macro_benchmark(
#     context_enumeration_benchmark.cpp
#     REPETITIONS 10 # optional, default: 3
# )
# target_use_datasources (context_enumeration_benchmark FILES sim_ref_10Kb_SNP_INDELs.jst
#                                                             sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz)
