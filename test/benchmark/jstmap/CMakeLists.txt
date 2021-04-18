
jstmap_macro_benchmark(
    search_benchmark.cpp
    REPETITIONS 10 # optional, default: 3
)
target_use_datasources (search_benchmark FILES sim_ref_10Kb.fasta.gz
                                               sim_ref_10Kb_SNP_INDELs.vcf
                                               sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz
                                               sim_reads_ref1x10.fa)


jstmap_macro_benchmark(
    context_enumeration_benchmark.cpp
    REPETITIONS 10 # optional, default: 3
)
target_use_datasources (context_enumeration_benchmark FILES sim_ref_10Kb.fasta.gz
                                                            sim_ref_10Kb_SNP_INDELs.vcf
                                                            sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz)
