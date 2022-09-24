cmake_minimum_required (VERSION 3.8)

include (app_datasources)

## Test sources

declare_datasource (FILE test_file.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_file.vcf
                    URL_HASH SHA256=ebce59acc2af48e8cb4258e43f8837073b3a0cac3b85cb0350a965b31d6be23d)

declare_datasource (FILE 1KGP.chr22.test.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/1KGP.chr22.test.vcf
                    URL_HASH SHA256=ed13488e95d800aa3a1ddd915c9d5872ae0af3f9abb71c2f557114c9703455b5)

declare_datasource (FILE 1KGP.chr22_20k.vcf
                    URL ${CMAKE_SOURCE_DIR}/../data/jstmap/new/1KGP.chr22_20k.vcf
                    URL_HASH SHA256=4152ee4f243159196885fc571256c698df48d6a56a1d8a023e04de6ff0b92a33)

declare_datasource (FILE 1KGP.chr22_100k.vcf
                    URL ${CMAKE_SOURCE_DIR}/../data/jstmap/new/1KGP.chr22_100k.vcf
                    URL_HASH SHA256=a6638912a752b160845e78b02c80df64cc9c2375d5e50d429ccbece5c85f0747)

declare_datasource (FILE 1KGP.chr22.vcf.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/1KGP.chr22.vcf.jst
                    URL_HASH SHA256=e7fb27a9a3d4ed18ff744f542613b33c3f5b9eae5726ff26fb8d4646db2d28a9)

# copies file to <build>/data/in.fasta
declare_datasource (FILE empty.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/empty.fa
                    URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)

declare_datasource (FILE in.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/in.fasta
                    URL_HASH SHA256=2c1ccd1b391c45cbbe1b4448584106d2ad2dc996a1636dcfd67342b7f943116a)

declare_datasource (FILE sim_refx5.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_refx5.fasta
                    URL_HASH SHA256=591ba084e2b903246a734699e361724ab97631cad07a4778afb82dad5c0cccd2)


# The jst was build from the following alignments, where the target sequences come from sim_refx5.fasta
# reference:
# TTCACGATGGAACCTAGGCAGAAATACGATCAGGTTGCGACCCTTGGCCGACTATC"
#
# alignment1:
# T-T-CACGA---T--GGAA-CCT-A----GGCA-G-A-AA-TA---C-GAT-CA--GG-T---TG--C-GACCC--TTGG-CC-G--A-CTATC
# TATGCACCAGAGTATGGAAGCATAAGCTCTGCATGCAAAGGTACATCAGATCCTGCGGTTGGGTGCCAACCCAAGTGTGTTCACGGGCGC----
#
# alignment2:
# T-T-C------ACGA---T--GGAA-CCT-A----GGCA-G-A-AA-TA---C-GATCA--GGTTG-C-GAC-CC--TTGG-CC-G--A-C-T-----A-T-C-
# TTGACAGACATCGGAGGATGGTGCACACTCACTCGACCAGCGCAAAGCACAGGATCTCACGG----GCGGACATCTCTTAGGTCAGTCATCGTGGAGGAATGCT
#
# alignment3:
# T-TCACGA---T--GGAA--CCT-A----GGCA-G-A-AA-TA---C-GAT-CA--GG-T---TG-C-GAC-CC--TTGG-CC-G--A-C-T-----A-T-C---------
# TGT-ACGTTCTTTTGGCTTCCCCTAACACGGCGGGCGTCTCCGGTACGTATCCTGTCGGTACACCCCTTAAGCCCCTAGGCCCGAAGAACATAGCGCATTTCACGCTCTCT
#
# alignment4:
# TTCACGA---T--GGAA-CCT-A----GGCA-G-A-AA-TA---C-GAT-CA--GG-T---TG--C-GAC-CC--TTGG-CC-G--A-CTATC
# ---ACGAATGACCGCAACGATCAAATGGGCGAGAACAACTAATTCCGATTCATGGGGTTTGTGGATTGTGACACAGCGCGCCCGCTAC-----
#
# alignment5:
# TT-CACGA--T--GGAA-CCT-A----GGCA-G-A-AA-TA---C-GAT-CA--GG-T---TG--C-GAC-CC--TTGG-CC-G--A-C-T-----A-T-C----
# -TGCGGGACGTGAGGACGCCCAATTCTGCCAAGGATTATTTAGGGTGTTTCACTAGAGTTATGCGCCGACCCCGGTTGGACCAGCTTGCATTCGAAACTGCGTTA

declare_datasource (FILE sim_refx5.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_refx5.jst
                    URL_HASH SHA256=ce5b3c884df5eaa39326b1270250946185422071b2f6d79f137ea184a2f66c32)

declare_datasource (FILE sim_refx5_p0.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_refx5_p0.jst
                    URL_HASH SHA256=0cc0b0c72f06ca54cb4c7afb2252d488c64804e8425ff7d39a58026eb64e77df)

declare_datasource (FILE sim_reads_ref1x10.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_reads_ref1x10.fa
                    URL_HASH SHA256=beabdfaf45099218936976f77aad763799cf6d8188f59c5763c82dc71b6f6f49)

declare_datasource (FILE sim_reads_ref2x10.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_reads_ref2x10.fa
                    URL_HASH SHA256=4f3293c25173344b2b41604cd4a06a25fa36e2fd0fdc55ff03060498931db27c)

declare_datasource (FILE sim_reads_ref3x10.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_reads_ref3x10.fa
                    URL_HASH SHA256=150fcb7d6e9401be5dc5ac26b1d46205b3f540a40cd2a416d74a721d66929e73)

declare_datasource (FILE sim_reads_ref4x10.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_reads_ref4x10.fa
                    URL_HASH SHA256=80fda2dbb2ee83f0d56baccf18200d98802982365a02facc01e79f5fe84cac2b)

declare_datasource (FILE sim_reads_ref5x10.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_reads_ref5x10.fa
                    URL_HASH SHA256=91d775282c7a04d7fd11b66e22c1b42fe947ac7948a35596e538982c95ea4820)

# Simulated vcf files

declare_datasource (FILE sim_ref_10Kb.fasta.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb.fasta.gz
                    URL_HASH SHA256=d4c76490c23668387228f5b9ada7d86212b38a1656f068f820b7d161ce8dce0a)

declare_datasource (FILE sim_ref_10Kb_SNPs.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNPs.vcf
                    URL_HASH SHA256=652694c7fc52fe77baeb6fcab8c4a36a41c8bdb40b98e48d27d2b73d913f2067)

declare_datasource (FILE sim_ref_10Kb_SNPs_haplotypes.fasta.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNPs_haplotypes.fasta.gz
                    URL_HASH SHA256=d9d37d7474884a7bc4bae99c5cdd6457cee6438bbf63ef48c0959dbff8fb1bca)

declare_datasource (FILE sim_ref_10Kb_SNPs.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNPs.jst
                    URL_HASH SHA256=8fe5dba7cbe8b906f3482b24d371e9f9919db0608aa95d5cac7052138331ac84)

declare_datasource (FILE sim_ref_10Kb_SNP_INDELs.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNP_INDELs.vcf
                    URL_HASH SHA256=3f9f65b60d160094f400a3695e6b3685d0275a0e46b5e11d3eed7846d28f088f)

declare_datasource (FILE sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz
                    URL_HASH SHA256=993400d2fd437cc49ddffb296d14cba505e85e5f5e13d3ab5609aa0c46d26d99)

declare_datasource (FILE sim_ref_10Kb_SNP_INDELs.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_SNP_INDELs.jst
                    URL_HASH SHA256=6144893bf71f0eb3fa9339ce7ea9339c392f1ddce705ad952e09d86e8c192a77)

declare_datasource (FILE sim_ref_10Kb_no_variants.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_ref_10Kb_no_variants.vcf
                    URL_HASH SHA256=fa0b191b7cba9e3da323bf7674767a4efce0830ee66750ad6c3b8d8952bbaef1)

# Data sources for benchmarking

declare_datasource (FILE Ash1_v2.2.fa.gz
                    URL ftp://ftp.ccb.jhu.edu/pub/data/Homo_sapiens/Ash1/v2.2/Assembly/Ash1_v2.2.fa.gz
                    URL_HASH "")

declare_datasource (FILE 1KGP.chr22.vcf.gz
                    URL /project/archive-index-data/rahn/1KGP/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
                    URL_HASH "")
