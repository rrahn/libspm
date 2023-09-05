cmake_minimum_required (VERSION 3.8)

include (app_datasources)

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

declare_datasource (FILE Ash1_v2.2.fa.gz USE_GUNZIP_EXTRACT
                    URL ftp://ftp.ccb.jhu.edu/pub/data/Homo_sapiens/Ash1/v2.2/Assembly/Ash1_v2.2.fa.gz)

declare_datasource (FILE ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/jstmap/current/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst
                    URL_HASH SHA256=9f787f1f0b14375a9da695d42457da13a7cf8a4ee977ff5b0d1867717b2acd17)

declare_datasource (FILE ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst.ibf
                    URL https://ftp.imp.fu-berlin.de//pub/rmaerker/jstmap/current/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst.ibf
                    URL_HASH SHA256=c1cc172b447bb0b4d95b62a1b948e5c23e8bb82b743c38aba0fac071b6399a33)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.gz
                    URL_HASH SHA256=907f7aa263e3357c6c6c70929ffac06268918ff4b29814f505bbfc2d01a2339b)

declare_datasource (FILE needle32.fa USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/needle32.fa.gz
                    URL_HASH SHA256=182ed3ee816c1b8a6fa88f7e6b09ae2d35b702b0a372628daeeae3fa6cadffde)

declare_datasource (FILE needle64.fa USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/needle64.fa.gz
                    URL_HASH SHA256=66705d6e5159e1d9675e83e430f1ecf853e28f811d38f77d3576262f9e7c9a2f)

declare_datasource (FILE needle128.fa USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/needle128.fa.gz
                    URL_HASH SHA256=0043b0dfa4a5418afdc9cc1ae3a4113b176cf40af47b5d871de343a0138fcaf0)

declare_datasource (FILE needle256.fa USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/needle256.fa.gz
                    URL_HASH SHA256=bf19935e362516d9046b252726baf769ab33b6a33e6bee0403905bf95398a7b8)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.b6k.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.b6k.k21.ibf.gz
                    URL_HASH SHA256=3e9e5eadd8ba599f4382de09a02530eac70b84505853e1358c38240c73461be9)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c64.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c64.k21.ibf.gz
                    URL_HASH SHA256=cd2b5d0b451a39644105cf4694a0bf176c5ab18a295bc81b7e452d1d8024252f)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c128.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c128.k21.ibf.gz
                    URL_HASH SHA256=64ba4b13f568fd9ea86c3d7fbb0671b994d50e947caaea1760c5e7f36cadaacc)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c256.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c256.k21.ibf.gz
                    URL_HASH SHA256=869b574fab0abdda63aeb96e3a60e6c296295f94d1c7c545b2a58af0d3229c6f)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c512.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c512.k21.ibf.gz
                    URL_HASH SHA256=bad977e0526cd4e37068cdea88ce42252e157ec1e6eece2e0d48f588c0799697)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c1024.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c1024.k21.ibf.gz
                    URL_HASH SHA256=959ee8a0f68428c921b18485cc13d161b37635f9a8221c10d27394e15af6f5c1)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c2048.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c2048.k21.ibf.gz
                    URL_HASH SHA256=874ae92db560e8259e93ff753bc7bdc8d07df2f463a601ba5c790d4db0ceaeea)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c4096.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c4096.k21.ibf.gz
                    URL_HASH SHA256=d76912394c16c6216e827bf60cf0344a1ddd0412e54bccf5094bf11634f01e34)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c8192.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c8192.k21.ibf.gz
                    URL_HASH SHA256=0aafb8a0b481c55a1338fd73fca69d3f6523f1bfb326895feac15d1264eb1bc1)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c16384.k21.ibf.gz
                    URL_HASH SHA256=29e343432f534dec113fe5feebd77d7f40e632bbe7dc38afe3f2eba3bff44e27)

declare_datasource (FILE ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf.gz
                    URL_HASH SHA256=1edac6e53ad441398ca3ff9016d26c9fbaa526ea9a302e79ee0e797f7141725e)

declare_datasource (FILE sim_reads_chr22_s100_c100K_e3.fa USE_GUNZIP_EXTRACT
                    URL https://ftp.imp.fu-berlin.de/pub/rmaerker/just_bench/v0.0.1/sim_reads_chr22_s100_c100K_e3.fa.gz
                    URL_HASH SHA256=8f4ee4d8ac66d2a72fd4da33591e39e2b6ad7e1fbbf258ad0708f5af255e8a22)
