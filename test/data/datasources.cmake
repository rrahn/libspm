cmake_minimum_required (VERSION 3.8)

include (app_datasources)

# copies file to <build>/data/in.fasta
declare_datasource (FILE in.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/in.fasta
                    URL_HASH SHA256=2c1ccd1b391c45cbbe1b4448584106d2ad2dc996a1636dcfd67342b7f943116a)

declare_datasource (FILE sim_refx5.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_refx5.fasta
                    URL_HASH SHA256=591ba084e2b903246a734699e361724ab97631cad07a4778afb82dad5c0cccd2)

declare_datasource (FILE sim_refx5.jst
                    URL ${CMAKE_SOURCE_DIR}/test/data/sim_refx5.jst
                    URL_HASH SHA256=dde7eeb54065a35105a6eacc379368ef1c85f095005eb2a59ee5ec1bb9a8eb7a)

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
