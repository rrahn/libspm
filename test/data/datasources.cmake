cmake_minimum_required (VERSION 3.8)

include (app_datasources)

# copies file to <build>/data/in.fasta
declare_datasource (FILE in.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/in.fasta
                    URL_HASH SHA256=2c1ccd1b391c45cbbe1b4448584106d2ad2dc996a1636dcfd67342b7f943116a)
