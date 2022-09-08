include(ExternalProject)

function(floatexpr expr output)
    execute_process(COMMAND awk "BEGIN {print ${expr}}" OUTPUT_VARIABLE __output)
    string(REGEX REPLACE "[ \t\r\n]" "" __trimmed_output ${__output})
    set(${output} ${__trimmed_output} PARENT_SCOPE)
endfunction()

# - call with parameters:
#   - create a custom simlation target: several steps to create a name etc.
#   -

function (simulate)

# #######################################################################
# Parse arguments
# #######################################################################

set(options "")
set(one_value_args
    SIMULATION_NAME #STEM name of the simulated files
    REFERENCE_FILE  #option: -ir
    READ_COUNT      #option: -n
    READ_LENGTH     #option: --illumina-read-length
    INDEL_PROB      #option: --illumina-prob-insert and --illumina-prob-deletion
    SNP_PROB        #option: --illumina-prob-mismatch
    THREAD_COUNT    #option: --num-threads
)
set(multi_value_args "")
cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

if ("${ARG_SIMULATION_NAME}" STREQUAL "")
    message (FATAL_ERROR "Please specify a name for the simulation target [SIMULATION_NAME]")
endif()

if ("${ARG_REFERENCE_FILE}" STREQUAL "")
    message (FATAL_ERROR "Please specify the name of the reference data source target [REFERENCE_FILE]")
endif()

if ("${ARG_READ_COUNT}" STREQUAL "")
    message (FATAL_ERROR "Please specify the number of reads to simulate [READ_COUNT]")
endif()

# #######################################################################
# Prepare reference file for simulation
# #######################################################################

# Prepare reference_target for the corresponding data source
datasource_target("${ARG_REFERENCE_FILE}" reference_target)

# Prepare reference file for simulation
ExternalProject_Get_property("${reference_target}" DOWNLOAD_NAME)
ExternalProject_Get_property("${reference_target}" INSTALL_DIR)
ExternalProject_Get_property("${reference_target}" URL_HASH)
string(REGEX REPLACE "(.*)\.gz$" "\\1" decompressed_file_name "${DOWNLOAD_NAME}")

set(reference_name "${DOWNLOAD_NAME}")
set(reference_hash "${URL_HASH}")
set(installed_file "${INSTALL_DIR}/${reference_name}")

# Prepare extraction command
string(REGEX MATCH ".*tar.gz$" has_tar_gz "${DOWNLOAD_NAME}")
string(REGEX MATCH ".*gz$" has_gz "${DOWNLOAD_NAME}")
set (needs_gunzip "$<AND:$<NOT:$<BOOL:${has_tar_gz}>>,$<BOOL:${has_gz}>>")

find_program(GUNZIP_COMMAND "gunzip")
set (has_gunzip "$<BOOL:${GUNZIP_COMMAND}>")

set (EXTRACT_GUNZIP_COMMAND "${GUNZIP_COMMAND}")
set (EXTRACT_TAR_COMMAND "${CMAKE_COMMAND}" "-E" "tar" "-xf")
set (use_gunzip "$<AND:${needs_gunzip},${has_gunzip}>")
set (EXTRACT_COMMAND "$<IF:${use_gunzip},${EXTRACT_GUNZIP_COMMAND},${EXTRACT_TAR_COMMAND}>")

# #######################################################################
# Setup simulation target
# #######################################################################

datasource_target("${ARG_SIMULATION_NAME}" simulation_target)
# string(TOLOWER "simsource--${ARG_SIMULATION_NAME}" simulation_target)

# Create local working directory for simulated files
set(SIMULATION_WORKING_DIR "${DATA_ROOT_DIR}/_simulation_sources")
file(MAKE_DIRECTORY "${SIMULATION_WORKING_DIR}")

# Prepare read file
set (read_file "${ARG_SIMULATION_NAME}.fq")
set (read_alignment_file "${ARG_SIMULATION_NAME}.bam")
# set (error_profile_file "${ARG_SIMULATION_NAME}_error_profile.txt")

# Prepare read length
set(read_length "$<IF:$<BOOL:${ARG_READ_LENGTH}>,${ARG_READ_LENGTH},100>")

# Prepare number of threads

set(thread_count "$<IF:$<BOOL:${ARG_THREAD_COUNT}>,${ARG_THREAD_COUNT},1>")

# Write error profile file into working directory
if ("${ARG_INDEL_PROB}" STREQUAL "")
    set (indel_prob "0.0001")
else()
    set (indel_prob "${ARG_INDEL_PROB}")
endif()

if ("${ARG_SNP_PROB}" STREQUAL "")
    set (snp_prob "0.004")
else()
    set (snp_prob "${ARG_SNP_PROB}")
endif()

# math(EXPR del_prob "${indel_prob} / 2.0" OUTPUT_FORMAT DECIMAL)
floatexpr("${indel_prob} / 2.0" del_prob)
floatexpr("${ARG_INDEL_PROB} - ${del_prob}" ins_prob)

# Prepare mason command
set(MASON_COMMAND "mason_simulator")
set(MASON_GLOBAL_OPTIONS "--illumina-read-length" "${read_length}" "--embed-read-info" "--num-threads" "${thread_count}")
set(MASON_SIM_OPTIONS "--illumina-prob-insert" "${ins_prob}" "--illumina-prob-deletion" "${del_prob}" "--illumina-prob-mismatch" "${snp_prob}")
set(MASON_ARGS "-n" "${ARG_READ_COUNT}" "-o" "${read_file}" "-oa" "${read_alignment_file}")

ExternalProject_Add(
    ${simulation_target}
    DEPENDS           "${reference_target}"
    DOWNLOAD_COMMAND  ${CMAKE_COMMAND} -E copy ${installed_file} <SOURCE_DIR>/
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> ${EXTRACT_COMMAND} ${reference_name}
    BUILD_COMMAND     ${CMAKE_COMMAND} -E echo Simulating reads: "${MASON_COMMAND} ${MASON_GLOBAL_OPTIONS} ${MASON_SIM_OPTIONS} -ir <SOURCE_DIR>/${decompressed_file_name} ${MASON_ARGS}"
    COMMAND           ${MASON_COMMAND} ${MASON_GLOBAL_OPTIONS} ${MASON_SIM_OPTIONS} -ir <SOURCE_DIR>/${decompressed_file_name} ${MASON_ARGS}
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E create_symlink <BINARY_DIR>/${read_file} <INSTALL_DIR>/${read_file}
    COMMAND           ${CMAKE_COMMAND} -E create_symlink <BINARY_DIR>/${read_alignment_file} <INSTALL_DIR>/${read_alignment_file}
    TEST_COMMAND      ""
    PREFIX            "${SIMULATION_WORKING_DIR}"
    INSTALL_DIR       "${DATA_DIR}"
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    DOWNLOAD_NO_EXTRACT TRUE
    ${ARG_UNPARSED_ARGUMENTS}
)
endfunction()


# create a simulation target
    # search for mason application
    # get list of parameters: read count | read length | snp error rate | indel error rate | reference sequence
    # invoke mason list
    # manage output fasta file as byproduct to depend on



# create a dependent target to get the file from somewhere
