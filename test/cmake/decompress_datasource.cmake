include (ExternalProject)

function (decompress)

# #######################################################################
# Parse arguments
# #######################################################################

set(options "")
set(one_value_args
    FILE            # name of the file to decompress including the compression ending
)
set(multi_value_args "")
cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

if ("${ARG_FILE}" STREQUAL "")
    message (FATAL_ERROR "Please specify the target file to decompress, e.g. seq.fa.gz [FILE]")
endif()

# #######################################################################
# Get datasource file target
# #######################################################################

# Prepare reference_target for the corresponding data source
datasource_target("${ARG_FILE}" file_target)

# Prepare reference file for simulation
ExternalProject_Get_property("${file_target}" DOWNLOAD_NAME)
ExternalProject_Get_property("${file_target}" INSTALL_DIR)
string(REGEX REPLACE "(.*)\.gz$" "\\1" decompressed_file_name "${DOWNLOAD_NAME}")

set(file_name "${DOWNLOAD_NAME}")
set(file_hash "${URL_HASH}")
set(installed_file "${INSTALL_DIR}/${file_name}")

message(STATUS "installed_file: ${installed_file}")
message(STATUS "decompressed_file_name: ${decompressed_file_name}")

# Prepare extraction command
string(REGEX MATCH ".*tar.gz$" has_tar_gz "${DOWNLOAD_NAME}")
string(REGEX MATCH ".*gz$" has_gz "${DOWNLOAD_NAME}")
set (needs_gunzip "$<AND:$<NOT:$<BOOL:${has_tar_gz}>>,$<BOOL:${has_gz}>>")

find_program(GUNZIP_COMMAND "gunzip")
set (has_gunzip "$<BOOL:${GUNZIP_COMMAND}>")

set (EXTRACT_GUNZIP_COMMAND "${GUNZIP_COMMAND}" "-f")
set (EXTRACT_TAR_COMMAND "${CMAKE_COMMAND}" "-E" "tar" "-xf")
set (use_gunzip "$<AND:${needs_gunzip},${has_gunzip}>")
set (EXTRACT_COMMAND "$<IF:${use_gunzip},${EXTRACT_GUNZIP_COMMAND},${EXTRACT_TAR_COMMAND}>")

datasource_target("${decompressed_file_name}" decompressed_file_target)

message(STATUS "decompressed_file_target: ${decompressed_file_target}")
message(STATUS "file_target: ${file_target}")
message(STATUS "file_name: ${file_name}")

ExternalProject_Add(
    ${decompressed_file_target}
    DEPENDS           "${file_target}"
    DOWNLOAD_COMMAND  ${CMAKE_COMMAND} -E copy ${installed_file} <SOURCE_DIR>/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> ${EXTRACT_COMMAND} ${file_name}
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/${decompressed_file_name} <INSTALL_DIR>/${decompressed_file_name}
    TEST_COMMAND      ""
    PREFIX            "${DATA_ROOT_DIR}/_datasources"
    INSTALL_DIR       "${DATA_DIR}"
    LOG_DOWNLOAD     TRUE
    LOG_BUILD        TRUE
    DOWNLOAD_NO_EXTRACT TRUE
    ${ARG_UNPARSED_ARGUMENTS}
)
endfunction()
