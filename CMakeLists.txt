cmake_minimum_required (VERSION 3.14)

## CUSTOMISE

# Define the application name and version.
project (jstmap VERSION 0.0.1 LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 REQUIRED HINTS lib/seqan3/build_system)

# Dependency: SeqAn.
set (SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include")
find_package (SeqAn REQUIRED HINTS "lib/seqan/util/cmake" NO_DEFAULT_PATH)

separate_arguments (SEQAN_CXX_FLAGS_LIST UNIX_COMMAND "${SEQAN_CXX_FLAGS}")

add_library (seqan2_seqan2 INTERFACE IMPORTED)
target_compile_definitions (seqan2_seqan2 INTERFACE ${SEQAN_DEFINITIONS})
target_compile_options (seqan2_seqan2 INTERFACE "${SEQAN_CXX_FLAGS_LIST}")
target_link_libraries (seqan2_seqan2 INTERFACE ${SEQAN_LIBRARIES})
target_include_directories (seqan2_seqan2 INTERFACE ${SEQAN_INCLUDE_DIRS})
add_library (seqan2::seqan2 ALIAS seqan2_seqan2)

get_target_property(seqan2_definitions seqan2_seqan2 COMPILE_DEFINITIONS)
get_target_property(seqan2_options seqan2_seqan2 COMPILE_OPTIONS)
get_target_property(seqan2_libs seqan2_seqan2 LINK_LIBRARIES)
get_target_property(seqan2_dirs seqan2_seqan2 INCLUDE_DIRECTORIES)

message(STATUS "seqan2_lib: ${SEQAN_FOUND}")
message(STATUS "seqan2_definitions: ${SEQAN_DEFINITIONS}")
message(STATUS "seqan2_options: ${SEQAN_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}")
message(STATUS "seqan2_libs: ${SEQAN_LIBRARIES} | ${SEQAN3_LIBRARIES}")
message(STATUS "seqan2_dirs: ${SEQAN_INCLUDE_DIRS}")

# Dependency: OpenMP
find_package(OpenMP REQUIRED)

# Dependency: PThreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

## APPLICATION

add_subdirectory (libcontrib)
add_subdirectory (libjst)
add_subdirectory (jstmap-global)
add_subdirectory (jstmap-create)
add_subdirectory (jstmap-index)
add_subdirectory (jstmap-search)
add_subdirectory (jstmap-simulate)
add_subdirectory (jstmap-view)
add_subdirectory (jstmap)

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

include (test/jstmap_test.cmake)

add_subdirectory (test/api EXCLUDE_FROM_ALL)
add_subdirectory (test/cli EXCLUDE_FROM_ALL)

add_subdirectory (test/benchmark EXCLUDE_FROM_ALL)
