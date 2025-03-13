# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# seqan3
set (LIBSPM_SEQAN3_VERSION 7e0d88d15fc82b8b8a5548f6eebea8602faf6446)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${LIBSPM_SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# seqan2
set (LIBSPM_SEQAN2_VERSION 7a8ef3cef61c57a4098018c3906daf9802cbfc4e)
CPMDeclarePackage (seqan2
                   NAME seqan2
                   GIT_TAG ${LIBSPM_SEQAN2_VERSION}
                   GITHUB_REPOSITORY rrahn/seqan
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "SEQAN_BUILD_SYSTEM DEVELOP" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# googletest
set (LIBSPM_GOOGLETEST_VERSION 1.16.0)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${LIBSPM_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_CXX_STANDARD 20"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

#googlebenchmark
set (LIBSPM_GOOGLEBENCHMARK_VERSION 1.9.1)
CPMDeclarePackage (googlebenchmark
                   NAME benchmark
                   VERSION ${LIBSPM_GOOGLEBENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF" "BENCHMARK_ENABLE_INSTALL OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING" "CMAKE_CXX_STANDARD 20"
)

# catch2
set (LIBSPM_CATCH2_VERSION 3.8.0)
CPMDeclarePackage (catch2
                   NAME Catch2
                   VERSION ${LIBSPM_CATCH2_VERSION}
                   GITHUB_REPOSITORY catchorg/Catch2
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (LIBSPM_USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${LIBSPM_USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
