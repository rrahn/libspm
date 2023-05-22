// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>

namespace just::bench
{

struct benchmark_configuration
{
    std::filesystem::path jst_file;
    std::filesystem::path needle_file;
    std::filesystem::path jst_ibf_file;
};

inline constexpr auto chr22_needle32_ibf32768 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle32.fa"},
        .jst_ibf_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf"}
    };
};

inline constexpr auto chr22_needle64_ibf32768 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle64.fa"},
        .jst_ibf_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf"}
    };
};

inline constexpr auto chr22_needle128_ibf32768 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle128.fa"},
        .jst_ibf_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf"}
    };
};

inline constexpr auto chr22_needle256_ibf32768 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle256.fa"},
        .jst_ibf_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst.c32768.k21.ibf"}
    };
};
}  // namespace just::bench
