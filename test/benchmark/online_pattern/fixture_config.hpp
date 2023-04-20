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
};

inline constexpr auto chr22_needle32 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle32.fa"}
    };
};

inline constexpr auto chr22_needle64 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle64.fa"}
    };
};

inline constexpr auto chr22_needle128 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle128.fa"}
    };
};

inline constexpr auto chr22_needle256 = [] () {
    return benchmark_configuration{
        .jst_file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"},
        .needle_file{DATADIR"needle256.fa"}
    };
};
}  // namespace just::bench
