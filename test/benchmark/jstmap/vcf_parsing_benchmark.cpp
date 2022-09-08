// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <vector>

#include <seqan/vcf_io.h>

#include <seqan3/test/performance/units.hpp>
#include <seqan3/std/filesystem>

template <typename ...args_t>
static void parse_vcf(benchmark::State & state, args_t && ...args)
{
    std::filesystem::path vcf_file_path{get<0>(std::tuple{args...})};

    seqan::VcfFileIn vcf_file{vcf_file_path.c_str()};
    seqan::VcfHeader vcf_header{};
    seqan::readHeader(vcf_header, vcf_file);

    seqan::VcfRecord _record{};

    std::vector<size_t> records_per_chr{};
    records_per_chr.resize(1000, 0);
    for (auto _ : state){
        while (!seqan::atEnd(vcf_file))
        {
            seqan::readRecord(_record, vcf_file);
            ++records_per_chr[_record.rID];
        }
    }

    size_t vcf_file_size = file_size(vcf_file_path);
    size_t record_count{};
    std::ranges::for_each(records_per_chr, [&] (size_t hits) { record_count += hits; });

    state.counters["file_size"] = vcf_file_size;
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(vcf_file_size);
    state.counters["#records"] = record_count;
}

BENCHMARK_CAPTURE(parse_vcf, vcf, DATADIR"1KGP.chr22.vcf");

BENCHMARK_MAIN();
