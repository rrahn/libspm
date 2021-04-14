// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/index/vcf_parser.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

auto create_jst_from_vcf(std::filesystem::path reference_file, std::filesystem::path vcf_file)
{
    return jstmap::construct_jst_from_vcf(reference_file, vcf_file);
}

template <typename ...args_t>
static void jst_search_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [reference_file, vcf_file, query_file] = std::tuple{args...};

    auto jst = create_jst_from_vcf(reference_file, vcf_file);
    auto queries = jstmap::load_queries(query_file);
    size_t total_bytes{};
    for (auto const & query : queries)
        total_bytes += std::ranges::size(query);

    size_t hit_count{};
    for (auto _ : state)
        hit_count += jstmap::search_queries(jst, queries).size();

    benchmark::DoNotOptimize(hit_count);
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(total_bytes);
}

// Register the function as a benchmark
BENCHMARK_CAPTURE(jst_search_benchmark,
                  vcf_indel_test,
                  DATADIR"sim_ref_10Kb.fasta.gz",
                  DATADIR"sim_ref_10Kb_SNP_INDELs.vcf",
                  DATADIR"sim_reads_ref1x10.fa");

BENCHMARK_MAIN();
