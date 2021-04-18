// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

#include "benchmark_utility.hpp"

template <typename ...args_t>
static void naive_search_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [haplotypes_file] = std::tuple{args...};

    auto haplotypes = jstmap::load_queries(haplotypes_file);
    sequence_t needle = generate_query(state.range(0));
    size_t const total_bytes = std::ranges::size(needle);

    size_t hit_count{};
    for (auto _ : state)
    {
        naive_traversal(haplotypes, [&] <std::ranges::random_access_range sequence_t>(sequence_t && sequence)
        {
            auto found_range = std::ranges::search(sequence, needle);
            hit_count += !std::ranges::empty(found_range);
            return std::ranges::subrange{std::ranges::next(found_range.begin(), 1, std::ranges::end(sequence)),
                                         std::ranges::end(sequence)};
        });
    }

    benchmark::DoNotOptimize(hit_count);
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(total_bytes);
    state.counters["#hits"] = hit_count;
}

template <typename ...args_t>
static void jst_search_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};
    auto jst = jstmap::load_jst(jst_file);
    std::vector<sequence_t> query{generate_query(state.range(0))};
    size_t const total_bytes = std::ranges::size(query[0]);

    size_t hit_count{};
    for (auto _ : state)
        hit_count += jstmap::search_queries(jst, query).size();

    benchmark::DoNotOptimize(hit_count);
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(total_bytes);
    state.counters["#hits"] = hit_count;
}

// Register the function as a benchmark
BENCHMARK_CAPTURE(naive_search_benchmark,
                  vcf_indel_test,
                  DATADIR"sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz")->Arg(64)->Arg(100)->Arg(150);

BENCHMARK_CAPTURE(jst_search_benchmark,
                  vcf_indel_test,
                  DATADIR"sim_ref_10Kb_SNP_INDELs.jst")->Arg(64)->Arg(100)->Arg(150);

BENCHMARK_MAIN();
