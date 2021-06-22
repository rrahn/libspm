// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <jstmap/search/load_queries.hpp>

#include "benchmark_utility.hpp"

template <typename ...args_t>
static void naive_context_enumerator_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [haplotypes_file] = std::tuple{args...};

    auto haplotypes = jstmap::load_queries(haplotypes_file);

    int32_t context_size = state.range(0);
    size_t context_count{};
    for (auto _ : state)
    {
        naive_traversal(haplotypes, [&] <std::ranges::random_access_range sequence_t>(sequence_t && sequence)
        {
            auto next_context_begin = std::ranges::next(std::ranges::begin(sequence), 1, std::ranges::end(sequence));
            if (std::ranges::distance(next_context_begin, std::ranges::end(sequence)) < context_size)
                next_context_begin = std::ranges::next(next_context_begin, std::ranges::end(sequence));

            ++context_count;
            return std::ranges::subrange(next_context_begin, std::ranges::end(sequence));
        });
    }

    benchmark::DoNotOptimize(context_count);
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(context_size);
}

template <typename ...args_t>
static void jst_context_enumerator_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};
    auto jst_handle = jstmap::load_jst(jst_file);

    size_t context_size = state.range(0);
    size_t context_count{};
    auto context_enumerator = jst_handle.first.context_enumerator(context_size);

    for (auto _ : state)
        std::ranges::for_each(context_enumerator, [&] (auto &&) { ++context_count; });

    benchmark::DoNotOptimize(context_count);
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(context_size);
}

// Register the function as a benchmark
BENCHMARK_CAPTURE(naive_context_enumerator_benchmark,
                  vcf_indel_test,
                  DATADIR"sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz")->Arg(64)->Arg(100)->Arg(150)->Arg(200);

BENCHMARK_CAPTURE(jst_context_enumerator_benchmark,
                  vcf_indel_test,
                  DATADIR"sim_ref_10Kb_SNP_INDELs.jst")->Arg(64)->Arg(100)->Arg(150)->Arg(200);

BENCHMARK_MAIN();
