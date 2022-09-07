// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <numeric>

#include <libjst/search/naive_search.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>

#include <libjst/matcher/horspool_matcher.hpp>
#include <libjst/matcher/shiftor_matcher.hpp>
#include <libjst/matcher/shiftor_matcher_restorable.hpp>
#include <libjst/search/polymorphic_sequence_searcher.hpp>

#include "benchmark_utility.hpp"

template <typename rcs_store_t>
inline size_t total_bytes(rcs_store_t const & rcs_store) noexcept {
    return std::ranges::size(rcs_store.source()) * rcs_store.size();
}

template <typename matcher_t, typename rcs_store_t>
static void run(benchmark::State & state, matcher_t && matcher, rcs_store_t const & rcs_store) {
    libjst::polymorphic_sequence_searcher searcher{rcs_store};

    size_t hit_count{};
    for (auto _ : state) {
        searcher(matcher, [&] (...) { ++hit_count; });
    }

    state.counters["bytes"] = total_bytes(rcs_store);
    state.counters["bytes_per_second"] =  seqan3::test::bytes_per_second(total_bytes(rcs_store));
    state.counters["#hits"] = hit_count;
}

template <typename ...args_t>
static void naive_search_benchmark(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};
    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t query{sample_query(rcs_store.source(), state.range(0))};

    size_t sequence_count = rcs_store.size();
    size_t batch_size = 16;
    size_t total_runs = sequence_count / batch_size;
    std::vector<sequence_t> batch{batch_size, rcs_store.source()};

    seqan::Pattern<sequence_t, seqan::Horspool> pattern{query};

    size_t hit_count{};
    for (auto _ : state)
    {
        for (size_t run = 0; run < total_runs; ++run) {
            std::ranges::for_each(batch, [&] (auto const & seq) {
                using range_t = std::remove_reference_t<decltype(seq)>;
                seqan::Finder<range_t> finder{seq};
                while (seqan::find(finder, pattern)) {
                    ++hit_count;
                }
            });
        }
    }
    state.counters["bytes"] = total_bytes(rcs_store);
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(total_bytes(rcs_store));
    state.counters["#hits"] = hit_count;
}

template <typename ...args_t>
static void jst_search_benchmark_horspool(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};
    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t query{sample_query(rcs_store.source(), state.range(0))};

    libjst::horspool_matcher matcher{query};
    run(state, matcher, rcs_store);
}

template <typename ...args_t>
static void jst_search_benchmark_shiftor(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};

    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t query{sample_query(rcs_store.source(), state.range(0))};

    libjst::shiftor_matcher matcher{query};
    run(state, matcher, rcs_store);
}

template <typename ...args_t>
static void jst_search_benchmark_restorable_shiftor(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};

    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t query{sample_query(rcs_store.source(), state.range(0))};

    libjst::restorable_shiftor_matcher matcher{query};
    run(state, matcher, rcs_store);
}
// Register the function as a benchmark

BENCHMARK_CAPTURE(jst_search_benchmark_horspool,
                  vcf_indel_test,
                  "/home/rahn/workspace/data/jstmap/new/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst")->Arg(30)->Arg(60)->Arg(120);
BENCHMARK_CAPTURE(jst_search_benchmark_shiftor,
                  vcf_indel_test,
                  "/home/rahn/workspace/data/jstmap/new/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst")->Arg(30)->Arg(60)->Arg(120);
BENCHMARK_CAPTURE(jst_search_benchmark_restorable_shiftor,
                  vcf_indel_test,
                  "/home/rahn/workspace/data/jstmap/new/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst")->Arg(30)->Arg(60)->Arg(120);
BENCHMARK_CAPTURE(naive_search_benchmark,
                  vcf_indel_test,
                  "/home/rahn/workspace/data/jstmap/new/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst")->Arg(30)->Arg(60)->Arg(120);

BENCHMARK_MAIN();
