// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <libjst/search/pigeonhole_filter.hpp>
#include <libjst/search/state_manager_stack.hpp>

#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

#include "benchmark_utility.hpp"

template <typename jst_t>
auto sample_reads(jst_t const & jst, size_t const read_count, size_t const read_size)
{
    auto enumerator = jst.context_enumerator(read_size);
    size_t const sample_rate = std::ceil(jst.reference().size() / read_count);
    size_t counter = 0;
    seqan::StringSet<sequence_t> sampled_reads{};

    for (auto it = enumerator.begin(); it != enumerator.end(); ++it)
    {
        if (++counter % sample_rate == 0)
        {
            auto read = *it;
            seqan::appendValue(sampled_reads, sequence_t{read.begin(), read.end()});
        }
    }
    return sampled_reads;
}

template <typename ...args_t>
static void pigeonhole_filter_bench(benchmark::State & state, args_t && ...args)
{
    auto [jst_file] = std::tuple{args...};
    size_t const read_count = state.range(0);
    size_t const read_length = state.range(1);
    auto [jst, partitioned_jst_handle] = jstmap::load_jst(jst_file);
    auto sampled_reads = sample_reads(jst, read_count, read_length);

    // Prepare pigeonhole filter.
    using state_t = typename decltype(libjst::pigeonhole_filter{sampled_reads})::state_type;
    libjst::pigeonhole_filter filter{sampled_reads, 0.0, libjst::search_state_manager_stack<state_t>{}};

    size_t const fragment_size = filter.qgram_size();

    size_t hit_count{};
    for (auto _ : state)
    {
        hit_count = 0;
        auto jst_range_agent = jst.range_agent(fragment_size, filter.state_manager()); // already pushing a branch.
        filter(jst_range_agent, [&] (auto const &, auto const &) { ++hit_count; });
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(jst.total_symbol_count());
    state.counters["#hits"] = hit_count;
}

static void CustomArguments(benchmark::internal::Benchmark* b)
{
    std::array range1{50, 100, 500, 1000};
    std::array range2{100, 150, 200, 250, 300};

    for (size_t i = 0; i < range1.size(); ++i)
        for (size_t j = 0; j < range2.size(); ++j)
            b->Args({range1[i], range2[j]});
}

BENCHMARK_CAPTURE(pigeonhole_filter_bench,
                  pigeonhole_filter,
                  DATADIR"sim_ref_10Kb_SNP_INDELs")->Apply(CustomArguments);

BENCHMARK_MAIN();
