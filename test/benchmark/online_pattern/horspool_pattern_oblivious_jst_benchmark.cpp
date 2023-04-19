// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <numeric>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>

#include <libjst/matcher/horspool_matcher.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>


template <typename tree_t>
inline size_t total_bytes(tree_t const & tree) noexcept {
    return libjst::stats(tree).symbol_count;
}

template <typename ...args_t>
static void bench(benchmark::State & state, args_t && ...args)
{
    using sequence_t = jstmap::reference_t;
    auto [jst_file, needle_file] = std::tuple{args...};
    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t needle = jstmap::load_queries(needle_file)[0].sequence();

    libjst::horspool_matcher matcher{needle};

    auto search_tree = libjst::make_volatile(rcs_store) | libjst::labelled()
                                                        | libjst::coloured()
                                                        | libjst::trim(libjst::window_size(matcher) - 1)
                                                        | libjst::prune()
                                                        | libjst::left_extend(libjst::window_size(matcher) - 1)
                                                        | libjst::merge();

    size_t hit_count{};
    for (auto _ : state)
    {
        hit_count = 0;
        libjst::tree_traverser_base oblivious_path{search_tree};
        for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
            auto && cargo = *it;
            matcher(cargo.sequence(), [&] (auto const &) { ++hit_count; });
        }
    }

    size_t processed_bytes = total_bytes(search_tree);
    state.counters["bytes"] = processed_bytes;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(processed_bytes);
    state.counters["#hits"] = hit_count;
}

BENCHMARK_CAPTURE(bench,
                  online_pattern_plain_needle32,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle32.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  online_pattern_plain_needle64,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle64.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  online_pattern_plain_needle128,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle128.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  online_pattern_plain_needle256,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle256.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_MAIN();
