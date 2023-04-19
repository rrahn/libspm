// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <numeric>
#include <stack>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>

#include <libjst/matcher/shiftor_matcher_restorable.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

template <typename matcher_t>
class state_manager {
private:

    using state_t = libjst::matcher_state_t<matcher_t>;
    using state_stack_t = std::stack<state_t>;

    matcher_t & _matcher;
    state_stack_t _states{};

public:

    constexpr explicit state_manager(matcher_t & matcher) noexcept :
        _matcher{matcher},
        _states{}
    {}

    constexpr void notify_push() {
        _states.push(_matcher.capture());
    }

    constexpr void notify_pop() {
        assert(!_states.empty());
        _matcher.restore(_states.top());
        _states.pop();
    }
};

// template <typename matcher_t>
// state_manager(matcher_t const &) -> state_capture_traverser::state_manager<matcher_t>;

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

    libjst::restorable_shiftor_matcher matcher{needle};

    auto search_tree = libjst::make_volatile(rcs_store) | libjst::labelled()
                                                        | libjst::coloured()
                                                        | libjst::trim(libjst::window_size(matcher) - 1)
                                                        | libjst::prune()
                                                        | libjst::merge();

    size_t hit_count{};
    for (auto _ : state)
    {
        hit_count = 0;
        libjst::tree_traverser_base resumable_path{search_tree};
        state_manager manager{matcher};
        resumable_path.subscribe(manager);
        for (auto it = resumable_path.begin(); it != resumable_path.end(); ++it) {
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
                  shiftor_pattern_resumable_jst_needle32,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle32.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  shiftor_pattern_resumable_jst_needle64,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle64.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  shiftor_pattern_resumable_jst_needle128,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle128.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(bench,
                  shiftor_pattern_resumable_jst_needle256,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle256.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_MAIN();
