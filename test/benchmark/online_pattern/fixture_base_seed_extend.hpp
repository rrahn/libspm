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

#include <benchmark/benchmark.h>

#include <omp.h>

#include <functional>
#include <ranges>
#include <stack>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/match_position.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/filter_queries.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include "fixture_config.hpp"

namespace just::bench
{

    template <auto _config_generator>
    struct capture {
        auto fixture() -> std::invoke_result_t<decltype(*_config_generator)>
        {
            return std::invoke(*_config_generator);
        }
    };

    template <typename capture_t>
    class fixture_base_seed_extend : public ::benchmark::Fixture,
                                     protected capture_t
    {
    protected:
        using sequence_t = jstmap::reference_t;
        using queries_t = std::vector<jstmap::sequence_record_t>;

        size_t processed_bytes;
    private:

        jstmap::rcs_store_t _rcs_store{};
        queries_t _queries{};

    public:

        fixture_base_seed_extend() = default;
        virtual ~fixture_base_seed_extend() = default;

        virtual void SetUp(::benchmark::State const &) override {
            auto [jst_file, needle_file] = this->fixture();

            _rcs_store = jstmap::load_jst(jst_file);
            _queries = jstmap::load_queries(needle_file);
        }

        virtual void TearDown(::benchmark::State & state) override {
            state.counters["bytes"] = processed_bytes;
            state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(state.counters["bytes"]);
        }

        sequence_t const & needle() const noexcept {
            return _queries[0].sequence();
        }

        jstmap::rcs_store_t const & store() const noexcept {
            return _rcs_store;
        }

        auto queries() const noexcept {
            return _queries | std::views::transform([] (auto const & record) -> sequence_t const & {
                return record.sequence();
            });
        }

        constexpr size_t get_chunk_size(size_t const thread_count) noexcept {
            return (std::ranges::size(store().source()) + thread_count - 1) / thread_count;
        }

        constexpr float to_error_rate(int32_t error_count) noexcept {
            return static_cast<float>(error_count)/100.0 + 0.0001;
        }

        inline size_t total_bytes() noexcept {
            size_t window_size = needle().size();

            auto tree = store() | libjst::make_volatile()
                                | libjst::labelled()
                                | libjst::coloured()
                                | libjst::trim(window_size - 1)
                                | libjst::prune()
                                | libjst::left_extend(window_size - 1)
                                | libjst::merge();

            return libjst::stats(tree).symbol_count * std::ranges::size(queries());
        }

        template <typename runner_creator_t>
        void run(::benchmark::State & state, runner_creator_t && make_runner)
        {
            auto trees = libjst::chunk(store(), get_chunk_size(state.range(0)));
            int32_t hit_count{};
            for (auto _ : state)
            {
                benchmark::DoNotOptimize(hit_count = execute(trees, make_runner, queries(), state.range(0)));
                benchmark::ClobberMemory();
            }
        }

        template <typename trees_t, typename runner_creator_t, typename search_queries_t>
        static int32_t execute(trees_t && trees,
                               runner_creator_t && make_runner,
                               search_queries_t && queries,
                               size_t const thread_count) {
            int32_t hit_count = 0;
            std::ptrdiff_t const chunk_count = std::ranges::ssize(trees);

            #pragma omp parallel for num_threads(thread_count), shared(trees), firstprivate(make_runner, queries), schedule(static), reduction(+:hit_count)
            for (std::ptrdiff_t chunk = 0; chunk < chunk_count; ++chunk) {
                auto runner = make_runner(trees[chunk], queries);
                runner([&] (std::ptrdiff_t, jstmap::match_position) { ++hit_count; });
            }
            return hit_count;
        }
    };
}  // namespace just::bench

