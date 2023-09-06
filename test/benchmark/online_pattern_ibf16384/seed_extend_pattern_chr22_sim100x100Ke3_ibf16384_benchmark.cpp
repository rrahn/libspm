// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <jstmap/search/bucket.hpp>
#include <jstmap/search/bucket_searcher.hpp>

#include "fixture_base_seed_extend.hpp"

namespace just::bench {

BENCHMARK_TEMPLATE_DEFINE_F(fixture_base_seed_extend, seed_extend, capture<&chr22_sim100x100Ke3_ibf16384>)(benchmark::State& state) {
    run(state,
        [&] (auto && tree, auto && needles) {
            jstmap::bucket tmp{.base_tree = std::move(tree), .needle_list = std::move(needles)};
            jstmap::bucket_searcher searcher{std::move(tmp), to_error_rate(state.range(1))};
            return searcher;
        });
    processed_bytes = total_bytes();
}

BENCHMARK_REGISTER_F(fixture_base_seed_extend, seed_extend)
    ->ArgsProduct({
        benchmark::CreateRange(1, std::thread::hardware_concurrency(), 2),
        {0, 1, 2, 3}
    })
    ->UseRealTime();
} // namespace just::bench

BENCHMARK_MAIN();
