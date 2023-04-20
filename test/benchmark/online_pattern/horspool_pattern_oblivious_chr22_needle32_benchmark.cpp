// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <libjst/matcher/horspool_matcher.hpp>

#include "fixture_oblivious_pattern.hpp"

namespace just::bench {

BENCHMARK_TEMPLATE_F(fixture_oblivious_pattern, horspool, capture<&chr22_needle32>)(benchmark::State& state) {
    run(state, libjst::horspool_matcher(needle()));
}
} // namespace just::bench

BENCHMARK_MAIN();
