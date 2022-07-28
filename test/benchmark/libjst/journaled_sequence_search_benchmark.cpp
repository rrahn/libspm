// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <libjst/journal.hpp>
#include <libjst/search/horspool_search.hpp>
#include <libjst/search/myers_search.hpp>
#include <libjst/search/naive_search.hpp>
#include <libjst/search/shift_or_search.hpp>

#include "sequence_variant_simulation.hpp"

static void benchmark_args(benchmark::internal::Benchmark* b) {

    for (int i = 7; i <= 20; i += 1) {
        // pattern sizes:
        int32_t const haystack_size = 1u << i;
        for (int k = 5; k <= 7; k += 1 ) {
            int32_t const needle_size = 1u << k;
            b->Args({haystack_size, needle_size, 0}); // 0 variance
            for (float j = 0.01; j <= 0.1; j += 0.01) {
                int64_t errors = std::ceil(haystack_size * j);
                b->Args({haystack_size, needle_size, errors});
            }
        }
    }
}

template<typename container_t, typename searcher_t>
class benchmark_search : public benchmark::Fixture
{
    using base_sequence_t = std::vector<char>;
    using variants_t = decltype(generate_variants(0, 0));

    base_sequence_t base_sequence{};
    variants_t sequence_variants{};
    container_t source{};
    base_sequence_t needle{};

public:

    benchmark_search() {
    }

    virtual void SetUp(benchmark::State & state) override {
        // Create source sequence!
        size_t const sequence_size = state.range(0);
        base_sequence.resize(sequence_size, 'A');
        sequence_variants = generate_variants(base_sequence.size(), state.range(2));
        source = generate_sequence<container_t>(base_sequence, sequence_variants);

        needle.resize(state.range(1));

        auto && haystack = get_haystack();

        // Create needle!
        std::mt19937 rng{sequence_size};
        std::uniform_int_distribution<size_t> needle_begin_dist{0, haystack.size() - needle.size()};
        size_t needle_begin = needle_begin_dist(rng);
        std::ranges::copy(haystack.begin() + needle_begin, haystack.begin() + needle_begin + needle.size(),
                          needle.begin());
    }

    decltype(auto) get_haystack() {
        if constexpr (std::same_as<container_t, base_sequence_t>) {
            return source;
        } else {
            return source.sequence();
        }
    }

    base_sequence_t get_needle() {
        return needle;
    }

    // template <typename searcher_t, typename haystack_t>
    void run(benchmark::State & state) { // , searcher_t && searcher, haystack_t && haystack) {
        auto haystack = get_haystack();
        auto needle = get_needle();

        size_t hit_count{};

        for (auto _ : state) {
            state.PauseTiming();
            searcher_t searcher{needle};
            state.ResumeTiming();
            searcher(haystack, [&] ([[maybe_unused]] auto && haystack_it) { ++hit_count; });
        }

        state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(haystack.size());
        state.counters["hit_count"] = hit_count;
    }
};

// Benchmark types
using source_t = std::vector<char>;
using journal_t = libjst::journal<uint32_t, source_t &>;

// ----------------------------------------------------------------------------
// Benchmark seach naive
// ----------------------------------------------------------------------------

using naive_t = decltype(libjst::naive_pattern_searcher{std::declval<source_t &>()});

BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, naive_vector, source_t, naive_t) (benchmark::State & state) { this->run(state); }
BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, naive_journal, journal_t, naive_t)(benchmark::State & state) { this->run(state); }

BENCHMARK_REGISTER_F(benchmark_search, naive_vector)->Apply(benchmark_args);
BENCHMARK_REGISTER_F(benchmark_search, naive_journal)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark seach horspool
// ----------------------------------------------------------------------------

using horspool_t = decltype(libjst::horspool_pattern_searcher{std::declval<source_t &>()});

BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, horspool_vector, source_t, horspool_t) (benchmark::State & state) { this->run(state); }
BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, horspool_journal, journal_t, horspool_t)(benchmark::State & state) { this->run(state); }

BENCHMARK_REGISTER_F(benchmark_search, horspool_vector)->Apply(benchmark_args);
BENCHMARK_REGISTER_F(benchmark_search, horspool_journal)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark seach shift-or/bitap
// ----------------------------------------------------------------------------

using bitap_t = decltype(libjst::shift_or_pattern_searcher{std::declval<source_t &>()});

BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, bitap_vector, source_t, bitap_t) (benchmark::State & state) { this->run(state); }
BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, bitap_journal, journal_t, bitap_t)(benchmark::State & state) { this->run(state); }

BENCHMARK_REGISTER_F(benchmark_search, bitap_vector)->Apply(benchmark_args);
BENCHMARK_REGISTER_F(benchmark_search, bitap_journal)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark seach shift-or/bitap
// ----------------------------------------------------------------------------

using edit_t = decltype(libjst::myers_pattern_searcher{std::declval<source_t &>()});

BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, edit_vector, source_t, edit_t) (benchmark::State & state) { this->run(state); }
BENCHMARK_TEMPLATE_DEFINE_F(benchmark_search, edit_journal, journal_t, edit_t)(benchmark::State & state) { this->run(state); }

BENCHMARK_REGISTER_F(benchmark_search, edit_vector)->Apply(benchmark_args);
BENCHMARK_REGISTER_F(benchmark_search, edit_journal)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
