// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/utility/container/aligned_allocator.hpp>
// #include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <libjst/journal.hpp>

#include "sequence_variant_simulation.hpp"

/* Benchmark:
 *
 * Data: sequence sizes small to mid 50 100 200 400 800 1600
 *                      large 10KiBi 100KiBi 1MiBi 10MiBi 100MiBi 1GiBi
 *
 * Variation: 0.5% 1% 2% 4% 8% 16% 32% 64%
 *
 * Access:
 *  * read left to right
 *  * random access
 *
 * Modify:
 *  * erase left to right
 *  * erase random access
 *  * insert left to right
 *  * insert random access
 *
 * Types:
 *  * std::string
 *  * journaled_sequence
 */

static void benchmark_args(benchmark::internal::Benchmark* b) {

    for (int i = 7; i <= 20; i+=1) {
        b->Args({1u << i, 0}); // 0 variance
        for (float j = 0.01; j <= 0.1; j += 0.01) {
            int64_t size = 1u << i;
            int64_t errors = std::ceil(size * j);
            b->Args({size, errors});
        }
    }

    // b->Args({64, 4});

    // for (uint32_t i = 10; i <= 1000; i*=10) {
    //     b->Args({i << 10, 0});
    //     for (int j = 5; j <= 160; j *= 2)
    //         b->Args({i << 10, j});
    // }
}

// ----------------------------------------------------------------------------
// Benchmark access forward
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_sequential_access(benchmark::State & state) {

    size_t const sequence_size = state.range(0);
    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));
    auto modified_seq = generate_sequence<container_t>(base_sequence, sequence_variants);

    auto as_sequence = [] (auto && source) {
        if constexpr (std::same_as<container_t, std::vector<char>>) {
            return std::forward<decltype(source)>(source);
        } else {
            return source.sequence();
        }
    };

    auto target_seq = as_sequence(modified_seq);

    size_t b_count{};
    size_t A_count{};
    for (auto _ : state) {
        b_count = 0;
        A_count = 0;
        std::ranges::for_each(target_seq, [&](char v) {
            benchmark::DoNotOptimize(b_count += (v == 'b'));
            benchmark::DoNotOptimize(A_count += (v == 'A'));
        });
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(target_seq.size());
    state.counters["A_count"] = A_count;
    state.counters["b_count"] = b_count;
}

BENCHMARK_TEMPLATE(benchmark_sequential_access, std::vector<char>)->Apply(benchmark_args);
BENCHMARK_TEMPLATE(benchmark_sequential_access, libjst::journal<uint32_t, std::vector<char> &>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark access random
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_random_access(benchmark::State & state) {

    size_t const sequence_size = state.range(0);
    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));
    auto generated_sequence = generate_sequence<container_t>(base_sequence, sequence_variants);

    auto as_sequence = [] (auto && source) {
        if constexpr (std::same_as<container_t, std::vector<char>>) {
            return std::forward<decltype(source)>(source);
        } else {
            return source.sequence();
        }
    };

    auto target_seq = as_sequence(generated_sequence);

    std::vector<size_t> positions{};
    positions.resize(10000);
    std::uniform_int_distribution<> pos_dist{0, static_cast<int>(target_seq.size()) - 1};
    std::mt19937 gen{};
    gen.seed(static_cast<int>(state.range(1)));
    std::ranges::generate(positions, [&](){ return pos_dist(gen); });

    size_t A_count{};
    size_t b_count{};
    for (auto _ : state) {
        b_count = 0;
        A_count = 0;
        std::ranges::for_each(positions, [&] (size_t const index) {
            auto it = std::ranges::next(std::ranges::begin(target_seq), index);
            benchmark::DoNotOptimize(b_count += (*it == 'b'));
            benchmark::DoNotOptimize(A_count += (*it == 'A'));
        });
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(positions.size());
    state.counters["A_count"] = A_count;
    state.counters["b_count"] = b_count;
}

BENCHMARK_TEMPLATE(benchmark_random_access, std::vector<char>)->Apply(benchmark_args);
BENCHMARK_TEMPLATE(benchmark_random_access, libjst::journal<uint32_t, std::vector<char> &>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark record back
// ----------------------------------------------------------------------------

template <typename container_t, typename sequence_t, typename variants_t>
static size_t run_record(benchmark::State & state, sequence_t && base_sequence, variants_t && variants) {

    size_t target_size{};
    for (auto _ : state) {
        container_t target_seq{base_sequence};
        // benchmark::DoNotOptimize(target_seq = container_t{base_sequence});
        std::ptrdiff_t offset{};

        for (auto && variant : variants) {
            record_variant(target_seq, offset, variant);
        }
        if constexpr (std::same_as<container_t, std::vector<char>>) {
            target_size = target_seq.size();
        }  else {
            target_size = target_seq.sequence().size();
        }
    }
    return target_size;
}

template <typename container_t>
void benchmark_sequential_record(benchmark::State & state) {
    size_t const sequence_size = state.range(0);
    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));

    std::stable_sort(sequence_variants.begin(), sequence_variants.end(), [] (auto const & var1, auto const & var2) {
        return (get<0>(var1) + get<1>(var1) <= get<0>(var2));
    });

    size_t target_size = run_record<container_t>(state, base_sequence, std::move(sequence_variants));

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(target_size);
    state.counters["size"] = target_size;
}

BENCHMARK_TEMPLATE(benchmark_sequential_record, std::vector<char>)->Apply(benchmark_args);
BENCHMARK_TEMPLATE(benchmark_sequential_record, libjst::journal<uint32_t, std::vector<char> &>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark record random
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_random_record(benchmark::State & state) {
    size_t const sequence_size = state.range(0);

    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));

    size_t target_size = run_record<container_t>(state, base_sequence, std::move(sequence_variants));

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(target_size);
    state.counters["size"] = target_size;
}

BENCHMARK_TEMPLATE(benchmark_random_record, std::vector<char>)->Apply(benchmark_args);
BENCHMARK_TEMPLATE(benchmark_random_record, libjst::journal<uint32_t, std::vector<char> &>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
