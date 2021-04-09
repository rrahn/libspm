// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstring>
#include <vector>
#include <memory_resource>

#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <sdsl/bit_vectors.hpp>

#include <libjst/utility/bit_vector_adaptor.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/utility/bit_vector_simd.hpp>

template <typename result_vector_t>
auto generate_bit_vector_pair(size_t const size)
{
    auto random_bit_sequence_first = seqan3::test::generate_numeric_sequence(size, 0u, 1u);
    auto random_bit_sequence_second = seqan3::test::generate_numeric_sequence(size, 0u, 1u, size);

    result_vector_t test_vector_first(size);
    result_vector_t test_vector_second(size);

    // std::cout << "After generation: " << test_vector_first.capacity() << "\n";
    if constexpr (std::ranges::output_range<result_vector_t, bool>)
    {
        std::ranges::copy(random_bit_sequence_first, test_vector_first.begin());
        std::ranges::copy(random_bit_sequence_second, test_vector_second.begin());
    }
    else
    {
        std::ranges::copy(random_bit_sequence_first, std::back_inserter(test_vector_first));
        std::ranges::copy(random_bit_sequence_second, std::back_inserter(test_vector_second));
    }

    return std::pair{std::move(test_vector_first), std::move(test_vector_second)};
}

template <typename bit_vector_t, typename operation_t>
void benchmark_bit_vector(benchmark::State & state, bit_vector_t, operation_t && operation)
{
    auto [test_vector_lhs, test_vector_rhs] = generate_bit_vector_pair<bit_vector_t>(state.range(0));

    // std::cout << "After copying: " << test_vector_lhs.capacity() << "\n";

    for (auto _ : state)
        benchmark::DoNotOptimize(operation(test_vector_lhs, test_vector_rhs));
}

// ----------------------------------------------------------------------------
// Benchmark operations
// ----------------------------------------------------------------------------

inline auto binary_and = [] <typename bv_t> (bv_t & lhs, bv_t & rhs) constexpr -> bv_t const &
{
    return lhs &= rhs;
};

inline auto binary_or = [] <typename bv_t> (bv_t & lhs, bv_t & rhs) constexpr -> bv_t const &
{
    return lhs |= rhs;
};

inline auto binary_xor = [] <typename bv_t> (bv_t & lhs, bv_t & rhs) constexpr -> bv_t const &
{
    return lhs ^= rhs;
};

inline auto binary_not = [] <typename bv_t> (bv_t & lhs, [[maybe_unused]] auto const & ...args) constexpr -> bv_t
{
    return ~lhs;
};

inline auto binary_flip = [] <typename bv_t> (bv_t & lhs, [[maybe_unused]] auto const & ...args) constexpr -> bv_t const &
{
    return lhs.flip();
};

inline auto none = [] (auto & lhs, [[maybe_unused]] auto const & ...args) constexpr -> bool
{
    return lhs.none();
};

inline auto all = [] (auto & lhs, [[maybe_unused]] auto const & ...args) constexpr -> bool
{
    return lhs.all();
};

inline auto any = [] (auto & lhs, [[maybe_unused]] auto const & ...args) constexpr -> bool
{
    return lhs.any();
};

static constexpr int32_t min_range = 32;
static constexpr int32_t max_range = 1'048'576;
static constexpr int32_t range_multiplier = 2;

// ----------------------------------------------------------------------------
// Benchmark libjst::bit_vector
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_and,
                  libjst::bit_vector<>{},
                  binary_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_or,
                  libjst::bit_vector<>{},
                  binary_or)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_xor,
                  libjst::bit_vector<>{},
                  binary_xor)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_not,
                  libjst::bit_vector<>{},
                  binary_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_flip,
                  libjst::bit_vector<>{},
                  binary_flip)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_none,
                  libjst::bit_vector<>{},
                  none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_all,
                  libjst::bit_vector<>{},
                  all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_any,
                  libjst::bit_vector<>{},
                  any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark libjst::bit_vector_simd
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_and,
                  libjst::bit_vector_simd<>{},
                  binary_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_or,
                  libjst::bit_vector_simd<>{},
                  binary_or)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_xor,
                  libjst::bit_vector_simd<>{},
                  binary_xor)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_not,
                  libjst::bit_vector_simd<>{},
                  binary_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_flip,
                  libjst::bit_vector_simd<>{},
                  binary_flip)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_none,
                  libjst::bit_vector_simd<>{},
                  none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_all,
                  libjst::bit_vector_simd<>{},
                  all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_simd_any,
                  libjst::bit_vector_simd<>{},
                  any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark std::vector<bool> adaptor
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_and,
                  libjst::utility::bit_vector_adaptor{},
                  binary_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_or,
                  libjst::utility::bit_vector_adaptor{},
                  binary_or)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_xor,
                  libjst::utility::bit_vector_adaptor{},
                  binary_xor)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_not,
                  libjst::utility::bit_vector_adaptor{},
                  binary_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_none,
                  libjst::utility::bit_vector_adaptor{},
                  none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_all,
                  libjst::utility::bit_vector_adaptor{},
                  all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  vector_bool_any,
                  libjst::utility::bit_vector_adaptor{},
                  any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark sdsl::bit_vector
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  sdsl_bv_and,
                  sdsl::bit_vector{},
                  binary_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  sdsl_bv_or,
                  sdsl::bit_vector{},
                  binary_or)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  sdsl_bv_xor,
                  sdsl::bit_vector{},
                  binary_xor)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
