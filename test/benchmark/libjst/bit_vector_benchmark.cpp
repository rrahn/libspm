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

#include <seqan3/utility/container/aligned_allocator.hpp>
// #include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <sdsl/bit_vectors.hpp>

#include <libjst/utility/bit_vector_adaptor.hpp>
#include <libjst/utility/bit_vector.hpp>

template <typename result_vector_t>
auto generate_bit_vector_pair(size_t const size)
{
    auto random_bit_sequence_first = seqan3::test::generate_numeric_sequence(size, 0u, 1u);
    auto random_bit_sequence_second = seqan3::test::generate_numeric_sequence(size, 0u, 1u, size);

    result_vector_t test_vector_first(size);
    result_vector_t test_vector_second(size);

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
    auto [lhs, rhs] = generate_bit_vector_pair<bit_vector_t>(state.range(0));

    bit_vector_t res{lhs};

    for (auto _ : state) {
        operation(res, lhs, rhs);
    }

    state.counters["#bits"] = std::ranges::count_if(res, [] (auto v) { return v; });
}

template <typename bit_vector_t, typename operation_t>
void benchmark_bit_vector_reduce(benchmark::State & state, bit_vector_t, operation_t && operation)
{
    auto [lhs, rhs] = generate_bit_vector_pair<bit_vector_t>(state.range(0));

    bool res{false};

    for (auto _ : state) {
        operation(res, lhs);
    }

    state.counters["#bits"] = std::ranges::count_if(lhs, [] (auto v) { return v; });
    state.counters["result"] = res;
}

template <typename bit_vector_t, typename operation_t>
void benchmark_bit_vector_reduce_zero(benchmark::State & state, bit_vector_t, operation_t && operation)
{
    bit_vector_t vec{};
    vec.resize(state.range(0), 0);

    bool res{false};

    for (auto _ : state) {
        operation(res, vec);
    }

    state.counters["#bits"] = std::ranges::count_if(vec, [] (auto v) { return v; });
    state.counters["result"] = res;
}

template <typename bit_vector_t, typename operation_t>
void benchmark_bit_vector_reduce_all(benchmark::State & state, bit_vector_t, operation_t && operation)
{
    bit_vector_t vec{};
    vec.resize(state.range(0), 1);

    bool res{false};

    for (auto _ : state) {
        operation(res, vec);
    }

    state.counters["#bits"] = std::ranges::count_if(vec, [] (auto v) { return v; });
    state.counters["result"] = res;
}

// ----------------------------------------------------------------------------
// Benchmark operations
// ----------------------------------------------------------------------------

inline auto bitparallel_and = [] <typename bv_t> (bv_t & res, bv_t const & lhs, bv_t const & rhs) constexpr
{
    res = lhs & rhs;
};

inline auto bitparallel_and_eq = [] <typename bv_t> (bv_t & lhs, bv_t const &, bv_t const & rhs) constexpr
{
    lhs &= rhs;
};

inline auto bitparallel_and_not = [] <typename bv_t> (bv_t & lhs, bv_t const &, bv_t const & rhs) constexpr
{
    lhs.and_not(rhs);
};

inline auto bitparallel_not = [] <typename bv_t> (bv_t & res, bv_t const & lhs, bv_t const &) constexpr
{
    res = ~lhs;
};

inline auto bitparallel_flip = [] <typename bv_t> (bv_t & res, ...) constexpr
{
    res.flip();
};

inline auto bitparallel_none = [] (auto & res, auto const & bv) constexpr
{
    res = bv.none();
};

inline auto bitparallel_all = [] (auto & res, auto const & bv) constexpr
{
    res = bv.all();
};

inline auto bitparallel_any = [] (auto & res, auto const & bv) constexpr
{
    res = bv.any();
};

static constexpr int32_t min_range = 1 << 5; // 2^5 = 32
static constexpr int32_t max_range = 1 << 22; // 1^22 = 4194304
static constexpr int32_t range_multiplier = 2;

// ----------------------------------------------------------------------------
// Benchmark libjst::bit_vector
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_and,
                  libjst::bit_vector<>{},
                  bitparallel_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_and_eq,
                  libjst::bit_vector<>{},
                  bitparallel_and_eq)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_and_not,
                  libjst::bit_vector<>{},
                  bitparallel_and_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_not,
                  libjst::bit_vector<>{},
                  bitparallel_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector,
                  libjst_bv_flip,
                  libjst::bit_vector<>{},
                  bitparallel_flip)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce,
                  libjst_bv_none,
                  libjst::bit_vector<>{},
                  bitparallel_none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce_zero,
                  libjst_bv_none,
                  libjst::bit_vector<>{},
                  bitparallel_none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce,
                  libjst_bv_all,
                  libjst::bit_vector<>{},
                  bitparallel_all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce_all,
                  libjst_bv_all,
                  libjst::bit_vector<>{},
                  bitparallel_all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce,
                  libjst_bv_any,
                  libjst::bit_vector<>{},
                  bitparallel_any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(benchmark_bit_vector_reduce_zero,
                  libjst_bv_any,
                  libjst::bit_vector<>{},
                  bitparallel_any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// // ----------------------------------------------------------------------------
// // Benchmark std::vector<bool> adaptor
// // ----------------------------------------------------------------------------

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_and,
//                   libjst::utility::bit_vector_adaptor{},
//                   bitparallel_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_and_eq,
//                   libjst::utility::bit_vector_adaptor{},
//                   bitparallel_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);


// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_not,
//                   libjst::utility::bit_vector_adaptor{},
//                   binary_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_none,
//                   libjst::utility::bit_vector_adaptor{},
//                   none)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_all,
//                   libjst::utility::bit_vector_adaptor{},
//                   all)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   vector_bool_any,
//                   libjst::utility::bit_vector_adaptor{},
//                   any)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// // ----------------------------------------------------------------------------
// // Benchmark sdsl::bit_vector
// // ----------------------------------------------------------------------------

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   sdsl_bv_and,
//                   sdsl::bit_vector{},
//                   bitparallel_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// BENCHMARK_CAPTURE(benchmark_bit_vector,
//                   sdsl_bv_and,
//                   sdsl::bit_vector{},
//                   bitparallel_and_eq)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
