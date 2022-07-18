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

#include <libjst/utility/sorted_vector.hpp>

static constexpr int32_t min_range = 1ull<<0;
static constexpr int32_t max_range = 1ull<<10;

template <typename container_t>
inline void insert_element(container_t & container, size_t const & value)
{
    if constexpr (std::ranges::random_access_range<container_t>)
    {
        // container.insert(std::ranges::upper_bound(container, value), value);
        container.push_back(value);
    }
    else
        container.insert(value);
}

template <typename container_t>
inline bool contains(container_t & container, size_t const & value)
{
    if constexpr (std::ranges::random_access_range<container_t>)
        return std::ranges::binary_search(container, value);
    else
        return container.contains(value);
}

// ----------------------------------------------------------------------------
// Benchmark insert random
// ----------------------------------------------------------------------------

template <typename dictionary_t, typename require_sort_t>
void benchmark_insert_random(benchmark::State & state, dictionary_t, require_sort_t)
{
    size_t const size = state.range(0);

    std::vector<size_t> elements;
    std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });

    dictionary_t dictionary{};

    for (auto _ : state)
    {
        state.PauseTiming();
        dictionary.clear();
        state.ResumeTiming();
        for (size_t element : elements) {
            dictionary.insert(std::ranges::end(dictionary), element);
        }

        if constexpr(require_sort_t::value)
            std::ranges::sort(dictionary);
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(elements.size());
    state.counters["size"] = dictionary.size();
}

BENCHMARK_CAPTURE(benchmark_insert_random, std_vector, std::vector<size_t>{}, std::true_type{})
    ->Range(min_range, max_range);
BENCHMARK_CAPTURE(benchmark_insert_random, sorted_vector, libjst::sorted_vector<size_t>{}, std::false_type{})
    ->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark insert back
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_insert_back(benchmark::State & state, container_t)
{
    size_t const size = state.range(0);

    std::vector<size_t> elements;
    std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });
    std::ranges::sort(elements);

    container_t cont{};

    for (auto _ : state)
    {
        cont.clear();
        for (size_t element : elements)
            cont.insert(cont.end(), element);

        // benchmark::DoNotOptimize(cont.size());
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(elements.size());
    state.counters["size"] = cont.size();
}

BENCHMARK_CAPTURE(benchmark_insert_back, std_vector, std::vector<size_t>{})
    ->Range(min_range, max_range);
BENCHMARK_CAPTURE(benchmark_insert_back, sorted_vector, libjst::sorted_vector<size_t>{})
    ->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark contains
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_contains(benchmark::State & state, container_t)
{
    size_t const size = state.range(0);

    std::vector<size_t> elements;
    std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });
    std::ranges::sort(elements);

    container_t cont{};
    std::ranges::for_each(elements, [&] (size_t const & element) { cont.insert(cont.end(), element); });

    bool found{};
    for (auto _ : state)
    {
        for (size_t pivot : elements)
            found = contains(cont, pivot);
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(elements.size() * std::log2l(cont.size()));
    state.counters["found"] = found;
}

BENCHMARK_CAPTURE(benchmark_contains, std_vector, std::vector<size_t>{})
    ->Range(min_range, max_range);
BENCHMARK_CAPTURE(benchmark_contains, sorted_vector, libjst::sorted_vector<size_t>{})
    ->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Benchmark read forward
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_access_linear(benchmark::State & state, container_t)
{
    size_t const size = state.range(0);

    std::vector<size_t> elements;
    std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });
    std::ranges::sort(elements);

    container_t cont{};
    std::ranges::for_each(elements, [&] (size_t const & element) { cont.insert(cont.end(), element); });

    size_t const pivot = elements[rand() % elements.size()];
    size_t hits{};
    for (auto _ : state)
    {
        hits = std::ranges::count(cont, pivot);
    }

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(elements.size());
    state.counters["hits"] = hits;
}

BENCHMARK_CAPTURE(benchmark_access_linear, std_vector, std::vector<size_t>{})
    ->Range(min_range, max_range);
BENCHMARK_CAPTURE(benchmark_access_linear, sorted_vector, libjst::sorted_vector<size_t>{})
    ->Range(min_range, max_range);


// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
