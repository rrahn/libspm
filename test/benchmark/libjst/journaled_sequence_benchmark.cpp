// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstring>
#include <memory_resource>
#include <random>
#include <string>
#include <vector>

#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <libjst/journal_decorator_revertable.hpp>
#include <libjst/journal_decorator.hpp>
#include <libjst/utility/position_map.hpp>

template <typename pos_map_t>
class sequence_iterator {

    using iterator_t = std::ranges::iterator_t<pos_map_t>;
    using reference_t = std::iter_reference_t<iterator_t>;
    using sequence_t = typename reference_t::second_type;
    using sequence_iterator_t = std::ranges::iterator_t<sequence_t>;

    iterator_t _dict_it{};
    size_t _position{};
    size_t _next_switch{};
    sequence_iterator_t _sequence_it{};
    // reference_t _active_node{};
    // size_t _end_position{};

public:

    using value_type = char;
    using reference = char const &;
    using difference_type = std::ptrdiff_t;
    using pointer_type = void;
    using iterator_category = std::forward_iterator_tag;

    sequence_iterator() = default;
    sequence_iterator(iterator_t dict_it, size_t position) :
        _dict_it{dict_it},
        _position{position}
    {
        if (position == _next_switch) {
            _next_switch = (*_dict_it).first + (*_dict_it).second.size();
            _sequence_it = (*_dict_it).second.begin();
        }
        // _active_node = *_dict_it;
    }

    reference operator*() const noexcept {
        return *_sequence_it;
    }

    sequence_iterator & operator++() noexcept {
        ++_sequence_it;
        if (++_position == _next_switch) [[unlikely]] {
            ++_dict_it;
            _next_switch = (*_dict_it).first + (*_dict_it).second.size();
            _sequence_it = (*_dict_it).second.begin();
            // _end_position = _active_node.first + _active_node.second.size();
        }
        return *this;
    }

    sequence_iterator operator++(int) noexcept {
        sequence_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr bool operator==(sequence_iterator const & rhs) const noexcept {
        return _position == rhs._position;
    }
};

// template <typename container_t>
// inline void insert_element(container_t & container, size_t const & value)
// {
//     if constexpr (std::ranges::random_access_range<container_t>)
//     {
//         // container.insert(std::ranges::upper_bound(container, value), value);
//         container.push_back(value);
//     }
//     else
//         container.insert(value);
// }

// template <typename container_t>
// inline bool contains(container_t & container, size_t const & value)
// {
//     if constexpr (std::ranges::random_access_range<container_t>)
//         return std::ranges::binary_search(container, value);
//     else
//         return container.contains(value);
// }

// template <typename container_t>
// void benchmark_insert_random(benchmark::State & state, container_t)
// {
//     size_t const size = state.range(0);

//     std::vector<size_t> elements;
//     std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });

//     container_t cont{};

//     for (auto _ : state)
//     {
//         cont.clear();
//         for (size_t element : elements)
//             insert_element(cont, element);

//         if constexpr (std::ranges::random_access_range<container_t>)
//             std::ranges::sort(cont);

//         benchmark::DoNotOptimize(cont.size());
//     }
// }

// template <typename container_t>
// void benchmark_insert_back(benchmark::State & state, container_t)
// {
//     size_t const size = state.range(0);

//     std::vector<size_t> elements;
//     std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });
//     std::ranges::sort(elements);

//     container_t cont{};

//     for (auto _ : state)
//     {
//         cont.clear();
//         for (size_t element : elements)
//             if constexpr (std::ranges::random_access_range<container_t>)
//                 cont.push_back(element);
//             else
//                 cont.insert(cont.end(), element);

//         benchmark::DoNotOptimize(cont.size());
//     }
// }

// template <typename container_t>
// void benchmark_contains(benchmark::State & state, container_t)
// {
//     size_t const size = state.range(0);

//     std::vector<size_t> elements;
//     std::ranges::generate_n(std::back_inserter(elements), size, [] () { return static_cast<size_t>(std::rand()); });
//     std::ranges::sort(elements);

//     container_t cont{};
//     std::ranges::for_each(elements, [&] (size_t const & element) { cont.insert(cont.end(), element); });

//     size_t const pivot = rand();
//     for (auto _ : state)
//     {
//         bool found = contains(cont, pivot);
//         benchmark::DoNotOptimize(found);
//     }
// }
enum struct variant_kind {
    snv,
    insertion,
    deletion
};

template <typename variant_t>
variant_kind kind(variant_t const & var) {
    auto const & [p, l, s] = var;
    if (l == s.size())
        return variant_kind::snv;
    else if (l == 0)
        return variant_kind::insertion;
    else
        return variant_kind::deletion;
}

auto generate_variants(size_t const source_size, int const variation) {
    using sequence_t = std::vector<char>;
    size_t variant_count = variation; //std::ceil(static_cast<float>(source_size) * (variation/100.0));
    // number of SNVs
    // number of small InDels
    // InDel Sizes
    // distribution of errors:
        // 99% = SNV
        // 1% = InDels
        //
    size_t max_indel_size = std::clamp(1, static_cast<int>(std::ceil(0.01*source_size)), 10);

    std::mt19937 gen(42);
    std::uniform_int_distribution<> position_dist(0, source_size - max_indel_size);
    std::uniform_int_distribution<> indel_size_dist(1, max_indel_size);

    std::vector<bool> free_positions{};
    free_positions.resize(source_size, true);

    size_t snv_count{};
    size_t indel_count{};
    if (variant_count <= 10) { // Only SNVs
        snv_count = variant_count;
    } else { // SNV and InDel
        snv_count = variant_count; //std::floor(variant_count * 1.00);
        // indel_count = std::max<int>(variant_count - snv_count, 0);
    }

    std::vector<std::tuple<size_t, size_t, std::vector<char>>> variants{};
    variants.reserve(variant_count);
    for (size_t i = 0; i < snv_count; ) {
        size_t pos = position_dist(gen);
        if (free_positions[pos] == true) {
            variants.emplace_back(pos, 1, sequence_t(1, 'b'));
            free_positions[pos] = false;
            ++i;
        }
    }

    using namespace std::literals;

    for (size_t i = 0; i < indel_count; ) {
        bool is_deletion = i & 1;
        size_t pos = position_dist(gen);
        size_t indel_size = indel_size_dist(gen);
        auto pos_it = free_positions.begin() + pos;
        auto pos_end = pos_it + (is_deletion * indel_size) + !is_deletion;
        if (std::all_of(pos_it, pos_end, [](bool value) { return value; })) {
            variants.emplace_back(pos, indel_size*(i&1), (i&1) ? sequence_t{} : sequence_t(indel_size, 'b'));
            for (; pos_it != pos_end; ++pos_it)
                *pos_it = false;

            ++i;
        }
    }
    return variants;
}

template <typename container_t, typename variant_t>
void record_variant(container_t & sequence, std::ptrdiff_t & offset, variant_t && variant) {

    auto && [p, l, s] = variant;
    size_t begin_pos = p + offset;

    switch(kind(variant)) {
        case variant_kind::snv: {
            if constexpr (std::same_as<container_t, std::vector<char>>) {
                benchmark::DoNotOptimize(std::ranges::copy(s, sequence.begin() + begin_pos));
                // sequence.erase(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + l));
                // sequence.insert(sequence.begin() + begin_pos, s.begin(), s.end());
            } else {
                // sequence.record_substitution(begin_pos, std::span{s});
                sequence.record_substitution(begin_pos, std::span{s});
            }
            break;
        } case variant_kind::deletion: {
            sequence.erase(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + l));
            break;
        } case variant_kind::insertion: {
            std::span s_proxy{s};
            sequence.insert(sequence.begin() + begin_pos, s_proxy.begin(), s_proxy.end());
            break;
        }
    }
    offset += std::ranges::ssize(s) - static_cast<std::ptrdiff_t>(l);
}

template <typename container_t, typename sequence_variants_t>
auto generate_sequence(std::vector<char> & base_sequence,
                       sequence_variants_t & sequence_variants) {

    // sort variants and apply in order!
    std::sort(sequence_variants.begin(), sequence_variants.end(), [] (auto const & var1, auto const & var2) {
        return (get<0>(var1) + get<1>(var1) <= get<0>(var2));
    });

    container_t target_seq{base_sequence};
    std::ptrdiff_t offset = 0;
    std::ranges::for_each(sequence_variants, [&] (auto && variant) {
        record_variant(target_seq, offset, variant);
    });

    return target_seq;
}

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

    for (int i = 6; i <= 10; i+=1) {
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
void benchmark_access_forward(benchmark::State & state) {

    size_t const sequence_size = state.range(0);
    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));
    auto target_seq = generate_sequence<container_t>(base_sequence, sequence_variants);

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

// BENCHMARK_TEMPLATE(benchmark_access_forward, std::vector<char>)->Apply(benchmark_args);
// BENCHMARK_TEMPLATE(benchmark_access_forward, libjst::journal_decorator<std::span<char>>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark access random
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_access_random(benchmark::State & state) {

    size_t const sequence_size = state.range(0);
    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));
    auto target_seq = generate_sequence<container_t>(base_sequence, sequence_variants);

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

// BENCHMARK_TEMPLATE(benchmark_access_random, std::vector<char>)->Apply(benchmark_args);
// BENCHMARK_TEMPLATE(benchmark_access_random, libjst::journal_decorator<std::span<char>>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark record back
// ----------------------------------------------------------------------------

template <typename container_t, typename sequence_t, typename variants_t>
static size_t run_record(benchmark::State & state, sequence_t && base_sequence, variants_t && variants) {

    size_t target_size{};
    for (auto _ : state) {
        state.PauseTiming();
        container_t target_seq{base_sequence};
        std::ptrdiff_t offset{};
        state.ResumeTiming();

        for (auto && variant : variants) {
            record_variant(target_seq, offset, variant);
        }
        target_size = target_seq.size();
    }
    return target_size;
}

template <typename container_t>
void benchmark_record_back(benchmark::State & state) {
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

BENCHMARK_TEMPLATE(benchmark_record_back, std::vector<char>)->Apply(benchmark_args);
BENCHMARK_TEMPLATE(benchmark_record_back, libjst::journal_decorator<std::span<char>>)->Apply(benchmark_args);
// BENCHMARK_TEMPLATE(benchmark_record_back, libjst::journal_decorator_revertable<libjst::journal_decorator<std::span<char>>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Benchmark record random
// ----------------------------------------------------------------------------

template <typename container_t>
void benchmark_record_random(benchmark::State & state) {
    size_t const sequence_size = state.range(0);

    std::vector<char> base_sequence{};
    base_sequence.resize(sequence_size, 'A');

    auto sequence_variants = generate_variants(base_sequence.size(), state.range(1));

    size_t target_size = run_record<container_t>(state, base_sequence, std::move(sequence_variants));

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(target_size);
    state.counters["size"] = target_size;
}


// BENCHMARK_TEMPLATE(benchmark_record_random, std::vector<char>)->Apply(benchmark_args);
// BENCHMARK_TEMPLATE(benchmark_record_random, libjst::journal_decorator<std::span<char>>)->Apply(benchmark_args);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
