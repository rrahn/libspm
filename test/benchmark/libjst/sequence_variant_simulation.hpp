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

#pragma once

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

auto generate_variants(size_t const source_size, int const variant_count) {
    using sequence_t = std::vector<char>;
    // number of SNVs
    // number of small InDels
    // InDel Sizes
    // distribution of errors:
        // 99% = SNV
        // 1% = InDels
        //
    size_t max_indel_size = 1; //std::clamp(1, static_cast<int>(std::ceil(0.01*source_size)), 10);

    std::mt19937 gen(42);
    std::uniform_int_distribution<> position_dist(0, source_size - max_indel_size);

    std::vector<bool> free_positions{};
    free_positions.resize(source_size, true);

    size_t snv_count{};
    size_t indel_count{};
    if (variant_count <= 10) { // Only SNVs
        snv_count = variant_count;
    } else { // SNV and InDel
        snv_count = std::floor(variant_count * 0.99);
        indel_count = std::max<int>(variant_count - snv_count, 0);
    }
    // std::cout << "Generate SNVs\n";
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

    // std::cout << "Generate indels\n";
    for (size_t i = 0; i < indel_count; ) {
        bool is_deletion = i & 1;
        size_t pos = position_dist(gen);
        if (free_positions[pos] == true) {
            variants.emplace_back(pos, is_deletion, (is_deletion) ? sequence_t{} : sequence_t(1, 'b'));
            free_positions[pos] = false;
            ++i;
        }
    }
    // std::cout << "Finished generation\n";
    return variants;
}

template <typename container_t, typename variant_t>
void record_variant(container_t & sequence, std::ptrdiff_t & offset, variant_t && variant) {

    auto && [p, l, s] = variant;
    size_t begin_pos = p + offset;

    // std::cout << "var: " << p << " " << l << " " << std::string{s.begin(), s.end()} << "\n";

    switch(kind(variant)) {
        case variant_kind::snv: {
            if constexpr (std::same_as<container_t, std::vector<char>>) {
                benchmark::DoNotOptimize(std::ranges::copy(s, sequence.begin() + begin_pos));
                // sequence.erase(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + l));
                // sequence.insert(sequence.begin() + begin_pos, s.begin(), s.end());
            } else {
                sequence.replace(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + std::ranges::size(s)), s);
            }
            break;
        } case variant_kind::deletion: {
            if constexpr (std::same_as<container_t, std::vector<char>>) {
                sequence.erase(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + l));
            } else {
                sequence.erase(sequence.begin() + begin_pos, sequence.begin() + (begin_pos + l));
            }
            break;
        } case variant_kind::insertion: {
            if constexpr (std::same_as<container_t, std::vector<char>>) {
                sequence.insert(sequence.begin() + begin_pos, s.begin(), s.end());
            } else {
                sequence.insert(sequence.begin() + begin_pos, s);
            }
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
