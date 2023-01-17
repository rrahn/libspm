// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function for loading in a reference sequence for the simulated alignment.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <ranges>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/simulate/read_sampler.hpp>

namespace jstmap
{

    read_sampler::read_sampler(rcs_store_t const & rcs_store) : _rcs_store{rcs_store}
    {}

    libjst::tree_stats
    read_sampler::compute_tree_stats(size_t const size) const {
        return libjst::volatile_tree(_rcs_store) | libjst::labelled<libjst::sequence_label_kind::root_path>()
                                                 | libjst::coloured()
                                                 | libjst::trim(size - 1)
                                                 | libjst::prune()
                                                 | libjst::merge()
                                                 | libjst::stats();
    }

    sampled_read_list_type
    read_sampler::sample_reads(libjst::tree_stats const stats,  size_t const count,  size_t const size) const {
        assert(size > 0);
        assert(count > 0);

        std::mt19937 random_engine{42};
        size_t const window_size = size - 1;
        size_t max_sample_size = stats.symbol_count;
        std::vector<size_t> sample_hints = generate_sample_hints(count, max_sample_size, random_engine);

        sampled_read_list_type sampled_reads{};
        sampled_reads.reserve(sample_hints.size());

        auto sample_tree = libjst::volatile_tree(_rcs_store) | libjst::labelled<libjst::sequence_label_kind::root_path>()
                                                             | libjst::coloured()
                                                             | libjst::trim(window_size)
                                                             | libjst::prune()
                                                             | libjst::left_extend(window_size)
                                                             | libjst::merge()
                                                             | libjst::seek();

        libjst::tree_traverser_base path{sample_tree};
        auto hint_it = sample_hints.begin();
        size_t last_hint{0};
        std::ptrdiff_t remaining = *hint_it;
        auto tree_it = path.begin();
        while (tree_it != path.end()) {
            auto cargo = *tree_it;
            auto label = cargo.sequence();
            remaining -= std::max<std::ptrdiff_t>(std::ranges::ssize(label) - window_size, 0);
            if (remaining <= 0) { // found label to sample from
                std::ptrdiff_t label_offset = (std::ranges::ssize(label) + remaining) - size;
                assert(label_offset >= 0);
                auto first = std::ranges::next(std::ranges::begin(label), label_offset, std::ranges::end(label));

                sampled_reads.emplace_back(read_type{first, std::ranges::next(first, size, std::ranges::end(label))},
                                           cargo.position(),
                                           label_offset);
                assert(std::ranges::size(get<0>(sampled_reads.back())) == size);

                last_hint = *hint_it;
                if (++hint_it == sample_hints.end())
                    break;

                remaining = (*hint_it - last_hint) + (label_offset + size); // TODO: fix me!
            } else {
                ++tree_it;
            }
        }
        return sampled_reads;
    }

    std::vector<size_t> read_sampler::generate_sample_hints(size_t const count,
                                                            size_t const max_sample_size,
                                                            std::mt19937 & random_engine) const {

        std::uniform_int_distribution<size_t> sample_dist{0, max_sample_size};
        std::vector<size_t> tmp{};
        tmp.resize(count);
        std::ranges::generate(tmp, [&] () { return sample_dist(random_engine); });
        std::ranges::sort(tmp);
        return tmp;
    }

    // read_type
    // read_sampler::sample_read(libjst::seek_position const pos,
    //                           std::ptrdiff_t const label_offset,
    //                           size_t const read_size) const {

    //     constexpr auto is_N = [] (auto const & symbol) constexpr {
    //         return symbol == seqan3::from_char(alphabet_t{}, 'N');
    //     };

    //     std::ptrdiff_t end_offset = std::ranges::ssize(label) + remaining;
    //     std::ptrdiff_t start_offset = end_offset - size;
    //     assert(end_offset >= size);
    //     assert(start_offset >= 0);
    //     auto first = std::ranges::next(std::ranges::begin(label), start_offset, std::ranges::end(label));

    //     sampled_reads.emplace_back();
    //     assert(std::ranges::size(sampled_reads.back()) == size);

    //     read_type read{first, std::ranges::next(first, size, std::ranges::end(label))}

    //     return sampled_reads;
    // }
} // namespace jstmap
