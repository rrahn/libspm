// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <atomic>
#include <functional>

#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include "fixture_base_ibf.hpp"

namespace just::bench
{
    template <typename capture_t>
    class fixture_pigeonhole_pattern_ibf : public fixture_base_ibf<capture_t>
    {
    private:
        using base_t = fixture_base_ibf<capture_t>;
    public:

        fixture_pigeonhole_pattern_ibf() = default;
        virtual ~fixture_pigeonhole_pattern_ibf() = default;

        template <typename pattern_creator_t>
        void run(::benchmark::State & state, pattern_creator_t && make_pattern)
        {
            std::atomic<size_t> _window_size{};
            auto tree_closure = [&] (size_t window_size) {
                _window_size.store(window_size, std::memory_order_relaxed);
                return libjst::labelled() | libjst::coloured()
                                          | libjst::trim(window_size - 1)
                                          | libjst::prune()
                                          | libjst::left_extend(window_size - 1)
                                          | libjst::merge();
            };

            base_t::run(state, make_pattern, tree_closure, [] (auto const & tree) {
                return libjst::tree_traverser_base{tree};
            });
            this->processed_bytes = base_t::total_bytes(_window_size.load());
        }
    };
}  // namespace just::bench

