// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides polymorphic sequence searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/matcher/concept.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/state_oblivious_traverser.hpp>

namespace libjst
{

    template <typename polymorphic_sequence_t>
    class polymorphic_sequence_searcher
    {
    private:

        using sequence_type = typename polymorphic_sequence_t::source_type;

        std::reference_wrapper<polymorphic_sequence_t const> _polymorphic_sequence;

    public:
        polymorphic_sequence_searcher() = delete;
        constexpr explicit polymorphic_sequence_searcher(polymorphic_sequence_t const & polymorphic_sequence) :
            _polymorphic_sequence{std::cref(polymorphic_sequence)}
        {}

        template <libjst::window_matcher pattern_t, typename callback_t>
            // requires libjst::online_matcher_for<pattern_t, sequence_type const &, callback_t>
        constexpr void operator()(pattern_t && pattern, callback_t && callback) const {

            auto tree = libjst::volatile_tree{_polymorphic_sequence.get()};
            auto traverser = make_traverser<pattern_t>(_polymorphic_sequence.get());

            traverser(tree, (pattern_t &&) pattern, (callback_t &&)callback);
        }

    private:

        template <typename pattern_t>
        static constexpr auto make_traverser([[maybe_unused]] polymorphic_sequence_t const & polymorphic_sequence) noexcept {
            if constexpr (libjst::restorable_matcher<pattern_t>) {
                if constexpr (libjst::reducable_state<libjst::matcher_state_t<pattern_t>>) {
                    return state_oblivious_traverser{};
                } else {
                    return state_oblivious_traverser{};
                }
            } else {
                return state_oblivious_traverser{};
            }
        }
    };

}  // namespace libjst
