// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>
#include <libjst/sequence_variant/variant_store_sorted.hpp>

namespace libjst
{

    template <typename jst_t>
    class journaled_sequence_tree_forward_ : private traversable_jst_base
    {
    private:
        using store_t = variant_store_t<jst_t const &>;
        using sorted_store_t = variant_store_sorted<store_t const>;

        class ordered_variant;

        [[no_unique_address]] jst_t const *_jst{};
        sorted_store_t _event_queue{};

    public:

        journaled_sequence_tree_forward_() = default;
        explicit journaled_sequence_tree_forward_(jst_t const & jst) :
            _jst{std::addressof(jst)},
            _event_queue{libjst::variant_store(jst)}
        {}

    private:
        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load>,
                                         journaled_sequence_tree_forward_ & me,
                                         archive_t & archive)
        {
            libjst::load_extern(archive, *me._jst);
            archive(me._event_queue);
        }

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save>,
                                         journaled_sequence_tree_forward_ const & me,
                                         archive_t & archive)
        {
            libjst::save_extern(archive, *me._jst);
            archive(me._event_queue);
        }

        template <typename cpo_t>
            requires std::invocable<cpo_t, jst_t const &>
        constexpr friend auto tag_invoke(cpo_t cpo, journaled_sequence_tree_forward_ const &me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, jst_t const &>)
            -> std::invoke_result_t<cpo_t, jst_t const &>
        {
            assert(me._jst != nullptr);
            return cpo(*me._jst);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::variant_store>,
                                         journaled_sequence_tree_forward_ const &me) noexcept
            -> sorted_store_t const &
        {
            return me._event_queue;
        }
    };
}  // namespace libjst
