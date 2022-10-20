// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides abstract base class for sequence strategy used by the node label.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/range/views/slice.hpp>

#include <libjst/traversal/sequence_strategy_abstract.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename journal_t>
    class sequence_strategy_prefix : sequence_strategy_abstract<journal_t>
    {
    private:
        using base_t = sequence_strategy_abstract<journal_t>;
    public:

        sequence_strategy_prefix() = default;
        template <typename source_sequence_t>
            requires (!std::same_as<std::remove_cvref_t<source_sequence_t>, sequence_strategy_prefix> &&
                      std::constructible_from<base_t, source_sequence_t>)
        sequence_strategy_prefix(source_sequence_t && source) : base_t{(source_sequence_t &&) source}
        {}

        auto sequence() const noexcept
        {
            return base_t::journal().sequence() | seqan3::views::slice(0, base_t::end_position());
        }
    };
}  // namespace libjst
