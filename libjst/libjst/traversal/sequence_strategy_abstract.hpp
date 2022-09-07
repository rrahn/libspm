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

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename journal_t>
    class sequence_strategy_abstract
    {
    private:
        journal_t _journal{};
        size_t _begin_position{};
        size_t _end_position{};

    public:

        sequence_strategy_abstract() = default;
        template <typename source_sequence_t>
            requires (!std::same_as<std::remove_cvref_t<source_sequence_t>, sequence_strategy_abstract> &&
                      std::constructible_from<journa_t, source_sequence_t>)
        sequence_strategy_abstract(source_sequence_t && source) :
            _journal{(source_sequence_t &&) source},
            _end_position{std::ranges::size(_journal.sequence())}
        {}

        template <sequence_variant variant_t>
        void record(variant_t && variant, size_t size)
        {
            if (libjst::is_insertion(variant)) {
                _journal.record_insertion(libjst::position(variant), libjst::insertion(variant));
            } else if (libjst::is_deletion(variant)) {
                _journal.record_deletion(libjst::position(variant), libjst::deletion(variant));
            } else {
                assert(libjst::is_replacement(variant));
                _journal.record_substitution(libjst::position(variant), libjst::insertion(variant));
            }

            _begin_position = libjst::position(variant);
            _end_position = _begin_position + size;
        }

    protected:
        // returns full sequence, but interface is protected
        journal_t const & journal() const noexcept
        {
            return _journal;
        }

        size_t begin_position() const noexcept
        {
            return _begin_position;
        }

        size_t end_position() const noexcept
        {
            return _end_position;
        }
    };
}  // namespace libjst
