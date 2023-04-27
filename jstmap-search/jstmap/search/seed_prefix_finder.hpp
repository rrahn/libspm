// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/sequence_tree/seek_position.hpp>

namespace jstmap
{
    template <typename finder_t>
    class seed_prefix_finder {
    private:
        using base_t = finder_t;
        finder_t const & _finder;
        std::ptrdiff_t _source_size{};

    public:

        explicit constexpr seed_prefix_finder() = delete;

        explicit constexpr seed_prefix_finder(finder_t const & finder, std::ptrdiff_t source_size) noexcept :
            _finder{finder},
            _source_size{source_size}
        {
        }

    private:

        finder_t const & base() const noexcept {
            return _finder;
        }

        friend auto beginPosition(seed_prefix_finder const & me) noexcept {
            return me._source_size - endPosition(me.base());
        }

        friend auto endPosition(seed_prefix_finder const & me) noexcept {
            return me._source_size - beginPosition(me.base());
        }
    };
}  // namespace jstmap
