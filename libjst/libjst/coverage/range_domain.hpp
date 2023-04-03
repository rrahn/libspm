// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a range domain for coverage objects.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <stdint.h>

#include <seqan3/core/concept/cereal.hpp>

namespace libjst
{

    template <std::unsigned_integral value_t = uint32_t>
    class range_domain {

        value_t _min{0};
        value_t _max{std::numeric_limits<value_t>::max()};
    public:

        using value_type = value_t;

        constexpr range_domain() noexcept = default;
        constexpr range_domain(value_type const min, value_type const max) noexcept :
            _min{min},
            _max{std::max(min, max)}
        {}

        constexpr size_t size() const noexcept {
            return _max - _min;
        }

        constexpr bool is_member(value_t elem) const noexcept {
            return elem == std::clamp(elem, _min, _max);
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_min, _max);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_min, _max);
        }

    private:

        friend constexpr bool operator==(range_domain const &, range_domain const &) noexcept = default;
    };
}  // namespace libjst
