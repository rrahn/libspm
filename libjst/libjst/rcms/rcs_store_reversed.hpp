// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides referentially compressed sequence store.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>

#include <libjst/rcms/compressed_multisequence_reversed.hpp>

namespace libjst
{
    // Client needs to make sure that the value types are compatible!?
    template <typename cms_t>
    class rcs_store_reversed
    {
    private:
        using wrapper_t = compressed_multisequence_reversed<cms_t>;

        wrapper_t _variant_map{};

    public:

        using variant_map_type = wrapper_t;
        using source_type = typename wrapper_t::source_type;
        using value_type = std::ranges::range_value_t<wrapper_t>;
        using reference = std::ranges::range_reference_t<wrapper_t const &>;
        using size_type = size_t;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr rcs_store_reversed() = delete; //!< Default.
        // Construct from dimensions?
        explicit constexpr rcs_store_reversed(cms_t const & wrappee) :
            _variant_map{wrappee}
        {}
        //!\}

        // ----------------------------------------------------------------------------
        // Accessor
        // ----------------------------------------------------------------------------

        constexpr source_type source() const noexcept
        {
            return variants().source();
        }

        constexpr variant_map_type const & variants() const noexcept
        {
            return _variant_map;
        }

        constexpr size_type size() const noexcept
        {
            return variants().coverage_domain().size();
        }
    };
}  // namespace libjst
