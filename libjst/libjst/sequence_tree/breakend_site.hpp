// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides breakend site for sequence graph.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename breakend_iterator>
    class breakend_site
    {
    private:

        breakend_iterator _breakend{};
        breakpoint_end _site{};

    public:

        using delta_reference = std::iter_reference_t<breakend_iterator>;
        using delta_value = std::iter_value_t<breakend_iterator>;
        using index_type = uint32_t;
        using value_type = uint32_t;

        breakend_site() = default;
        breakend_site(breakend_iterator breakend, breakpoint_end site) :
            _breakend{std::move(breakend)},
            _site{site}
        {
        }

        constexpr delta_reference operator*() const noexcept {
            return *_breakend;
        }

        constexpr breakend_iterator get_breakend() const noexcept {
            return _breakend;
        }

        constexpr bool is_high_end() const noexcept {
            return _site == breakpoint_end::high;
        }

        constexpr bool is_low_end() const noexcept {
            return _site == breakpoint_end::low;
        }

    private:

        friend constexpr libjst::variant_position_t<delta_reference>
        tag_invoke(std::tag_t<libjst::position>, breakend_site const & me) noexcept
        {
            if (me.is_low_end()) {
                return libjst::low_breakend(*me);
            } else {
                return libjst::high_breakend(*me);
            }
        }

        constexpr friend bool operator==(breakend_site const & lhs, breakend_site const & rhs) noexcept {
            return (lhs._breakend == rhs._breakend) && (lhs._site == rhs._site);
        }
    };

}  // namespace libjst
