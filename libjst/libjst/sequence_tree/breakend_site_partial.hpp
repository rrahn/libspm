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

#include <type_traits>

#include <libjst/sequence_tree/breakend_site.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename breakend_t>
        requires (std::is_object_v<breakend_t>)
    class breakend_site_partial
    {
    public:

        using delta_reference = std::iter_reference_t<breakend_t>;
        using delta_value = std::iter_value_t<breakend_t>;
        using index_type = uint32_t;
        using value_type = uint32_t;

    private:

        breakend_t _original{};
        breakend_t _bound{};
        breakpoint_end _site{};

    public:

        breakend_site_partial() = default;


        explicit constexpr breakend_site_partial(breakend_t original,
                                                 breakend_t bound,
                                                 breakpoint_end site) :
            _original{std::move(original)},
            _bound{std::move(bound)},
            _site{site}
        {
        }

        explicit constexpr breakend_site_partial(breakend_site<breakend_t> other) :
            breakend_site_partial{other.get_breakend(), other.get_breakend(), other.get_breakend_site()}
        {
        }

        constexpr delta_reference operator*() const noexcept {
            return *_bound;
        }

        constexpr breakend_t get_breakend() const noexcept {
            return _original;
        }

        constexpr breakend_t get_bound() const noexcept {
            return _bound;
        }

        constexpr breakpoint_end get_breakend_site() const noexcept {
            return _site;
        }

        constexpr bool is_high_end() const noexcept {
            return get_breakend_site() == breakpoint_end::high;
        }

        constexpr bool is_low_end() const noexcept {
            return get_breakend_site() == breakpoint_end::low;
        }

    private:

        friend constexpr libjst::variant_position_t<delta_reference>
        tag_invoke(std::tag_t<libjst::position>, breakend_site_partial const & me) noexcept
        {
            if (me.is_low_end()) {
                return libjst::low_breakend(*(me.get_breakend()));
            } else {
                return libjst::high_breakend(*(me.get_breakend()));
            }
        }

        constexpr friend bool operator==(breakend_site_partial const & lhs, breakend_site_partial const & rhs) noexcept {
            return (lhs._original == rhs._original) && (lhs._bound == rhs._bound) && (lhs._site == rhs._site);
        }
    };

}  // namespace libjst
