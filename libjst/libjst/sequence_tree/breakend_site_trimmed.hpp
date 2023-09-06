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

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename wrapped_breakend_site_t>
        requires (std::is_object_v<wrapped_breakend_site_t>)
    class breakend_site_trimmed
    {
    public:

        using base_site_t = std::remove_cvref_t<wrapped_breakend_site_t>;

        using delta_reference = typename base_site_t::delta_reference;
        using delta_value = typename base_site_t::delta_value;
        using index_type = typename base_site_t::index_type;
        using value_type = typename base_site_t::value_type;
        using position_value_type = libjst::variant_position_t<delta_reference>;

    private:

        wrapped_breakend_site_t _wrappee{};
        position_value_type _max_position{};

    public:

        breakend_site_trimmed() = default;

        explicit constexpr breakend_site_trimmed(wrapped_breakend_site_t breakend_site,
                                                 position_value_type max_position = std::numeric_limits<position_value_type>::max()) :
            _wrappee{std::move(breakend_site)},
            _max_position{max_position}
        {
        }

        constexpr delta_reference operator*() const noexcept {
            return *_wrappee;
        }

        constexpr auto get_breakend() const noexcept {
            return _wrappee.get_breakend();
        }

        constexpr auto get_breakend_site() const noexcept {
            return _wrappee.get_breakend_site();
        }

        constexpr bool is_high_end() const noexcept {
            return _wrappee.is_high_end();
        }

        constexpr bool is_low_end() const noexcept {
            return _wrappee.is_low_end();
        }

        constexpr wrapped_breakend_site_t const & base() const noexcept {
            return _wrappee;
        }

    private:

        friend constexpr position_value_type
        tag_invoke(std::tag_t<libjst::position>, breakend_site_trimmed const & me) noexcept
        {
            return std::min(libjst::position(me._wrappee), me._max_position);
        }

        constexpr friend bool operator==(breakend_site_trimmed const &, breakend_site_trimmed const &) noexcept = default;
    };

}  // namespace libjst
