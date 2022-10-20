// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides proxy variant setting an additional position offset.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <sequence_variant variant_t>
    class offset_variant
    {
    private:
        using position_t = variant_position_t<variant_t>;

        variant_t  _wrappee;
        position_t _offset;
    public:
        template <typename _variant_t>
        offset_variant(_variant_t && variant, position_t offset) :
            _wrappee{(_variant_t &&)variant},
            _offset{offset}
        {}
    private:

        constexpr friend position_t tag_invoke(std::tag_t<libjst::position>, offset_variant const & me) noexcept
        {
            return libjst::position(me._wrappee) + me._offset;
        }

        template <typename cpo_t, typename this_t>
            requires (std::invocable<cpo_t, this_t> && std::same_as<std::remove_cvref_t<this_t>, offset_variant>)
        constexpr friend auto tag_invoke(cpo_t cpo, this_t && me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, this_t>)
            -> std::invoke_result_t<cpo_t, this_t>
        {
            return cpo((this_t &&)me);
        }
    };

    template <sequence_variant variant_t>
    offset_variant(variant_t &&, variant_position_t<variant_t> const) -> offset_variant<variant_t>;
}  // namespace libjst
