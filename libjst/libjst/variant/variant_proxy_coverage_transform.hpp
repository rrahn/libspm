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

#include <type_traits>

#include <libcontrib/copyable_box.hpp>

#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <covered_sequence_variant variant_t, std::invocable<variant_coverage_t<variant_t>> transform_fn_t>
        requires std::is_object_v<transform_fn_t>
    class coverage_transform_variant
    {
    private:
        using coverage_t = variant_coverage_t<variant_t>;

        variant_t _wrappee;
        jst::contrib::copyable_box<transform_fn_t> _transform_fn;
    public:
        template <covered_sequence_variant _variant_t, std::invocable<variant_coverage_t<_variant_t>> _transform_fn_t>
            requires (!std::same_as<std::remove_cvref_t<_varaint>, coverage_transform_variant>)
        coverage_transform_variant(_variant_t && variant, _transform_fn_t && transform_fn) :
            _wrappee{(_variant_t &&)variant},
            _transform_fn{(_transform_fn_t &&)transform_fn}
        {}
    private:

        constexpr friend coverage_t tag_invoke(std::tag_t<libjst::coverage>, coverage_transform_variant const & me) noexcept
        {
            return std::invoke(*me._transform_fn, libjst::coverage(me._wrappee));
        }

        template <typename cpo_t, typename this_t>
            requires (std::invocable<cpo_t, this_t> && std::same_as<std::remove_cvref_t<this_t>, coverage_transform_variant>)
        constexpr friend auto tag_invoke(cpo_t cpo, this_t && me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, this_t>)
            -> std::invoke_result_t<cpo_t, this_t>
        {
            return cpo((this_t &&)me);
        }
    };

    template <covered_sequence_variant variant_t, std::invocable<variant_coverage_t<variant_t>> transform_fn_t>
    coverage_transform_variant(variant_t &&, transform_fn_t &&)
        -> coverage_transform_variant<variant_t, std::remove_cvref_t<transform_fn_t>>;
}  // namespace libjst
