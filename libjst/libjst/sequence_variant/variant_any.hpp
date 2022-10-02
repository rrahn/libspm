// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a type-erased sequence variant.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <memory>
#include <ranges>

#include <libjst/sequence_variant/concept.hpp>

namespace libjst
{
    template <std::integral position_t, std::ranges::forward_range insertion_t, std::integral deletion_t>
    class any_variant
    {
    private:

        template <typename variant_t>
        static constexpr bool is_compatible_v =
            std::constructible_from<position_t, libjst::variant_position_t<variant_t>> &&
            std::constructible_from<insertion_t, libjst::variant_insertion_t<variant_t>> &&
            std::constructible_from<deletion_t, libjst::variant_deletion_t<variant_t>>;
    public:
        constexpr any_variant() = default;

        // TODO: we could even make it more DRY but reusing the move assignment after default construcing this.
        // template <sequence_variant variant_t>
        // constexpr explicit any_variant(variant_t &&var) : any_variant{}//_var{std::make_unique<impl<variant_t>>((variant_t &&)var)}
        // {
        //     *this = (variant_t &&)var;
        // }

        template <sequence_variant variant_t>
            requires is_compatible_v<variant_t>
        constexpr explicit any_variant(variant_t &&var) : _var{std::make_unique<impl<variant_t>>((variant_t &&)var)}
        {
        }

        template <sequence_variant variant_t>
            requires is_compatible_v<variant_t>
        constexpr any_variant &operator=(variant_t &&var)
        {
            if constexpr (std::is_lvalue_reference_v<variant_t>) { // deep copy, achieved by stacking
                _var = std::make_unique<impl<variant_t>>((variant_t &&)var);
            } else { // move assignment!
                _var = std::move(var._var);
            }
            return *this;
        }

    private:
        template <typename cpo_t>
        constexpr friend decltype(auto) tag_invoke(cpo_t cpo, any_variant const &me) noexcept
        {
            return me._var->apply(cpo);
        }

        struct base
        {
            virtual ~base() = default;
            virtual position_t position() const noexcept = 0;
            virtual deletion_t deletion() const noexcept = 0;
            virtual insertion_t insertion() const noexcept = 0;

            template <typename cpo_t>
            decltype(auto) apply(cpo_t const &) const noexcept
            {
                if constexpr (std::same_as<cpo_t, std::tag_t<libjst::position>>)
                    return position();
                else if constexpr (std::same_as<cpo_t, std::tag_t<libjst::deletion>>)
                    return deletion();
                else if constexpr (std::same_as<cpo_t, std::tag_t<libjst::insertion>>)
                    return insertion();
            }
        };

        template <typename var_t>
        struct impl final : public base
        {
            [[no_unique_address]] var_t _var;

            impl(var_t var) : _var{(var_t &&)var}
            {
            }

            position_t position() const noexcept override
            {
                return libjst::position(_var);
            }

            deletion_t deletion() const noexcept override
            {
                return libjst::deletion(_var);
            }

            insertion_t insertion() const noexcept override
            {
                return libjst::insertion(_var);
            }
        };

        std::unique_ptr<base> _var{};
    };
} // namespace libjst
