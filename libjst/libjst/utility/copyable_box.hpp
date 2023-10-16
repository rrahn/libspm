// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides copyable box implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <optional>
#include <type_traits>

namespace libjst
{

    namespace detail
    {

        template <typename derived_t, typename type>
        class make_copyable : private std::optional<type>
        {
        private:
            using base_t = std::optional<type>;

        public:

            using base_t::base_t; // inherit all constuctors.

            constexpr derived_t& operator=(derived_t const & other)
                noexcept(std::is_nothrow_copy_assignable_v<type>)
                requires (std::copyable<type>)
            {
                static_cast<base_t &>(*as_derived()) = static_cast<base_t const &>(other);
            }

            constexpr derived_t& operator=(derived_t const & other)
                noexcept(std::is_nothrow_copy_constructible_v<type>)
                requires (!std::copyable<type>)
            {
                if (as_derived() != std::addressof(other)) {
                    if (other) emplace(*other);
                    else base_t::reset();
                }
                return *as_derived();
            }

            constexpr derived_t& operator=(derived_t && other)
                noexcept(std::is_nothrow_move_assignable_v<type>)
                requires (std::movable<type>)
            {
                static_cast<base_t &>(*as_derived()) = static_cast<base_t &&>(other);
            }

            constexpr derived_t& operator=(derived_t && other)
                noexcept(std::is_nothrow_move_constructible_v<type>)
                requires (!std::movable<type>)
            {
                if (as_derived() != std::addressof(other)) {
                    if (other) emplace(std::move(*other));
                    else base_t::reset();
                }
                return *as_derived();
            }

            using base_t::operator bool;
            using base_t::operator*;
            using base_t::operator->;
            using base_t::has_value;
            using base_t::emplace;
            using base_t::reset;
            using base_t::value;

        private:

            derived_t * as_derived() noexcept
            {
                return static_cast<derived_t *>(this);
            }

            derived_t const * as_derived() const noexcept
            {
                return static_cast<derived_t const *>(this);
            }
        };
    } // namespace detail

    // we want to inherit certain operations only if applicable!
    template <typename type>
        requires (std::copy_constructible<type> && std::is_object_v<type>)
    class copyable_box : public detail::make_copyable<copyable_box<type>, type> // we can make a decision based on the features we want to implement.
    {
        using base_t = detail::make_copyable<copyable_box<type>, type>;

    public:
        // conditionally initialisable
        constexpr copyable_box() noexcept(std::is_nothrow_default_constructible_v<type>)
            requires std::default_initializable<type>
            : copyable_box(std::in_place) // default construction of base class.
        {}

        constexpr copyable_box(copyable_box const &)
            noexcept(std::is_nothrow_copy_constructible_v<base_t>) = default;
        constexpr copyable_box(copyable_box &&)
            noexcept(std::is_nothrow_move_constructible_v<base_t>) = default;

        constexpr copyable_box & operator=(copyable_box const &)
            noexcept(std::is_nothrow_copy_assignable_v<base_t>) = default;
        constexpr copyable_box & operator=(copyable_box &&)
            noexcept(std::is_nothrow_move_assignable_v<base_t>) = default;

        template<typename ...args_t>
            requires std::constructible_from<type, args_t...>
        constexpr explicit copyable_box(std::in_place_t, args_t &&...args)
            noexcept(std::is_nothrow_constructible_v<type, args_t...>)
            : base_t{std::in_place, (args_t &&)args...}
        {}

        template <typename _type = type>
            requires (std::constructible_from<type, _type> &&
                      !std::same_as<std::remove_cvref_t<_type>, std::in_place_t> &&
                      !std::same_as<std::remove_cvref_t<_type>, copyable_box>)
        explicit constexpr copyable_box(_type && value)
            noexcept(std::is_nothrow_constructible_v<base_t, _type>)
            : base_t{(_type &&)value}
        {}

        template <typename _type = type>
            requires (!std::same_as<std::remove_cvref_t<_type>, copyable_box>  &&
                      std::is_constructible_v<type, _type> &&
                      std::is_assignable_v<type &, _type>)
        constexpr copyable_box &operator=(_type && value)
            noexcept(std::is_nothrow_assignable_v<base_t &, _type>)
        {
            static_cast<base_t &>(*this) = (_type &&)value;
            return *this;
        }
        // all other operations are inherited?

        template <typename ...args_t>
            requires std::invocable<type &, args_t...>
        constexpr auto operator()(args_t && ...args)
            noexcept(std::is_nothrow_invocable_v<type &, args_t...>)
            -> std::invoke_result_t<type &, args_t...> {
            return std::invoke(base_t::value(), (args_t &&) args...);
        }

        template <typename ...args_t>
            requires std::invocable<type const &, args_t...>
        constexpr auto operator()(args_t && ...args) const
            noexcept(std::is_nothrow_invocable_v<type const &, args_t...>)
            -> std::invoke_result_t<type const &, args_t...> {
            return std::invoke(base_t::value(), (args_t &&) args...);
        }
    };
}  // namespace libjst
