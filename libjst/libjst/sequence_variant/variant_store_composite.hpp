// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides composite structure to combine multiple distinct variant stores.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <functional>
#include <optional>
#include <ranges>
#include <type_traits>
#include <tuple>
#include <utility>
#include <variant>

#include <libjst/sequence_variant/concept.hpp>
#include <libjst/sequence_variant/variant_any.hpp>
#include <libjst/sequence_variant/variant_store_iterator.hpp>

namespace libjst
{
    namespace detail
    {
        template <typename ...ts>
        auto make_variant_impl(ts &...) -> std::variant<std::monostate, std::reference_wrapper<ts>...>;

        template <typename ...ts>
        auto make_variant_impl(ts const &...) -> std::variant<std::monostate, std::reference_wrapper<ts const>...>;

        template <typename ...ts>
        auto make_variant_impl(ts &&...) -> std::variant<std::monostate, ts...>;

        template <typename ...ts>
        auto make_variant_impl(ts const &&...) -> std::variant<std::monostate, ts...>;

        template <typename ...ts>
        using make_variant_t = decltype(make_variant_impl(std::declval<ts>()...));
    } // namespace detail

    template <sequence_variant ...variants_t>
        requires (sizeof...(variants_t) > 0)
    class composite_proxy
    {
    private:
        template <typename u>
        using pure_t = std::decay_t<std::unwrap_reference_t<u>>;

        struct monostate_fail {

            template <typename cpo_t>
                requires (std::tag_invocable<cpo_t, variants_t> && ...)
            constexpr friend auto tag_invoke(cpo_t, monostate_fail const &)
                -> std::common_reference_t<std::tag_invoke_result_t<cpo_t, variants_t>...>
            {
                throw std::bad_variant_access{};
            }
        };

        using variant_t = detail::make_variant_t<variants_t...>;
        variant_t _value{};
    public:

        composite_proxy() = default;
        composite_proxy(composite_proxy const & other) = default;
        composite_proxy(composite_proxy &&) = default;
        composite_proxy &operator=(composite_proxy const &) = default;
        composite_proxy &operator=(composite_proxy &&) = default;

        template <typename ...qs>
            requires (!std::same_as<composite_proxy<qs...>, composite_proxy> &&
                     (std::assignable_from<variant_t &, qs const &> || ...))
        composite_proxy(composite_proxy<qs...> const & other)
        {
            std::visit([&] (auto const & arg) { _value = arg; }, other._value);
        }

        template <typename ...qs>
            requires (!std::same_as<composite_proxy<qs...>, composite_proxy> &&
                     (std::assignable_from<variant_t &, qs &> || ...))
        composite_proxy(composite_proxy<qs...> & other)
        {
            std::visit([&] (auto & arg) { _value = arg; }, other._value);
        }

        template <typename ...qs>
            requires (!std::same_as<composite_proxy<qs...>, composite_proxy> &&
                     (std::assignable_from<variant_t &, qs> || ...))
        composite_proxy(composite_proxy<qs...> const && other)
        {
            std::visit([&] <typename arg_t>(arg_t && arg) { _value = (arg_t &&)arg; }, std::move(other._value));
        }

        // Construction from variant type
        template <typename arg_t>
            requires (!std::same_as<std::remove_cvref_t<arg_t>, composite_proxy> &&
                      std::constructible_from<variant_t, arg_t>)
        composite_proxy(arg_t && arg) : _value{(arg_t &&)arg}
        {}

        variant_t const & get() const noexcept
        {
            return _value;
        }

    private:
        // template <typename arg_t>
        // constexpr static decltype(auto) wrap(arg_t && arg) noexcept
        // {
        //     if constexpr (std::same_as<arg_t, std::remove_cvref_t<arg_t> &>)
        //         return std::ref(arg);
        //     else if constexpr (std::same_as<arg_t, std::remove_cvref_t<arg_t> const &>)
        //         return std::cref(arg);
        //     else
        //         return (arg_t &&) arg;
        // }

        template <typename arg_t>
        constexpr static decltype(auto) unwrap(arg_t && arg) noexcept
        {
            using u = std::unwrap_reference_t<std::remove_cvref_t<arg_t>>;
            if constexpr (std::same_as<u, std::monostate>)
                return monostate_fail{};
            else if constexpr (std::same_as<u, std::remove_cvref_t<arg_t>>)
                return (arg_t &&)arg;
            else
                return arg.get();
        }

        template <typename cpo_t>
            requires (std::tag_invocable<cpo_t, variants_t> && ...)
        constexpr friend auto tag_invoke(cpo_t cpo, composite_proxy const &me)
            -> std::common_reference_t<std::tag_invoke_result_t<cpo_t, variants_t>...>
        {
            return std::visit([&] <typename arg_t> (arg_t && arg)
                -> std::tag_invoke_result_t<cpo_t, decltype(unwrap(std::declval<arg_t &&>()))> {
                return std::tag_invoke(cpo, unwrap((arg_t &&)arg));
            }, me._value);
        }
    };

    template <std::ranges::random_access_range ...variant_stores_t>
        requires ((sizeof...(variant_stores_t) > 0) &&
                  (sequence_variant<std::ranges::range_value_t<variant_stores_t>> && ...))
    class variant_store_composite
    {
    protected:

        using value_type = composite_proxy<std::ranges::range_value_t<variant_stores_t>...>;
        using reference = composite_proxy<std::ranges::range_reference_t<variant_stores_t const &>...>;
        using iterator = variant_store_iterator<variant_store_composite>;

        friend iterator;

    private:
        using composite_t = std::tuple<variant_stores_t...>;

        composite_t _variants{};

    public:

        variant_store_composite() = default;

        // ----------------------------------------------------------------------------
        // Access
        // ----------------------------------------------------------------------------

        constexpr reference operator[](size_t offset) const noexcept {
            assert(offset < size());
            bool is_assigned{false};
            std::optional<reference> any_ref_opt{};
            for_each_store(*this, [&] (auto const & store)
            {
                if (is_assigned) return;
                if (std::ranges::size(store) <= offset) {
                    offset -= std::ranges::size(store); // update offset
                } else {
                    is_assigned = true; // signal early termination of subsequent iterations.
                    any_ref_opt = reference{store[offset]};
                }
            });
            assert(any_ref_opt.has_value());
            return *any_ref_opt;
        }

        constexpr size_t size() const {
            size_t _size{};
            for_each_store(*this, [&] (auto const & store) { _size += store.size(); });
            return _size;
        }

        // ----------------------------------------------------------------------------
        // Modifier
        // ----------------------------------------------------------------------------

        iterator insert(value_type any_variant)
        {
            return std::visit([&] <typename arg_t> (arg_t && concrete_variant) {
                return insert((arg_t &&)concrete_variant);
            }, any_variant.get());
        }

        template <sequence_variant variant_t>
            requires (std::constructible_from<std::ranges::range_value_t<variant_stores_t>, variant_t> || ...)
        iterator insert(variant_t && variant)
        {
            bool assigned{false};
            size_t offset{};
            iterator inserted = end();
            for_each_store(*this, [&] <typename store_t> (store_t & store) {
                if (!assigned) {
                    offset += store.size();
                    if constexpr (std::constructible_from<std::ranges::range_value_t<store_t>, variant_t>) {
                        store.reserve(store.size() + 1); // enough memory for push_back.
                        store.push_back(variant);
                        inserted = std::ranges::next(begin(), offset);
                        assigned = true;
                    }
                }
            });
            return inserted;
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        iterator begin() const noexcept { return iterator{*this, 0u}; }
        iterator end() const noexcept { return iterator{*this, size()}; }

    private:
        template <typename this_t, typename fn_t>
        constexpr static void for_each_store(this_t && me, fn_t && fn) noexcept
        {
            [] <typename _this_t, typename _fn_t, size_t ...idx> (_this_t & me, _fn_t && fn, std::index_sequence<idx...>)
            {
                (fn(get<idx>(me._variants)), ...);
            } (me, (fn_t &&) fn, std::make_index_sequence<std::tuple_size_v<composite_t>>());
        }
    };

}  // namespace libjst

namespace std
{
    template <typename ...ts, typename ...us,
              template<typename> typename qualifier_t, template<typename> typename qualifier_u>
        requires (sizeof...(ts) == sizeof...(us)) && requires {
            typename libjst::composite_proxy<std::common_reference_t<qualifier_t<ts>, qualifier_u<us>>...>;
        }
    struct basic_common_reference<libjst::composite_proxy<ts...>, libjst::composite_proxy<us...>, qualifier_t, qualifier_u>
    {
        using type = libjst::composite_proxy<std::common_reference_t<qualifier_t<ts>, qualifier_u<us>>...>;
    };
} // namespace std
