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

#include <libcontrib/type_traits.hpp>

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

        template <typename var_t, typename ...args_t>
            requires (!std::same_as<std::remove_cvref_t<var_t>, composite_proxy> &&
                      std::constructible_from<variant_t, std::in_place_type_t<var_t>, args_t...>)
        composite_proxy(std::in_place_type_t<var_t>, args_t &&...args) :
            _value{std::in_place_type<var_t>, (args_t &&)args...}
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

        template <typename this_t, typename member_t>
        using fwd_t = std::conditional_t<std::is_lvalue_reference_v<member_t>,
                                         member_t,
                                         jst::contrib::member_type_t<this_t, std::remove_cvref_t<member_t>>>;

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

        template <typename cpo_t, typename this_t>
            requires std::same_as<std::remove_cvref_t<this_t>, composite_proxy> &&
                     (std::tag_invocable<cpo_t, fwd_t<this_t, variants_t>> && ...)
        constexpr friend auto tag_invoke(cpo_t cpo, this_t &&me)
            -> std::common_reference_t<std::tag_invoke_result_t<cpo_t, fwd_t<this_t, variants_t>>...>
        {
            return std::visit([&] <typename arg_t> (arg_t && arg)
                -> std::tag_invoke_result_t<cpo_t, decltype(unwrap(std::declval<arg_t>()))>
            {
                return std::tag_invoke(cpo, unwrap((arg_t &&)arg));
            }, (fwd_t<this_t, variant_t> &&)me._value);
        }
    };

    template <std::ranges::random_access_range ...variant_stores_t>
        requires ((sizeof...(variant_stores_t) > 0) &&
                  (sequence_variant<std::ranges::range_value_t<variant_stores_t>> && ...))
    class variant_store_composite
    {
    protected:

        using value_type = composite_proxy<std::ranges::range_value_t<variant_stores_t>...>;
        using reference = composite_proxy<std::ranges::range_reference_t<variant_stores_t &>...>;
        using const_reference = composite_proxy<std::ranges::range_reference_t<variant_stores_t const &>...>;
        using iterator = variant_store_iterator<variant_store_composite>;
        using const_iterator = variant_store_iterator<variant_store_composite const>;

        friend iterator;
        friend const_iterator;

    private:
        using composite_t = std::tuple<variant_stores_t...>;

        composite_t _variants{};

    public:

        variant_store_composite() = default;

        // ----------------------------------------------------------------------------
        // Access
        // ----------------------------------------------------------------------------

        constexpr reference operator[](size_t offset) noexcept {
            assert(offset < size());
            return access_impl<0>(_variants, offset);
        }

        constexpr const_reference operator[](size_t offset) const noexcept {
            assert(offset < size());
            return access_impl<0>(_variants, offset);
        }

        constexpr size_t size() const noexcept {
            size_t _size{};
            for_each_store(*this, [&] (auto const & store) { _size += store.size(); });
            return _size;
        }

        // ----------------------------------------------------------------------------
        // Modifier
        // ----------------------------------------------------------------------------
        // template <std::unsigned_integral ...sizes_t>
        //     requires (sizeof...(sizes_t) == sizeof...(std::tuple_size_v<composite_t>))
        // constexpr void resize(sizes_t const...sizes) {
        //     auto size_tpl = std::make_tuple(sizes...);

        //     auto resize_impl = [&] <size_t idx> () {
        //         if constexpr (idx >= std::tuple_size_v<composite_t>)
        //             return; // break condition
        //         else {
        //             get<idx>(_variants).resize(get<idx>(size_tpl));
        //             resize_impl<idx + 1>();
        //         }
        //     };

        //     resize_impl<0>();
        // }

        iterator insert(value_type any_variant)
        {
            return std::visit([&] <typename arg_t> (arg_t && concrete_variant) {
                return insert((arg_t &&)concrete_variant);
            }, any_variant.get());
        }

        template <sequence_variant variant_t>
            requires (std::constructible_from<std::ranges::range_value_t<variant_stores_t>, variant_t> || ...)
        iterator insert(variant_t variant)
        {
            bool assigned{false};
            size_t offset{};
            iterator inserted = end();
            for_each_store(*this, [&] <typename store_t> (store_t & store) {
                if (!assigned) {
                    offset += store.size();
                    if constexpr (std::constructible_from<std::ranges::range_value_t<store_t>, variant_t>) {
                        store.reserve(store.size() + 1); // enough memory for push_back.
                        store.push_back(std::move(variant));
                        inserted = std::ranges::next(begin(), offset);
                        assigned = true;
                    }
                }
            });
            return inserted;
        }

        template <typename ...args_t>
            // requires (std::constructible_from<std::ranges::range_value_t<variant_stores_t>, args_t...> || ...)
        const_iterator emplace(args_t &&...args)
        {
            return emplace_impl<0>(0, (args_t &&)args...);
            // for_each_store(*this, [&, tpl = std::forward_as_tuple((args_t &&)args...)]
            // <typename store_t> ([[maybe_unused]] store_t & store)
            // {
            //     using local_value_t = std::ranges::range_value_t<store_t>;
            //     if (!assigned) {
            //         if constexpr (std::constructible_from<local_value_t, args_t...>) {
            //             std::apply([&] (auto &&..._args) {
            //                 inserted = insert(local_value_t{(decltype(_args) &&)_args...});
            //             }, std::move(tpl));
            //             assigned = true;
            //         }
            //     }
            // });
            // return inserted;
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        iterator begin() noexcept { return iterator{*this, 0u}; }
        const_iterator begin() const noexcept { return const_iterator{*this, 0u}; }
        iterator end() noexcept { return iterator{*this, size()}; }
        const_iterator end() const noexcept { return const_iterator{*this, size()}; }

        // template <seqan3::cereal_archive archive_t>
        // void serialize(archive_t & ioarchive)
        // {
        //     for_each_store(*this, [&] <typename store_t> (store_t & store) {
        //         ioarchive(store);
        //     });
        // }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            for_each_store(*this, [&] <typename store_t> (store_t & store) {
                iarchive(store);
            });
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            for_each_store(*this, [&] <typename store_t> (store_t & store) {
                oarchive(store);
            });
        }

    private:
        template <typename this_t, typename fn_t>
        constexpr static void for_each_store(this_t && me, fn_t && fn) noexcept
        {
            [] <typename _this_t, typename _fn_t, size_t ...idx> (_this_t & me, _fn_t && fn, std::index_sequence<idx...>)
            {
                (fn(get<idx>(me._variants)), ...); // but then we can not return the correct type
            } (me, (fn_t &&) fn, std::make_index_sequence<std::tuple_size_v<composite_t>>());
        }

        // template <size_t idx>
        // static constexpr bool assert_out_of_range = idx >= std::tuple_size_v<composite_t>;

        template <size_t idx, typename variants_t>
            requires std::same_as<std::remove_cvref_t<variants_t>, composite_t>
        static auto access_impl(variants_t & variants, size_t offset) noexcept ->
            std::conditional_t<std::is_const_v<variants_t>, const_reference, reference>
        {
            if constexpr (idx >= std::tuple_size_v<composite_t>)
                // static_assert(assert_out_of_range<idx>, "Out of range!");
                return get<0>(variants)[0]; // dummy return to make it compile.
            else {
                if (size_t store_size = std::ranges::size(get<idx>(variants)); offset < store_size)
                    return get<idx>(variants)[offset];
                else
                    return access_impl<idx + 1>(variants, offset - store_size);
            }
        }

        template <size_t idx, typename ...args_t>
        const_iterator emplace_impl(size_t offset, args_t &&...args)
        {
            if constexpr (idx < std::tuple_size_v<composite_t>) {
                using store_t = std::tuple_element_t<idx, composite_t>;
                using local_value_t = std::ranges::range_value_t<store_t>;
                store_t & store = get<idx>(_variants);
                offset += store.size();
                if constexpr (std::constructible_from<local_value_t, args_t...>) {
                    store.emplace_back((args_t &&)args...);
                    return const_iterator{*this, offset};
                } else {
                    return emplace_impl<idx + 1>(offset, (args_t &&)args...);
                }
            } else {
                return const_iterator{*this, size()};
            }
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
