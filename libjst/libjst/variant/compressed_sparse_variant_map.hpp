// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a compressed sequence variant matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>

#include <cereal/types/vector.hpp>

#include <libcontrib/type_traits.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/variant/concept.hpp>
#include <libjst/variant/variant_store_iterator.hpp>

namespace libjst
{
    namespace detail {

        enum struct proxy_type_category {
                object = 0,
                reference = 1,
                const_reference = 2
        };

        template <template <typename> typename qualifier_t1,
                  template <typename> typename qualifier_t2>
        struct common_category_selector {
        private:

            template <template <typename> typename qualifier_t>
            static constexpr proxy_type_category get_target_category_for() noexcept {
                if constexpr (std::is_lvalue_reference_v<qualifier_t<int>>) {
                    if constexpr (std::is_const_v<std::remove_reference_t<qualifier_t<int>>>)
                        return proxy_type_category::const_reference;
                    else
                        return proxy_type_category::reference;
                } else {
                    return proxy_type_category::object;
                }
            }

            template <proxy_type_category category1, proxy_type_category category2>
            static constexpr proxy_type_category common_reference_category() noexcept {
                using rank_t = std::underlying_type_t<proxy_type_category>;

                auto to_rank = [] (proxy_type_category category) constexpr noexcept -> rank_t {
                    return static_cast<rank_t>(category);
                };

                constexpr proxy_type_category target_category1 = get_target_category_for<qualifier_t1>();
                constexpr proxy_type_category target_category2 = get_target_category_for<qualifier_t2>();
                constexpr rank_t qualified_category1 = std::max(to_rank(category1), to_rank(target_category1));
                constexpr rank_t qualified_category2 = std::max(to_rank(category2), to_rank(target_category2));
                return static_cast<proxy_type_category>(std::max(qualified_category1, qualified_category2));
            }
        public:

            template <proxy_type_category category1, proxy_type_category category2>
            static constexpr proxy_type_category value = common_reference_category<category1, category2>();
        };

        template <std::integral key_t, typename alternative_t, typename coverage_t, proxy_type_category category>
        class compressed_sparse_variant_map_proxy
        {
        private:

            static constexpr proxy_type_category category_v = category;

            using ref_position_t = key_t;

            template <std::integral t0, typename t1, typename t2, proxy_type_category>
            friend class compressed_sparse_variant_map_proxy;

            ref_position_t _ref_position;
            alternative_t _alternative;
            coverage_t _coverage;

        public:

            compressed_sparse_variant_map_proxy(ref_position_t ref_position,
                                                alternative_t alternative,
                                                coverage_t coverage) noexcept :
                _ref_position{(ref_position_t &&) ref_position},
                _alternative{(alternative_t &&) alternative},
                _coverage{(coverage_t &&) coverage}
            {}

            compressed_sparse_variant_map_proxy() = default;
            compressed_sparse_variant_map_proxy(compressed_sparse_variant_map_proxy const &) = default;
            compressed_sparse_variant_map_proxy(compressed_sparse_variant_map_proxy &&) = default;
            compressed_sparse_variant_map_proxy &operator=(compressed_sparse_variant_map_proxy const &) = default;
            compressed_sparse_variant_map_proxy &operator=(compressed_sparse_variant_map_proxy &&) = default;

            //  - from reference to value <true, false> -> <false, false>
            //  - from reference to const-reference <true, false> -> <true, true>
            template <typename other_proxy_t>
                requires (!std::same_as<other_proxy_t, compressed_sparse_variant_map_proxy> &&
                          other_proxy_t::category_v == proxy_type_category::reference &&
                          category_v != proxy_type_category::reference)
            compressed_sparse_variant_map_proxy(other_proxy_t const & other) :
                _ref_position(other._ref_position),
                _alternative{other._alternative},
                _coverage{other._coverage}
            {}

            //  - from const-reference to value <true, true> -> <false, false>
            template <typename other_proxy_t>
                requires (!std::same_as<other_proxy_t, compressed_sparse_variant_map_proxy> &&
                          other_proxy_t::category_v == proxy_type_category::const_reference &&
                          category_v == proxy_type_category::object)
            compressed_sparse_variant_map_proxy(other_proxy_t const & other) :
                _ref_position(other._ref_position),
                _alternative{other._alternative},
                _coverage{other._coverage}
            {}

            //  - from value & to reference || const reference: <false, false> -> <true, false> || <true, true>
            template <typename other_proxy_t>
                requires (!std::same_as<std::remove_cvref_t<other_proxy_t>, compressed_sparse_variant_map_proxy> &&
                          other_proxy_t::category_v == proxy_type_category::object &&
                          category_v != proxy_type_category::object)
            compressed_sparse_variant_map_proxy(other_proxy_t & other) :
                _ref_position(other._ref_position),
                _alternative{other._alternative},
                _coverage{other._coverage}
            {}

            //  - from value const & to const reference: <false, false> -> <true, true>
            template <typename other_proxy_t>
                requires (!std::same_as<std::remove_cvref_t<other_proxy_t>, compressed_sparse_variant_map_proxy> &&
                          other_proxy_t::category_v == proxy_type_category::object &&
                          category_v != proxy_type_category::object)
            compressed_sparse_variant_map_proxy(other_proxy_t const & other) :
                _ref_position(other._ref_position),
                _alternative{other._alternative},
                _coverage{other._coverage}
            {}

            explicit constexpr operator std::remove_reference_t<alternative_t> &()
                requires (category_v == proxy_type_category::object) {
                return _alternative;
            }

            explicit constexpr operator std::remove_cvref_t<alternative_t> const &() const
                requires (category_v == proxy_type_category::object) {
                return _alternative;
            }

        private:
            // no equality!
            constexpr bool operator==(compressed_sparse_variant_map_proxy const &) const = delete;

            template <typename other_alternative_t, typename other_coverage_t, proxy_type_category other_category>
            constexpr friend std::weak_ordering operator<=>(compressed_sparse_variant_map_proxy const & lhs,
                                                            compressed_sparse_variant_map_proxy<key_t, other_alternative_t, other_coverage_t, other_category> const & rhs) noexcept {
                if (auto cmp = libjst::position(lhs) <=> libjst::position(rhs); cmp == 0)
                    return libjst::alt_kind(lhs) <=> libjst::alt_kind(rhs);
                else
                    return cmp;
            }

            template <typename this_t, typename member_t>
            using fwd_t = std::conditional_t<std::is_lvalue_reference_v<member_t>,
                                             member_t,
                                             jst::contrib::member_type_t<this_t, std::remove_cvref_t<member_t>>>;

            template <typename this_t>
                requires std::same_as<std::remove_cvref_t<this_t>, compressed_sparse_variant_map_proxy>
            constexpr friend auto tag_invoke(std::tag_t<libjst::position>, this_t &&me) noexcept
                -> fwd_t<this_t, ref_position_t>
            {
                return (fwd_t<this_t, ref_position_t> &&)me._ref_position;
            }

            template <typename this_t>
                requires std::same_as<std::remove_cvref_t<this_t>, compressed_sparse_variant_map_proxy>
            constexpr friend auto tag_invoke(std::tag_t<libjst::coverage>, this_t &&me) noexcept
                -> fwd_t<this_t, coverage_t>
            {
                return (fwd_t<this_t, coverage_t> &&)me._coverage;
            }

            // Forward variant specific calls to stored variant.
            template <typename cpo_t, typename this_t>
                requires std::same_as<std::remove_cvref_t<this_t>, compressed_sparse_variant_map_proxy> &&
                        std::tag_invocable<cpo_t, fwd_t<this_t, alternative_t>>
            constexpr friend auto tag_invoke(cpo_t cpo, this_t &&me)
                noexcept(std::is_nothrow_tag_invocable_v<cpo_t, fwd_t<this_t, alternative_t>>)
                -> std::tag_invoke_result_t<cpo_t, fwd_t<this_t, alternative_t>>
            {
                return std::tag_invoke(cpo, (fwd_t<this_t, alternative_t> &&)me._alternative);
            }
        };

    } // namespace detail

    // we could extract this into a variant store tuple type.
    template <typename alternate_store_t, typename variant_coverage_t>
    class compressed_sparse_variant_map
    {
    private:
        using internal_map_value_type = std::pair<int32_t, std::size_t>;
        using ref_position_map_type = std::vector<internal_map_value_type>; // better exclude into separate class!
        using coverage_store_type = std::vector<variant_coverage_t>;

    public:
        using key_type = typename internal_map_value_type::first_type;
        using mapped_type = detail::compressed_sparse_variant_map_proxy<
                                key_type,
                                std::ranges::range_value_t<alternate_store_t>,
                                std::ranges::range_value_t<coverage_store_type>,
                                detail::proxy_type_category::object>;
        using value_type = mapped_type;
        using reference = detail::compressed_sparse_variant_map_proxy<
                                key_type,
                                std::ranges::range_reference_t<alternate_store_t>,
                                std::ranges::range_reference_t<coverage_store_type>,
                                detail::proxy_type_category::reference>;
        using const_reference = detail::compressed_sparse_variant_map_proxy<
                                    key_type,
                                    std::ranges::range_reference_t<alternate_store_t const>,
                                    std::ranges::range_reference_t<coverage_store_type const>,
                                    detail::proxy_type_category::const_reference>;

        using iterator = variant_store_iterator<compressed_sparse_variant_map>;
        using const_iterator = variant_store_iterator<compressed_sparse_variant_map const>;

        using size_type = typename internal_map_value_type::second_type;
        using difference_type = std::iter_difference_t<const_iterator>;

    private:

        friend iterator;
        friend const_iterator;

        ref_position_map_type _ref_position_map{};
        alternate_store_t _alternatives{};
        coverage_store_type _coverages{};

        static constexpr auto _map_project = [] (internal_map_value_type const & value) constexpr noexcept
            -> key_type const & {
            return value.first;
        };

        static constexpr auto _map_compare = [] (key_type const & lhs, key_type const & rhs) constexpr noexcept
            -> bool {
            return lhs < rhs;
        };

        // ----------------------------------------------------------------------------
        // Access for random access iterator implementation.
        // ----------------------------------------------------------------------------

        constexpr reference operator[](difference_type const offset) noexcept {
            assert(static_cast<size_type>(offset) < size());
            auto [key, index] = _ref_position_map[offset];
            return reference{key, _alternatives[index], _coverages[index]};
        }

        constexpr const_reference operator[](difference_type const offset) const noexcept {
            assert(static_cast<size_type>(offset) < size());
            auto [key, index] = _ref_position_map[offset];
            return const_reference{key, _alternatives[index], _coverages[index]};
        }

    public:
        // ----------------------------------------------------------------------------
        // Construction, assignment and destruction
        // ----------------------------------------------------------------------------
        compressed_sparse_variant_map() = default;

        // ----------------------------------------------------------------------------
        // Capacity
        // ----------------------------------------------------------------------------

        constexpr size_type size() const noexcept {
            return std::ranges::size(_ref_position_map);
        }

        // ----------------------------------------------------------------------------
        // Modifier
        // ----------------------------------------------------------------------------

        constexpr iterator insert(value_type variant) {
            return insert_impl(end(), std::move(variant));
        }
// constexpr libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::iterator
//             libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::insert_impl(
//                 libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::const_iterator,
//                 libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::value_type)

//     alternate_store_t = libjst::single_base_replacement_store<seqan3::dna4>;
//     variant_coverage_t = libjst::bit_vector<>;
//     libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::iterator =
//         libjst::variant_store_iterator<
//             libjst::compressed_sparse_variant_map<libjst::single_base_replacement_store<seqan3::dna4>, libjst::bit_vector<> > >;
//     libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::const_iterator =
//         libjst::variant_store_iterator<
//             const libjst::compressed_sparse_variant_map<libjst::single_base_replacement_store<seqan3::dna4>, libjst::bit_vector<> > >;
//     libjst::compressed_sparse_variant_map<alternate_store_t, variant_coverage_t>::value_type =
//         libjst::detail::compressed_sparse_variant_map_proxy<int, libjst::single_base_replacement_store<seqan3::dna4>::element_type, libjst::bit_vector<>, libjst::detail::proxy_type_category::object>


        // insert with hint!
        constexpr iterator insert(const_iterator hint, value_type variant) {
            return insert_impl(hint, std::move(variant));
        }

        template <typename ...args_t>
            requires std::constructible_from<mapped_type, args_t...>
        iterator emplace_hint(const_iterator hint, args_t &&...args) {
            return insert_hint(hint, mapped_type{(args_t &&)args...});
        }

        template <typename ...args_t>
            requires std::constructible_from<mapped_type, args_t...>
        iterator emplace(args_t &&...args) {
            return insert(mapped_type{(args_t &&)args...});
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        iterator begin() noexcept { return iterator{*this, 0u}; }
        const_iterator begin() const noexcept { return const_iterator{*this, 0u}; }
        iterator end() noexcept { return iterator{*this, size()}; }
        const_iterator end() const noexcept { return const_iterator{*this, size()}; }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_ref_position_map, _alternatives, _coverages);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_ref_position_map, _alternatives, _coverages);
        }

    private:

        constexpr iterator insert_impl(const_iterator hint, value_type variant) {
            auto const mapped_store_idx = size();

            _ref_position_map.reserve(size() + 1);
            _coverages.reserve(size() + 1);
            _alternatives.push_back(std::move(static_cast<std::ranges::range_reference_t<alternate_store_t>>(variant)));
            //------------------------------------------------------------------ noexcept here after
            _coverages.push_back(libjst::coverage(std::move(variant)));
            auto map_it = find_insert_position_near_hint(hint, variant);
            map_it = _ref_position_map.emplace(map_it, libjst::position(std::move(variant)), mapped_store_idx);

            assert(std::ranges::size(_ref_position_map) == size());
            assert(std::ranges::size(_coverages) == size());
            assert(std::ranges::size(_alternatives) == size());

            return std::ranges::next(begin(), std::ranges::distance(std::ranges::begin(_ref_position_map), map_it));
        }

        constexpr auto find_insert_position(value_type const & variant) const noexcept
        {
            key_type const key = libjst::position(variant);
            auto it_key = std::ranges::lower_bound(_ref_position_map, key, _map_compare, _map_project);

            while (it_key != std::ranges::end(_ref_position_map) && _map_project(*it_key) == key) {
                if (libjst::alt_kind(variant) <= libjst::alt_kind(_alternatives[it_key->second]))
                    break;
                ++it_key;
            }
            return it_key;
        }

        constexpr auto find_insert_position_near_hint(const_iterator hint, value_type const & variant) const noexcept {
            bool is_lower_bound = hint == end() || variant <= *hint;
            if (hint != begin())
                is_lower_bound &= *std::ranges::prev(hint) < variant;

            if (is_lower_bound)
                return std::ranges::next(std::ranges::begin(_ref_position_map), std::ranges::distance(begin(), hint));

            return find_insert_position(variant);
        }
    };
}  // namespace libjst

namespace std
{
    template <integral key_t1, typename alternative_t1, typename coverage_t1, libjst::detail::proxy_type_category category1,
              integral key_t2, typename alternative_t2, typename coverage_t2, libjst::detail::proxy_type_category category2,
              template<typename> typename qualifier_t1,
              template<typename> typename qualifier_t2>
        requires requires {
            typename libjst::detail::compressed_sparse_variant_map_proxy<
                        common_type_t<key_t1, key_t2>,
                        common_reference_t<qualifier_t1<alternative_t1>, qualifier_t2<alternative_t2>>,
                        common_reference_t<qualifier_t1<coverage_t1>, qualifier_t2<coverage_t2>>,
                        libjst::detail::common_category_selector<qualifier_t1, qualifier_t2>::template value<category1, category2>>;
        }
    struct basic_common_reference<libjst::detail::compressed_sparse_variant_map_proxy<key_t1, alternative_t1, coverage_t1, category1>,
                                  libjst::detail::compressed_sparse_variant_map_proxy<key_t2, alternative_t2, coverage_t2, category2>,
                                  qualifier_t1,
                                  qualifier_t2>
    {
        using type = libjst::detail::compressed_sparse_variant_map_proxy<
                        common_type_t<key_t1, key_t2>,
                        common_reference_t<qualifier_t1<alternative_t1>, qualifier_t2<alternative_t2>>,
                        common_reference_t<qualifier_t1<coverage_t1>, qualifier_t2<coverage_t2>>,
                        libjst::detail::common_category_selector<qualifier_t1, qualifier_t2>::template value<category1, category2>>;
    };
} // namespace std
