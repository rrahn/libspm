// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a union over multiple variant stores.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>

#include <cereal/types/vector.hpp>

#include <libjst/sequence_variant/concept.hpp>
#include <libjst/sequence_variant/variant_store_iterator.hpp>

namespace libjst
{
    template <typename _variant_t, typename _coverage_t>
    class variant_proxy
    {
    private:
        _variant_t _variant;
        _coverage_t _coverage;

    public:
        variant_proxy() = default;
        variant_proxy(variant_proxy const &) = default;
        variant_proxy(variant_proxy &&) = default;
        variant_proxy &operator=(variant_proxy const &) = default;
        variant_proxy &operator=(variant_proxy &&) = default;

        // explicit construction from value types
        explicit variant_proxy(_variant_t variant, _coverage_t coverage) noexcept :
            _variant{(_variant_t &&) variant},
            _coverage{(_coverage_t &&) coverage}
        {}

        template <typename other_variant_t, typename other_coverage_t>
            requires (!std::same_as<variant_proxy<other_variant_t, other_coverage_t>, variant_proxy> &&
                      std::constructible_from<_variant_t, other_variant_t const &> &&
                      std::constructible_from<_coverage_t, other_coverage_t const &>)
        variant_proxy(variant_proxy<other_variant_t, other_coverage_t> const &other) :
            _variant{other._variant},
            _coverage{other._coverage}
        {}

        template <typename other_variant_t, typename other_coverage_t>
            requires (!std::same_as<variant_proxy<other_variant_t, other_coverage_t>, variant_proxy> &&
                      std::constructible_from<_variant_t, other_variant_t &> &&
                      std::constructible_from<_coverage_t, other_coverage_t &>)
        variant_proxy(variant_proxy<other_variant_t, other_coverage_t> &other) :
            _variant{other._variant},
            _coverage{other._coverage}
        {}

        template <typename other_variant_t, typename other_coverage_t>
            requires (!std::same_as<variant_proxy<other_variant_t, other_coverage_t>, variant_proxy> &&
                      std::constructible_from<_variant_t, other_variant_t> &&
                      std::constructible_from<_coverage_t, other_coverage_t>)
        variant_proxy(variant_proxy<other_variant_t, other_coverage_t> const &&other) :
            _variant{std::move(other._variant)},
            _coverage{std::move(other._coverage)}
        {}

        _variant_t & get() noexcept
        {
            return _variant;
        }

        _variant_t const & get() const noexcept
        {
            return _variant;
        }

    private:
        constexpr friend _coverage_t const & tag_invoke(std::tag_t<libjst::coverage>, variant_proxy const &me) noexcept
        {
            return me._coverage;
        }

        template <typename cpo_t>
            requires std::tag_invocable<cpo_t, _variant_t const &>
        constexpr friend auto tag_invoke(cpo_t cpo, variant_proxy const &me)
            noexcept(std::is_nothrow_tag_invocable_v<cpo_t, _variant_t const &>)
            -> std::tag_invoke_result_t<cpo_t, _variant_t const &>
        {
            return std::tag_invoke(cpo, me._variant);
        }
    };

    // we could extract this into a variant store tuple type.
    template <typename variant_store_t, typename coverage_t>
    class variant_store_covered : public variant_store_t
    {
    private:
        using base_t = variant_store_t;

        using value_type = variant_proxy<std::ranges::range_value_t<base_t>, coverage_t>; // TODO: maybe some thing different?
        using reference = variant_proxy<std::ranges::range_reference_t<base_t const &>, coverage_t const &>;
        using iterator = variant_store_iterator<variant_store_covered>;

        friend iterator;

        std::vector<coverage_t> _coverage{};

    public:
        // ----------------------------------------------------------------------------
        // Construction, assignment and destruction
        // ----------------------------------------------------------------------------
        variant_store_covered() = default;

        // ----------------------------------------------------------------------------
        // Access
        // ----------------------------------------------------------------------------

        constexpr reference operator[](size_t offset) const noexcept {
            assert(offset < _coverage.size());
            return reference{static_cast<base_t const &>(*this)[offset], _coverage[offset]};
        }

        // ----------------------------------------------------------------------------
        // Modifier
        // ----------------------------------------------------------------------------

        iterator insert(value_type covered_variant) // may throw
        {
            // we do not know the insert position in the underlying, so how do we connect it?
            _coverage.reserve(_coverage.size() + 1); // may throw
            auto && tmp_coverage = std::move(libjst::coverage(covered_variant));
            auto inserted_it = base_t::insert(std::move(covered_variant.get())); // may throw
            assert(inserted_it != base_t::end());
            //------------------------------------------------------------------ noexcept here after
            auto position = inserted_it - base_t::begin();
            _coverage.insert(std::ranges::next(_coverage.begin(), position), std::move(tmp_coverage));
            assert(variant_store_t::size() == _coverage.size());
            return std::ranges::next(begin(), position);
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        iterator begin() const noexcept { return iterator{*this, 0u}; }
        iterator end() const noexcept { return iterator{*this, base_t::size()}; }

        template <seqan3::cereal_archive archive_t>
        void serialize(archive_t & ioarchive)
        {
            ioarchive(cereal::base_class<base_t>( this ), _coverage);
        }
    };
}  // namespace libjst

namespace std
{
    template <typename variant_t, typename coverage_t,
              typename variant_u, typename coverage_u,
              template<typename> typename qualifier_t, template<typename> typename qualifier_u>
        requires requires {
            typename libjst::variant_proxy<std::common_reference_t<qualifier_t<variant_t>, qualifier_u<variant_u>>,
                                           std::common_reference_t<qualifier_t<coverage_t>, qualifier_u<coverage_u>>>;
        }
    struct basic_common_reference<libjst::variant_proxy<variant_t, coverage_t>,
                                  libjst::variant_proxy<variant_u, coverage_u>,
                                  qualifier_t, qualifier_u>
    {
        using type = libjst::variant_proxy<std::common_reference_t<qualifier_t<variant_t>, qualifier_u<variant_u>>,
                                           std::common_reference_t<qualifier_t<coverage_t>, qualifier_u<coverage_u>>>;
    };
} // namespace std
