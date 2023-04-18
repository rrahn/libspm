// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bit coverage.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/core/concept/cereal.hpp>

#include <libjst/coverage/concept.hpp>
#include <libjst/coverage/range_domain.hpp>
#include <libjst/utility/bit_vector.hpp>

namespace libjst
{

    template <std::unsigned_integral value_t>
    class bit_coverage {

        using coverage_domain_t = range_domain<value_t>;
        using domain_value_type = typename coverage_domain_t::value_type;
        using data_type = bit_vector<>;

        data_type _data{};
        [[no_unique_address]] coverage_domain_t _domain{};

        // needs its own iterator and output iterator?

        explicit bit_coverage(data_type data, coverage_domain_t domain) :
            _data{std::move(data)},
            _domain{std::move(domain)}
        {}

    public:
        using value_type = domain_value_type;
        using iterator = std::ranges::iterator_t<data_type const>;

        constexpr bit_coverage() = default;

        explicit constexpr bit_coverage(coverage_domain_t domain) :
            _data{},
            _domain{std::move(domain)}
        {
            if (max_size() > _data.max_size()) {
                throw std::runtime_error{"The domain size is too large for the container type."};
            }
            _data.resize(max_size(), false);
        }

        template <typename elem_range_t>
            requires (!std::same_as<std::remove_cvref_t<elem_range_t>, bit_coverage>) &&
                     (!std::same_as<std::remove_cvref_t<elem_range_t>, std::initializer_list<value_type>>) &&
                      std::integral<std::ranges::range_value_t<elem_range_t>>
        explicit constexpr bit_coverage(elem_range_t && from_list, coverage_domain_t domain) :
            bit_coverage{std::move(domain)}
        {
            std::ranges::for_each(from_list, [&] (auto const & elem) {
                if (!get_domain().is_member(elem))
                    throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};
                _data[elem] = true;
            });
        }

        explicit constexpr bit_coverage(std::initializer_list<value_type> from_list, coverage_domain_t domain) :
            bit_coverage{std::move(domain)}
        {
            std::ranges::for_each(from_list, [&] (auto const & elem) {
                if (!get_domain().is_member(elem))
                    throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};
                _data[elem] = true;
            });
        }



        // explicit constexpr bit_coverage(std::initializer_list<value_type> init_list, coverage_domain_t domain) :
        //     bit_coverage{std::move(domain)}
        // {
        //     std::ranges::for_each(init_list, [&] (value_type const & elem) {
        //         if (!get_domain().is_member(elem))
        //             throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};
        //         _data[elem] = true;
        //     });
        // }


        constexpr std::iter_reference_t<iterator> operator[](std::ptrdiff_t idx) const noexcept {
            return _data[idx];
        }

        constexpr iterator insert(value_type elem) {
            if (!get_domain().is_member(elem))
                throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};
            _data[elem] = true;
            return std::ranges::begin(_data)[elem];
        }

        constexpr void clear() noexcept {
            _data.clear();
        }

        constexpr iterator erase(iterator first) noexcept {
            return erase(first, std::ranges::next(first));
        }

        constexpr iterator erase(iterator first, iterator last) noexcept {
            for (auto it =  first; it != last; ++it) {
                *it = false;
            }
            return last;
        }

        constexpr value_type front() const noexcept {
            assert(!empty());
            return *begin();
        }

        constexpr value_type back() const noexcept {
            assert(!empty());
            return *std::ranges::prev(end());
        }

        constexpr bool empty() const noexcept {
            return _data.none();
        }

        constexpr size_t size() const noexcept {
            return _data.size();
        }

        constexpr size_t max_size() const noexcept {
            return get_domain().size();
        }

        constexpr coverage_domain_t const & get_domain() const noexcept {
            return _domain;
        }

        constexpr iterator begin() const noexcept {
            return _data.begin();
        }

        constexpr iterator end() const noexcept {
            return _data.end();
        }

        constexpr bool any() const noexcept {
            return _data.any();
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_data, _domain);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_data, _domain);
        }

    private:

        constexpr friend bool operator==(bit_coverage const &, bit_coverage const &) noexcept = default;

        constexpr friend bit_coverage
        tag_invoke(std::tag_t<coverage_intersection>, bit_coverage const & first, bit_coverage const & second) {
            // if (first.get_domain() != second.get_domain())
            //     throw std::domain_error{"Trying to intersect elements from different coverage domains."};

            return bit_coverage{first._data & second._data, first.get_domain()};
        }

        constexpr friend bit_coverage
        tag_invoke(std::tag_t<coverage_difference>, bit_coverage first, bit_coverage const & second) {
            // if (first.get_domain() != second.get_domain())
            //     throw std::domain_error{"Trying to intersect elements from different coverage domains."};

            return bit_coverage{first._data.and_not(second._data), first.get_domain()};
        }

    };
}  // namespace libjst
