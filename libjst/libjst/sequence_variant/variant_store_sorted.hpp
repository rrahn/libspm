// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>
#include <vector>

#include <cereal/types/vector.hpp>

#include <libcontrib/type_traits.hpp>

#include <libjst/sequence_variant/variant_store_iterator.hpp>

namespace libjst
{
    struct variant_less
    {
        template <sequence_variant lhs_t, sequence_variant rhs_t>
        constexpr bool operator()(lhs_t const &lhs, rhs_t const &rhs) const noexcept
        {
            auto effective_size = [](auto const &variant) {
                return std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
            };

            return (libjst::position(lhs) < libjst::position(rhs)) ||
                    (libjst::position(lhs) == libjst::position(rhs) && effective_size(lhs) > effective_size(rhs));
        }
    };

    template <typename variant_store_t, typename compare_t = variant_less>
    class variant_store_sorted
    {
    private:
        using store_t = variant_store_t;

        using value_type = std::ranges::range_value_t<variant_store_t>;
        using reference = std::ranges::range_reference_t<variant_store_t>;
        using const_reference = std::ranges::range_reference_t<variant_store_t const>;

        using position_t = variant_position_t<value_type>;
        using events_t = std::vector<position_t>;
        using iterator = variant_store_iterator<variant_store_sorted>;
        using const_iterator = variant_store_iterator<variant_store_sorted const>;

        friend iterator;
        friend const_iterator;

        std::reference_wrapper<variant_store_t> _store{};
        events_t _event_queue{};

    public:

        variant_store_sorted() = delete;
        explicit variant_store_sorted(variant_store_t &store) : _store{store}
        {
            auto & variant_store = _store.get();
            position_t store_size = std::ranges::size(variant_store);
            _event_queue.resize(store_size);
            std::ranges::copy(std::views::iota(0u, store_size), std::ranges::begin(_event_queue));

            std::ranges::stable_sort(_event_queue.begin(), _event_queue.end(), compare_t{}, [&] (position_t const idx) {
                return variant_store[idx];
            }); // now sort based on the indices.
        }

        // ----------------------------------------------------------------------------
        // Access
        // ----------------------------------------------------------------------------

        constexpr reference operator[](size_t offset) noexcept {
            assert(offset < _event_queue.size());
            return _store.get()[_event_queue[offset]];
        }

        constexpr const_reference operator[](size_t offset) const noexcept {
            assert(offset < _event_queue.size());
            return _store.get()[_event_queue[offset]];
        }

        constexpr size_t size() const noexcept {
            return _event_queue.size();
        }

        // ----------------------------------------------------------------------------
        // Iterator
        // ----------------------------------------------------------------------------

        iterator begin() noexcept { return iterator{*this, 0u}; }
        const_iterator begin() const noexcept { return const_iterator{*this, 0u}; }
        iterator end() noexcept { return iterator{*this, size()}; }
        const_iterator end() const noexcept { return const_iterator{*this, size()}; }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_event_queue);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_event_queue);
        }
    };

}  // namespace libjst
