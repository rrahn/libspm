// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides referentially compressed sequence store.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithms>
#include <seqan3/range/concept.hpp>

#include <libjst/utility/bit_vector.hpp>
#include <libjst/variant/alternate_sequence_kind.hpp>
#include <libjst/variant/compressed_sparse_variant_map.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{
    // Client needs to make sure that the value types are compatible!?
    template <seqan3::sequence source_sequence_t, typename alternate_sequence_store_t>
        // requires sequence_alternative<std::ranges::range_value_t<alternate_sequence_store_t>>
        // requires std::ranges::sized_range<source_sequence_t>
    class rcs_store
    {
    public:

        using coverage_type = bit_vector<>;
        using variant_map_type = compressed_sparse_variant_map<alternate_sequence_store_t, coverage_type>;
        using source_type = source_sequence_t;

        using key_type = typename variant_map_type::key_type;

        using size_type = decltype(std::ranges::size(std::devlcal<coverage_type const &>()));
        using alternate_sequence_type = std::ranges::range_value_t<alternate_sequence_store_t>;

    private:

        source_sequence_t _source{};
        variant_map_type _variant_map{};
        size_type _row_count{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr rcs_store() = default; //!< Default.
        // Construct from dimensions?
        constexpr rcs_store(source_sequence_t source, size_type initial_row_count) :
            _source{std::move(source)},
            _row_count{initial_row_count}
        {}
        //!\}

        // what can we do to add some information
        // we may add a sequence and then perform some alignment steps -> via a proxy type
        // we may add a variant using vertical information -> like a list of indices for the coverage?
        // we may add a variant for a single sequence.
        // initial interface: position, alternate_sequence, coverage<>
        template <std::forward_range haplotypes_t>
            requires std::integral<std::ranges::range_value_t<haplotypes_t>>
        constexpr bool add(size_type const src_position,
                           alternate_sequence_type alt_sequence,
                           haplotypes_t && haplotypes) {
            assert(src_position + libjst::breakpoint_span(alt_sequence) <= std::ranges::size(_source_sequence));
            assert([&] { return std::ranges::all_of(haplotypes, [&] (auto const id) { return id < _row_count; }); }());

            // TODO: coverage interface:
                // .set(position [size/iterator]) -> rather same as select.
                // .unset(position [size/iterator]) -> maybe first some other format before constructing select/ranl support?
                // construction from list of ids.
            coverage_t coverage{};
            coverage.resize(_row_count, false);
            std::ranges::for_each(haplotypes, [&] (auto const id) { coverage[id] = true; });
            // How to handle certain types of collisions? -> provide various sets of proxies -> or as a policy design.
            // default is hard fail in case variant is there and
            // how can we make sure that the same variant is not invalidated by some other variant?
            // e.g. long deletion covering the current variant at the same row?

            // we can search the map from outside since we offer a RA iterator, but then we load the coverage everytime we access an element.
            // using search_key = std::pair<size_type, alternate_sequence_kind>;
            // auto map_it = std::ranges::lower_bound(_variant_map,
            //                                        search_key{src_position, libjst::alt_kind(alt_sequence)},
            //                                        [] (auto && variant, search_key const & key) {
            //                                             if (key.first == libjst::position(variant))
            //                                                 return key.second < libjst::alt_kind(variant);
            //                                             else
            //                                                 return key.first < libjst::position(variant);
            //                                        });
            _variant_map.emplace(src_position, std::move(alt_sequence), std::move(coverage));
            return true;
        }

        // ----------------------------------------------------------------------------
        // Accessor
        // ----------------------------------------------------------------------------

        constexpr source_sequence_t const & source_sequence() const noexcept
        {
            return _source;
        }

        constexpr variant_map_type const & sequence_variants() const noexcept
        {
            return _variant_map;
        }

        constexpr size_type size() const noexcept
        {
            return _row_count;
        }

    };
}  // namespace libjst
