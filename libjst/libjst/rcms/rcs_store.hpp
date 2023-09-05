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

#include <algorithm>

#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/core/concept/cereal.hpp>

#include <libjst/variant/concept.hpp>

namespace libjst
{
    // Client needs to make sure that the value types are compatible!?
    template <seqan3::sequence source_sequence_t, typename cms_t>
    class rcs_store
    {
    private:

        using coverage_domain_type = std::remove_cvref_t<decltype(std::declval<cms_t const &>().coverage_domain())>;

        cms_t _variant_map{};

    public:

        using variant_map_type = cms_t;
        using source_type = typename cms_t::source_type;
        using value_type = std::ranges::range_value_t<cms_t>;
        using reference = std::ranges::range_reference_t<cms_t const &>;
        using size_type = typename coverage_domain_type::value_type;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr rcs_store() = default; //!< Default.
        // Construct from dimensions?
        constexpr rcs_store(source_sequence_t source, size_type initial_row_count) :
            _variant_map{std::move(source), coverage_domain_type{0, initial_row_count}}
        {}
        //!\}

        // what can we do to add some information
        // we may add a sequence and then perform some alignment steps -> via a proxy type
        // we may add a variant using vertical information -> like a list of indices for the coverage?
        // we may add a variant for a single sequence.
        // initial interface: position, alternate_sequence, coverage<>
        constexpr bool add(value_type value) {
        // key_type const src_position,
        //                    alternate_sequence_type variant,
        //                    coverage_type coverage) {
            // assert(src_position.value() + libjst::breakpoint_span(variant) <= std::ranges::size(source()));
            // assert(std::ranges::size(coverage) == size());

            // TODO: coverage interface:
                // .set(position [size/iterator]) -> rather same as select.
                // .unset(position [size/iterator]) -> maybe first some other format before constructing select/ranl support?
                // construction from list of ids.
            // coverage_t coverage{};
            // coverage.resize(_row_count, false);
            // std::ranges::for_each(haplotypes, [&] (auto const id) { coverage[id] = true; });
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
            // mapped_type value{src_position, std::move(variant), std::move(coverage)};
            _variant_map.insert(std::move(value));
            return true;
        }

        constexpr void reserve(size_type const new_capacity)
        {
            _variant_map.reserve(new_capacity);
        }

        // ----------------------------------------------------------------------------
        // Accessor
        // ----------------------------------------------------------------------------

        constexpr source_type source() const noexcept
        {
            return variants().source();
        }

        constexpr variant_map_type const & variants() const noexcept
        {
            return _variant_map;
        }

        constexpr size_type size() const noexcept
        {
            return variants().coverage_domain().size();
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_variant_map);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_variant_map);
        }
    };
}  // namespace libjst
