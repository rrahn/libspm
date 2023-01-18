// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides haplotype viewer for rcs store
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

#include <libjst/journal.hpp>
#include <libjst/referentially_compressed_sequence_store/rcs_store.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <typename rcs_store_t>
    class haplotype_viewer
    {
    private:

        using source_type = typename rcs_store_t::source_type;

        using variant_map_type = typename rcs_store_t::variant_map_type;
        using variant_type = std::ranges::range_value_t<variant_map_type>;
        using position_type = std::remove_cvref_t<libjst::variant_position_t<variant_type>>;
        using difference_type = typename position_type::value_type;

        using journal_type = journal<difference_type, source_type const &>;

        class proxy;

        std::reference_wrapper<rcs_store_t const>  _wrappee;

    public:

        using reference = proxy;

        haplotype_viewer() = delete;
        explicit haplotype_viewer(rcs_store_t const & wrappee) : _wrappee{wrappee}
        {}

        constexpr reference operator[](difference_type const offset) const noexcept {
            return proxy{*this, offset};
        }

        constexpr rcs_store_t const & base() const noexcept {
            return _wrappee.get();
        }
    };

    template <typename rcs_store_t>
    haplotype_viewer(rcs_store_t const &) -> haplotype_viewer<rcs_store_t>;

    template <typename rcs_store_t>
    class haplotype_viewer<rcs_store_t>::proxy : std::ranges::view_base {
    private:

        friend haplotype_viewer;

        journal_type _journal{};

        constexpr explicit proxy(haplotype_viewer const & host, difference_type const offset) :
            _journal{host.base().source()}
        {

            if (offset >= 0 && offset < host.base().size()) {
                std::ptrdiff_t journal_offset{};
                std::ranges::for_each(host.base().variants(), [&] (auto const & variant) {
                    if (libjst::coverage(variant)[offset]) {
                        record(variant, journal_offset + libjst::position(variant));
                        journal_offset += libjst::effective_size(variant);
                    }
                });
            }
        }
    public:
        proxy() = default;

        constexpr auto begin() const noexcept {
            return std::ranges::begin(_journal.sequence());
        }

        constexpr auto end() const noexcept {
            return std::ranges::end(_journal.sequence());
        }

    private:

        constexpr void record(variant_type const & variant, std::size_t position) {
            switch (libjst::alt_kind(variant)) {
                case alternate_sequence_kind::replacement: {
                    _journal.record_substitution(position, libjst::alt_sequence(variant));
                    break;
                } case alternate_sequence_kind::deletion: {
                    _journal.record_deletion(position, libjst::breakpoint_span(variant));
                    break;
                } case alternate_sequence_kind::insertion: {
                    _journal.record_insertion(position, libjst::alt_sequence(variant));
                    break;
                } default: {
                    throw std::runtime_error{"Unknown alternate sequence kind."};
                }
            }
        }
    };

}  // namespace libjst


