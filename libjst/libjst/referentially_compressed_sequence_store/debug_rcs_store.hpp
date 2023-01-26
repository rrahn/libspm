// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides rooted rcs store implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

#include <seqan3/utility/views/slice.hpp>

namespace libjst
{

    template <typename rcs_store_t>
        requires std::common_reference_with<
            std::ranges::range_value_t<typename rcs_store_t::variant_map_type> &,
            std::ranges::range_reference_t<typename rcs_store_t::variant_map_type>>
    class debug_rcs_store
    {
    private:
        using base_variant_map_type = typename rcs_store_t::variant_map_type;
        using sliced_variant_map_t = decltype(std::declval<base_variant_map_type const &>() | seqan3::views::slice(0, 1));

        using base_source_type = typename rcs_store_t::source_type;
        using sliced_source_type = decltype(std::declval<base_source_type const &>() | seqan3::views::slice(0, 1));

        rcs_store_t const & _wrappee;
        sliced_source_type _sliced_source;
        sliced_variant_map_t _sliced_variants;

    public:

        using variant_map_type = sliced_variant_map_t;
        using source_type = sliced_source_type;

        debug_rcs_store() = delete;
        explicit debug_rcs_store(rcs_store_t const & wrappee,
                                 std::pair<size_t, size_t> source_slice,
                                 std::pair<size_t, size_t> variants_slice) :
            _wrappee{wrappee},
            _sliced_source{_wrappee.source() | seqan3::views::slice(source_slice.first, source_slice.second)},
            _sliced_variants{_wrappee.variants() | seqan3::views::slice(variants_slice.first, variants_slice.second)}
        {
            for (auto it = _sliced_variants.begin(); it != _sliced_variants.end() - 1; ++it) {
                auto next = it + 1;

                if (libjst::right_breakpoint(*it) > libjst::left_breakpoint(*next)) {
                    if (auto cov = libjst::coverage(*it) & libjst::coverage(*next); cov.any()) {
                        std::cout << "POS: " << (it - _sliced_variants.begin()) << "\n";
                        std::cout << "LBP fst: " << libjst::left_breakpoint(*it) << "\n";
                        std::cout << "LBP snd: " << libjst::left_breakpoint(*next) << "\n";
                        seqan3::debug_stream << "Cov: " << cov << "\n";
                        throw std::runtime_error{"Variants not properly sorted!\n"};
                    }
                }
            }
        }

        constexpr rcs_store_t const & base() const noexcept {
            return _wrappee;
        }

        constexpr auto size() const noexcept -> decltype(_wrappee.size()){
            return _wrappee.size();
        }

        constexpr sliced_source_type const & source() const noexcept {
            return _sliced_source;
        }

        constexpr variant_map_type const & variants() const noexcept {
            return _sliced_variants;
        }

    };

    template <typename rcs_store_t>
    debug_rcs_store(rcs_store_t const &, size_t const, size_t const) -> debug_rcs_store<rcs_store_t>;

}  // namespace libjst
