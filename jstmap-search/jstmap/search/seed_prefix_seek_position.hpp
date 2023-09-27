// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/sequence_tree/seek_position.hpp>
#include <libjst/utility/multi_invocable.hpp>

namespace jstmap
{
    class seed_prefix_seek_position : public libjst::seek_position {
    private:
        using base_t = libjst::seek_position;

    public:

        explicit constexpr seed_prefix_seek_position() = default;

        explicit constexpr seed_prefix_seek_position(base_t seed_position, size_t breakends_count) noexcept : base_t{}
        {
            size_t reverse_breakend_idx = breakends_count - seed_position.get_variant_index();
            seed_position.visit(libjst::multi_invocable{
                [&, reverse_breakend_idx] (libjst::breakpoint_end site) {
                    base_t::reset(reverse_breakend_idx, site);
                },
                [&, reverse_breakend_idx] (libjst::alternate_path_descriptor const &) {
                    base_t::initiate_alternate_node(reverse_breakend_idx + 1);
                }
            });
        }

    private:

        constexpr friend bool operator==(seed_prefix_seek_position const &, seed_prefix_seek_position const &) noexcept = default;

        constexpr friend std::strong_ordering operator<=>(seed_prefix_seek_position const &, seed_prefix_seek_position const &) noexcept = default;
    };

    // template <typename char_t, typename char_traits_t, typename seed_prefix_seek_position>
    //     requires std::same_as<std::remove_cvref_t<seed_prefix_seek_position>, seed_prefix_seek_position>
    // inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
    //                                                               seed_prefix_seek_position && position)
    // {
    //     stream << "<";
    //     position.visit([&] <typename descriptor_t>(descriptor_t const & descriptor) {
    //         if constexpr (std::same_as<descriptor_t, breakpoint_end>) {
    //             stream << "ref = " << ((descriptor == breakpoint_end::low) ? "low" : "high");
    //         } else {
    //             stream << "alt = " << descriptor;
    //         }
    //     });
    //     stream << " variant_idx = " << position.get_variant_index() << ">";
    //     return stream;
    // }
}  // namespace jstmap
