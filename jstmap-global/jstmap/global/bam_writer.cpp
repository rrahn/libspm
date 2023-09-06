// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of bam writer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <iostream>
#include <string>

#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/sequence_tree/path_descriptor.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>

#include <jstmap/global/bam_writer.hpp>
#include <jstmap/global/search_query.hpp>

namespace seqan3 {
    template <> struct sam_tag_type<"ad"_tag> { using type = std::vector<std::byte>; };
    template <> struct sam_tag_type<"rd"_tag> { using type = std::vector<std::byte>; };
    template <> struct sam_tag_type<"lo"_tag> { using type = int32_t; };
} // namespace seqan3

namespace jstmap {


    bam_writer::bam_writer(rcs_store_t const & rcs_store, std::filesystem::path file_name)
        noexcept(std::is_nothrow_move_constructible_v<std::filesystem::path>)
        : _rcs_store{rcs_store},
          _output_file{create_output_file(std::move(file_name))}
    {
        write_program_info();
    }

    bam_writer::output_file_type bam_writer::create_output_file(std::filesystem::path file_name)
    {
        using namespace std::literals;
        reference_names_type reference_names{"referentially compressed sequence store"s};
        std::vector<std::size_t> reference_lengths{_rcs_store.variants().size()};
        return output_file_type{std::move(file_name), std::move(reference_names), std::move(reference_lengths)};
    }

    void bam_writer::write_matches(search_matches const & query_matches)
    {
        using namespace std::literals;
        size_t cnt{};
        for (auto const & match : query_matches.matches()) {
            _output_file.emplace_back(query_matches.query().value().id(),                   /*QNAME*/
                                      _output_file.header().ref_ids()[0],                   /*RNAME*/
                                      match.position().tree_position.get_variant_index(),   /*POS*/
                                      match.get_cigar(),                                    /*CIGAR*/
                                      query_matches.query().value().sequence(),             /*SEQ*/
                                      encode_position(match.position())                     /*OPTIONAL TAGS*/
                                    );
        }
    }

    seqan3::sam_tag_dictionary bam_writer::encode_position(match_position const & position) const noexcept {
        using namespace seqan3::literals;

        seqan3::sam_tag_dictionary dict{};

        // conversion to byte value?
        position.tree_position.visit(seqan3::detail::multi_invocable {
            [&](libjst::alternate_path_descriptor descriptor) {
                size_t max_byte_count = (descriptor.size() + 7) / 8;

                std::vector<std::byte> tmp(max_byte_count, std::byte{});
                std::memcpy(tmp.data(), descriptor.data(), tmp.size());
                dict.get<"ad"_tag>() = std::move(tmp);
            },
            [&](libjst::breakpoint_end descriptor) {
                using data_t = std::underlying_type_t<libjst::breakpoint_end>;
                static constexpr size_t max_byte_count = sizeof(data_t);

                std::vector<std::byte> tmp(max_byte_count, std::byte{});
                data_t _data = static_cast<data_t>(static_cast<libjst::breakpoint_end>(descriptor));
                std::memcpy(tmp.data(), &_data, tmp.size());
                dict.get<"rd"_tag>() = std::move(tmp);
            }
        });
        dict.get<"lo"_tag>() = position.label_offset;
        return dict;
    }

    void bam_writer::write_program_info() noexcept
    {
        using namespace std::literals;
        using header_t = std::remove_reference_t<decltype(_output_file.header())>;
        using program_info_t = typename header_t::program_info_t;

        _output_file.header().program_infos.push_back(program_info_t{
            .name = "jst tools"s,
            .command_line_call = "add program call"s,
            .description = "Generated from the jst tools"s,
            .version = "0.0.1"s
        });
    }

}  // namespace jstmap
