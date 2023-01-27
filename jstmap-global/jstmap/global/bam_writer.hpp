// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bam writer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp> // TODO: fix this!
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

#include <jstmap/global/search_matches.hpp>
#include <jstmap/global/jstmap_types.hpp>

namespace jstmap
{

    class bam_writer {

        using field_ids_type = seqan3::fields<seqan3::field::id,            /*QNAME*/
                                              seqan3::field::ref_id,        /*RNAME*/
                                              seqan3::field::ref_offset,    /*POS*/
                                              seqan3::field::cigar,         /*CIGAR*/
                                              seqan3::field::seq,           /*SEQ*/
                                              seqan3::field::tags           /*OPTIONAL TAGS*/
                                            >;

        using valid_format_type = seqan3::type_list<seqan3::format_bam, seqan3::format_sam>;
        using reference_names_type = std::vector<std::string>;
        using reference_lengths_type = std::vector<std::size_t>;
        using output_file_type = seqan3::sam_file_output<field_ids_type, valid_format_type, reference_names_type>;

        rcs_store_t const & _rcs_store;
        output_file_type _output_file;

    public:
        explicit bam_writer(rcs_store_t const &, std::filesystem::path)
            noexcept(std::is_nothrow_move_constructible_v<std::filesystem::path>);

        void write_matches(search_matches const &);

    private:
        output_file_type create_output_file(std::filesystem::path);
        seqan3::sam_tag_dictionary encode_position(match_position const &) const noexcept;
        void write_program_info() noexcept;
    };
}  // namespace jstmap
