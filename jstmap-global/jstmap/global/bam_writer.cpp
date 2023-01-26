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

#include <jstmap/global/bam_writer.hpp>
#include <jstmap/global/search_query.hpp>

namespace jstmap
{

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

    void bam_writer::write_matches(all_matches const & query_matches)
    {
        using namespace std::literals;
        size_t cnt{};
        for (match_position const & pos : query_matches.matches()) {
            _output_file.emplace_back(query_matches.query().value().id(),        /*QNAME*/
                                      _output_file.header().ref_ids()[0],        /*RNAME*/
                                      pos.tree_position.get_variant_index(),     /*POS*/
                                      query_matches.query().value().sequence()   /*SEQ*/
                                    );
        }
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
