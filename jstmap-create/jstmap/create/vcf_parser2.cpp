// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <algorithm>
#include <charconv>
#include <chrono>
#include <ranges>
#include <span>
#include <string>

#include <seqan/vcf_io.h>

#include <seqan3/io/record.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/create/vcf_parser.hpp>
#include <jstmap/create/stripped_vcf_record.hpp>

#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_direct.hpp>
#include <libjst/journaled_sequence_tree/serialiser_delegate.hpp>
#include <libjst/sequence_variant/variant_store_sorted.hpp>

namespace jstmap
{
raw_sequence_t load_base_sequence(std::filesystem::path const & reference_file, std::string_view contig_name)
{
    using namespace std::literals;
    seqan3::sequence_file_input<sequence_input_traits> record_contig_names{reference_file};

    std::string chr_id = /*"chr"s + */std::string{contig_name};
    for (auto && record : record_contig_names) {
        if (seqan3::get<seqan3::field::id>(record).starts_with(chr_id))
            return seqan3::get<seqan3::field::seq>(record);
    }
    throw std::runtime_error{"Could not find a contig with the name <"s + chr_id + ">!"s};
    return {};
}

void construct_jst_from_vcf2(std::filesystem::path const & reference_file,
                             std::filesystem::path const & vcf_file_path,
                             std::filesystem::path const & out_file_path)
{
    // Get the application logger.
    auto & log = get_application_logger();

    auto duration = [] (auto const & start) {
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() -
                    start).count();
    };

    // ----------------------------------------------------------------------------
    // Parse the vcf file.
    // ----------------------------------------------------------------------------

    log(verbosity_level::verbose, logging_level::info, "Initialise parsing vcf file ", vcf_file_path);

    // ----------------------------------------------------------------------------
    // Open vcf file handle.
    seqan::VcfFileIn vcf_file{vcf_file_path.c_str()};
    seqan::VcfHeader vcf_header{};
    seqan::readHeader(vcf_header, vcf_file);

    if (seqan::atEnd(vcf_file))
    {
        log(verbosity_level::standard,
            logging_level::warning,
            "The vcf file ", vcf_file_path, " does not contain any records!");
        return;
    }

    // ----------------------------------------------------------------------------
    // Transform vcf record into intermediate variants.

    auto start = std::chrono::high_resolution_clock::now();

    // we first build the composite variants!
    // std::vector<stripped_vcf_record::intermediate_variant_t> intermediate_variant_store{};
    variant_store_t intermediate_variant_store{};

    seqan::VcfRecord record{};
    while (!seqan::atEnd(vcf_file))
    {
        seqan::readRecord(record, vcf_file);
        stripped_vcf_record tmp_record{record, seqan::context(vcf_file)};
        tmp_record.alternatives(intermediate_variant_store);
        // std::ranges::copy(tmp_record.alternatives(), std::back_inserter(intermediate_variant_store));
    }

    log(verbosity_level::verbose,
        logging_level::info,
        "Time parsing vcf: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Sort the variants.

    // different model representations -> so the forward_jst is over the reference and the sorted objects model
    start = std::chrono::high_resolution_clock::now();
    libjst::variant_store_sorted<variant_store_t> sorted_store{intermediate_variant_store};
    log(verbosity_level::verbose,
        logging_level::info,
        "Time sorting jst: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Resolve ambiguities.
    start = std::chrono::high_resolution_clock::now();

    auto first = std::ranges::begin(sorted_store);
    auto last = std::ranges::end(sorted_store);

    auto resolve_coverage = [] (auto first, auto last, auto && cmp)
    {
        auto overlap_end = std::ranges::find_if(first, last, std::forward<decltype(cmp)>(cmp));
        for (; first != overlap_end; ++first) {
            for (auto it = std::ranges::next(first, 1, overlap_end); it != overlap_end; ++it) {
                    (libjst::coverage(*it)).and_not(libjst::coverage(*first));
            }
        }
        return first;
    };

    auto effective_size = [](auto const &variant) {
        return std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
    };

    while (first != last) {
        // find first adjacent position overlap
        first = std::ranges::adjacent_find(first, last, [] (auto const & lhs, auto const & rhs) {
            return libjst::position(lhs) == libjst::position(rhs);
        });

        // Resolve overlaps of insertions at current position overlap
        first = resolve_coverage(first, last, [&](auto const & variant) {
            return libjst::position(variant) != libjst::position(*first) || effective_size(variant) <= 0;
        });

        // Resolve overlaps of remaining variants at current position overlap
        first = resolve_coverage(first, last, [&](auto const & variant) {
            return libjst::position(variant) != libjst::position(*first);
        });
    }

    log(verbosity_level::verbose,
        logging_level::info,
        "Time resolving ambiguities in store: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Generate jst model
    // ----------------------------------------------------------------------------

    start = std::chrono::high_resolution_clock::now();
    // Load reference sequence!
    auto base_sequence = load_base_sequence(reference_file,
                                            seqan::toCString(seqan::contigNames(seqan::context(vcf_file))[record.rID]));

    jst_model_t jst{base_sequence, std::move(intermediate_variant_store)};

    log(verbosity_level::verbose,
        logging_level::info,
        "Time generating jst: ", duration(start), " s");

    start = std::chrono::high_resolution_clock::now();
    fwd_jst_t fwd_jst{jst};

    log(verbosity_level::verbose,
        logging_level::info,
        "Time building forward layer: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Serialise jst
    // ----------------------------------------------------------------------------

    start = std::chrono::high_resolution_clock::now();

    std::ofstream archive_stream{out_file_path, std::ios_base::out | std::ios_base::binary};
    {
        cereal::BinaryOutputArchive output_archive(archive_stream);
        auto arch = output_archive | libjst::direct_serialiser(base_sequence)
                                   | libjst::delegate_serialiser(jst);
        libjst::save(fwd_jst, arch);
    }

    log(verbosity_level::verbose,
        logging_level::info,
        "Time serialising jst: ", duration(start), " s");
}

}  // namespace jstmap
