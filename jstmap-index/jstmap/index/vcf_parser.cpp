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
#include <ranges>
#include <span>

#include <seqan/vcf_io.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/to.hpp>

#include <jstmap/index/application_logger.hpp>
#include <jstmap/index/vcf_parser.hpp>

#include <libjst/detail/delta_event_shared.hpp>

namespace jstmap
{

//!\brief An augmented vcf record which extracts the genotype infos as a libjst::detail::delta_event_shared.
class augmented_vcf_record
{
    //!\brief The internally used shared delta event type.
    using shared_event_type = libjst::detail::delta_event_shared<seqan3::dna5>;
    //!\brief The pure delta event type.
    using event_type = typename shared_event_type::delta_event_type;
    //!\brief The substitution type.
    using substitution_type = typename shared_event_type::substitution_type;
    //!\brief The deletion type.
    using deletion_type = typename shared_event_type::deletion_type;
    //!\brief The insertion type.
    using insertion_type = typename shared_event_type::insertion_type;
    //!\brief The coverage type.
    using coverage_type = typename shared_event_type::coverage_type;

    //!\brief The header stored with the record.
    seqan::VcfHeader _header{};
    //!\brief The actual record which is augmented.
    seqan::VcfRecord _record{};

    //!\brief How many samples are represented.
    size_t _sample_count{0};
    //!\brief How many haplotypes per sample are present.
    size_t _haplotype_per_sample_count{0};
    //!\brief How total haplotype count.
    size_t _haplotype_count{0};
    //!\brief If this record was initialised already.
    bool _is_initialised{false};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    augmented_vcf_record() = default;
    augmented_vcf_record(seqan::VcfHeader && header) : _header{std::move(header)}
    {}
    //!\}

    //!\brief Returns a reference to the read record.
    seqan::VcfRecord & seqan_record() noexcept
    {
        return _record;
    }

    //!\brief overload
    seqan::VcfRecord const & seqan_record() const noexcept
    {
        return _record;
    }

    //!\brief Initialises the sample and haplotype counts from the first record.
    void initialise_counts()
    {
        if (!_is_initialised)
        {
            _sample_count = sample_count();
            _haplotype_per_sample_count = haplotypes_per_sample_count();
            _haplotype_count = _sample_count * _haplotype_per_sample_count;
            _is_initialised = true;
        }
    }

    //!\brief Returns the total number of haplotypes.
    size_t haplotype_count() const noexcept
    {
        assert(_is_initialised);
        return _haplotype_count;
    }

    //!\brief Returns the reference position stored inside of the record.
    size_t reference_position() const noexcept
    {
        return _record.beginPos;
    }

    //!\brief Returns the variant identifier stored for the record.
    auto variant_identifier() const noexcept
    {
        auto * data_ptr = std::addressof(_record.id[0]);
        return std::span{data_ptr, seqan::length(_record.id)};
    }

    //!\brief Returns the chromosome id stored for the record.
    size_t chromosome_id() const noexcept
    {
        return _record.rID;
    }

    //!\brief Generates and returns the delta events for this record.
    [[nodiscard]] auto generate_delta_events() const
        -> std::pair<bool, std::vector<shared_event_type>>
    {
        if (!genotype_info_given())
            return {false, {}}; // break here if there is no genotype info.

        auto [has_valid_alternatives, delta_events] = extract_delta_events();
        size_t const record_alternative_count = delta_events.size();

        bool is_invalid_record = !has_valid_alternatives;

        // Create vector with coverages to fill.
        std::vector<coverage_type> coverage_per_alternative(record_alternative_count, coverage_type(_haplotype_count));

        size_t haplotype_idx{};
        // Iterate over each contained sample.
        for (auto const & sample_genotype_info : _record.genotypeInfos)
        {
            // Extract the sample genotype info.
            auto sample_genotype = split_by_delimiter(to_span(sample_genotype_info), ':');
            // Iterate over the haplotypes per sample.
            for (size_t const alt_idx : extract_alternative_indices(sample_genotype.front()))
            {
                if (haplotype_idx > _haplotype_count)
                    return {false, {}};

                assert(alt_idx <= record_alternative_count); // Must be a valid id into the delta events vector.
                assert(haplotype_idx < _haplotype_count); // Must not exceed the global haplotype count.

                if (alt_idx > 0) // skip if it is a reference allele
                    coverage_per_alternative[alt_idx - 1][haplotype_idx] = true;

                ++haplotype_idx; // Go to next global haplotype.
            }
        }

        if (haplotype_idx < _haplotype_count)
            return {false, {}};

        std::vector<shared_event_type> shared_events(record_alternative_count);
        for (size_t idx = 0; idx < shared_events.size(); ++idx)
        {
            shared_events[idx] = shared_event_type{std::move(delta_events[idx]),
                                                   std::move(coverage_per_alternative[idx])};
        }

        return {!is_invalid_record, shared_events};
    }

private:

    //!\brief Helper function to split some structural information by a given delimiter.
    auto split_by_delimiter(std::span<const char> source_range, char const delimiter) const
        -> std::vector<std::span<const char>>
    {
        std::vector<std::span<const char>> split_range{};

        auto segment_it = std::ranges::begin(source_range);
        auto source_end = std::ranges::end(source_range);
        while (segment_it != source_end)
        {
            auto segment_end = std::find(segment_it, source_end, delimiter);
            const char * segment_ptr = std::addressof(*segment_it);
            split_range.emplace_back(segment_ptr, segment_ptr + std::distance(segment_it, segment_end));
            segment_it = std::ranges::next(segment_end, 1, source_end);
        }

        return split_range;
    }

    //!\brief Extracts the indices of the alternatives stored inside of one genotype information.
    auto extract_alternative_indices(std::span<const char> genotype) const
        -> std::vector<size_t>
    {
        std::vector<size_t> alt_ids{};
        const char * first = genotype.data();
        const char * last = genotype.data() + genotype.size();

        while (first != last)
        {
            size_t alt_id;
            if (auto [next, ec] = std::from_chars(first, last, alt_id); ec == std::errc())
            {
                alt_ids.push_back(alt_id);
                assert(next == last || *next == '|');
                first = std::ranges::next(next, 1, last);
            }
            else
            {
                break;
            }
        }
        return alt_ids;
    }

    //!\brief Extracts the delta events of the stored genotype infos.
    auto extract_delta_events() const
        -> std::pair<bool, std::vector<event_type>>
    {
        using namespace std::literals;

        std::vector<event_type> delta_events;
        // ALT either a sequence or a id?
        // ALT value: comma separated list, each value can represent "ACGTN*",
            // * means missing allele due to upstream deletion. What is that?
            // we skip this but record it depending on the warning level.
            // ACGTN: means we use dna5 as refernce symbols. ?How to process N?

        // REF: "ACGTN" -> if we read something else, then we need to throw.
            // POS is reference_position of first base in reference string
            // always prepends at least some bases to avoid empty indels in alternatives
            // if event occurs at reference_position one
                // must include base after the event;
                // padding base not required for complex substitutions where all alleles have at least one base represented in their strings!
            // if ALT allele is symbolic? allele ('<ID>') then padding base is required and POS denotes the
                    // coordinate of the base preceding the polymorphsim
            // case must not be preserved

        // Do not process SVs or invalid alternative data.
        if (seqan::empty(_record.alt) || seqan::empty(_record.ref) || _record.alt[0] == '*' || _record.alt[0] == '<')
            return {false, delta_events};

        auto alternatives = split_by_delimiter(to_span(_record.alt), ',');
        assert(!alternatives.empty());

        std::span reference_segment = to_span(_record.ref);

        for (auto const & alternative : alternatives) // we can now process all alternatives.
        {
            // Alternative with no allele information
            if (std::ranges::equal(alternative, "*"sv))
                continue;

            if (alternative.front() == '<')
            {
                if (alternative.back() != '>')
                    throw std::domain_error{"Unexpected symbolic allele terminator: "s +
                                            std::string{alternative.data(), alternative.size()} + "!"s};

                // TODO: Process symbolic allele
                continue;
            }

            // TODO: handle breakend replacement string

            // TCG -> TG,T,TCAG
                // TCG vs TG -> single base deletion: T-G
                // TCG - T -> two base deletion: T--
                // TCG - TCAG -> single base insertion: TCAG
                // substitution? size must be equal?

            // Find the padded leading and trailing region.
            auto [ref_it, alt_it] = std::ranges::mismatch(reference_segment, alternative);
            size_t const ref_prefix_offset = std::ranges::distance(reference_segment.begin(), ref_it);
            size_t const alt_prefix_offset = std::ranges::distance(alternative.begin(), alt_it);

            auto [ref_it_rev, alt_it_rev] =
                std::ranges::mismatch(reference_segment | seqan3::views::drop(ref_prefix_offset) | std::views::reverse,
                                      alternative | seqan3::views::drop(alt_prefix_offset) | std::views::reverse);
            // base returns the iterator to the actual end.
            assert(ref_it <= ref_it_rev.base());
            assert(alt_it <= alt_it_rev.base());

            size_t delta_position = reference_position() + ref_prefix_offset;

            if (alternative.size() < reference_segment.size()) // Deletion.
            {
                size_t const deletion_size = std::ranges::distance(ref_it, ref_it_rev.base());
                delta_events.emplace_back(delta_position, deletion_type{deletion_size});
            }
            else // Substitution or Insertion.
            {
                size_t const variant_size = std::ranges::distance(alt_it, alt_it_rev.base());
                auto variant = alternative.subspan(alt_prefix_offset, variant_size)
                             | seqan3::views::char_to<seqan3::dna5>
                             | seqan3::views::to<std::vector>;

                if (alternative.size() == reference_segment.size()) // Substitution.
                    delta_events.emplace_back(delta_position, substitution_type{std::move(variant)});
                else // Insertion.
                    delta_events.emplace_back(delta_position, insertion_type{std::move(variant)});
            }
        }

        return {true, delta_events};
    }

    //!\brief Determines the sample count.
    size_t sample_count() const noexcept
    {
        return seqan::length(_record.genotypeInfos);
    }

    //!\brief Determines the number of haplotypes per sample.
    size_t haplotypes_per_sample_count() const
    {
        if (sample_count() == 0 || !genotype_info_given())
            return 0;

        auto first_genotype = to_span(_record.genotypeInfos[0]);
        auto genotype_sentinel = std::ranges::find(first_genotype, ':'); // genotype is always first if present

        if (first_genotype.begin() == genotype_sentinel) // expect a non-empty genotype here.
            throw std::length_error{"Expected genotype information for the first sample but couldn't find any!"};

        // Count number of haplotype delimiters and add one for the last haplotype.
        return std::ranges::count_if(first_genotype.begin(), genotype_sentinel, [] (char const c)
        {
            return c == '|' || c == '/';
        }) + 1;
    }

    //!\brief Checks wether genotype infos are given for the current record.
    bool genotype_info_given() const noexcept
    {
        return !seqan::empty(_record.format) && seqan::startsWith(_record.format, "GT");
    }

    //!\brief Returns the seqan char string as a `std::span<char const>`.
    template <typename seqan_range_t>
    //!\cond
        requires std::same_as<std::ranges::range_value_t<seqan_range_t>, char>
    //!\endcond
    std::span<const char> to_span(seqan_range_t const & range) const noexcept
    {
        if (seqan::empty(range))
            return std::span<const char>{};

        return std::span{std::addressof(range[0]), seqan::length(range)};
    }
};

//!\brief Formatted output operator for the augmented vcf record.
template <typename char_t, typename char_traits_t>
std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                       augmented_vcf_record const & record)
{
    auto const & r = record.seqan_record();
    stream << r.rID << "\t"
           << r.beginPos << "\t"
           << r.id << "\t"
           << r.ref << "\t"
           << r.alt << "\t"
           << r.qual << "\t"
           << r.filter << "\t"
           << r.info << "\t"
           << r.format << "\t";
    for (auto genotype : r.genotypeInfos)
        stream << genotype << "\t";

    return stream;
}

jst_t construct_jst_from_vcf(std::filesystem::path const & reference_file, std::filesystem::path const & vcf_file_path)
{
    // Get the application logger.
    auto & log = get_application_logger();

    // ----------------------------------------------------------------------------
    // Parse the reference file.
    // ----------------------------------------------------------------------------

    seqan3::sequence_file_input ref_file{reference_file};
    auto ref_sequence_it = ref_file.begin();

    seqan3::dna5_vector reference = seqan3::get<seqan3::field::seq>(*ref_sequence_it);

    // ----------------------------------------------------------------------------
    // Parse the vcf file.
    // ----------------------------------------------------------------------------

    seqan::VcfFileIn vcf_file{vcf_file_path.c_str()};

    seqan::VcfHeader vcf_header{};
    seqan::readHeader(vcf_header, vcf_file);

    // TODO: lets assume range based interface.
    // When we construct the file the first record could be already processed and the header as well as some auxiliary
    // infos are already stored and accessible.

    // If the file is empty (no record stored), return an empty JST.
    if (seqan::atEnd(vcf_file))
        return jst_t{std::move(reference)};

    // Read first record:
    augmented_vcf_record record{std::move(vcf_header)};
    seqan::readRecord(record.seqan_record(), vcf_file);
    record.initialise_counts(); // initialise counts with first record.
    // Initial jst with known haplotype count.
    jst_t jst{std::move(reference), record.haplotype_count()};

    log(verbosity_level::verbose, logging_level::info, "Processing records:\n", "-----------------------------");
    // Insert the events generated from the record into the jst.
    auto insert_events_from_record = [&] (auto && record)
    {
        log(verbosity_level::verbose, logging_level::info, record);

        if (auto [is_valid_record, delta_events] = record.generate_delta_events(); is_valid_record)
        {
            for (auto shared_event : delta_events)
                if (!jst.insert(std::move(shared_event)))
                    log(verbosity_level::standard, logging_level::error, "Event was not inserted: ", shared_event);
        }
        else
        {
            log(verbosity_level::standard, logging_level::error, "Skipping invalid record: ", record);
        }
    };

    insert_events_from_record(record);

    while (!seqan::atEnd(vcf_file))
    {
        seqan::readRecord(record.seqan_record(), vcf_file);
        insert_events_from_record(record);
    }
    log(verbosity_level::verbose, logging_level::info, "Done\n", "-----------------------------");

    return jst;
}

}  // namespace jstmap
