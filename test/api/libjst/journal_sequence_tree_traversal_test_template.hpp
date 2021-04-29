// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <ranges>
#include <set>
#include <span>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/journaled_sequence_tree.hpp>

#include "test_utility.hpp" // make_gaped

using namespace std::literals;

namespace libjst::test
{

using alphabet_t = char;
using shared_event_t = libjst::detail::delta_event_shared<alphabet_t>;
using delta_event_t = typename shared_event_t::delta_event_type;
using substitution_t = typename shared_event_t::substitution_type;
using insertion_t = typename shared_event_t::insertion_type;
using deletion_t = typename shared_event_t::deletion_type;
using coverage_t = typename shared_event_t::coverage_type;
using jst_events_t = std::vector<shared_event_t>;

struct traversal_fixture
{
    std::string reference{};
    size_t sequence_count{};
    jst_events_t events{};
    uint32_t context_size{};
    uint32_t bin_count{};

    template <typename char_t>
    friend seqan3::debug_stream_type<char_t> & operator<<(seqan3::debug_stream_type<char_t> & stream,
                                                          traversal_fixture const & fixture)
    {
        return stream << "["
                      << "reference: " << fixture.reference << ", "
                      << "sequence_count: " << fixture.sequence_count << ", "
                      << "events: " << fixture.events << ", "
                      << "context_size: " << fixture.context_size << ", "
                      << "bin_count: " << fixture.bin_count << "]";
    }
};

struct traversal_fixture_base : public ::testing::TestWithParam<traversal_fixture>,
                                public libjst::test::jst_context_map_fixture
{
    using aligned_sequence_t = std::vector<seqan3::gapped<alphabet_t>>;
    using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;

    std::vector<std::string> sequences{}; // The generated sequences from the delta events
    std::vector<alignment_t> alignments{}; // The alignments against the reference sequence generated from the delta events.

    void SetUp() override
    {
        generate_alignments();
        generate_context_map(GetParam().context_size, sequences);
    }

    auto construct_jst() const
    {
        libjst::journaled_sequence_tree<std::string> jst{std::string{GetParam().reference}};

        std::ranges::for_each(alignments, [&] (alignment_t const alignment)
        {
            jst.add(alignment);
        });

        return jst;
    }

private:

    void generate_alignments()
    {
        // We generate all sequences here from the reference and the events.
        // Then we create a map with position and context.
        sequences.resize(GetParam().sequence_count);
        alignments.resize(GetParam().sequence_count);

        // We construct the sequences from the reference and the given delta map.
        for (unsigned i = 0; i < GetParam().sequence_count; ++i)
        {
            std::string first_sequence = GetParam().reference;
            std::string second_sequence = first_sequence;

            int32_t virtual_offset = 0;
            std::ranges::for_each(GetParam().events, [&] (shared_event_t const & event)
            {
                EXPECT_EQ(event.coverage().size(), GetParam().sequence_count);

                // Continue only if the coverage is true for this sequence.
                if (!event.coverage()[i])
                    return;

                std::visit([&] (auto const & event_kind)
                {
                    size_t const event_position = event.position() + virtual_offset;
                    size_t const insertion_size = event.insertion_size();
                    size_t const deletion_size = event.deletion_size();

                    EXPECT_LE(event_position, first_sequence.size());
                    EXPECT_LE(event_position, second_sequence.size());

                    seqan3::detail::multi_invocable
                    {
                        [&] <typename alphabet_t> (libjst::detail::delta_kind_substitution<alphabet_t> const & s)
                        {
                            // aaaaaaaaa
                            // aaaabbbaa
                            second_sequence.replace(event_position, insertion_size, s.value().data(), insertion_size);
                        },
                        [&] <typename alphabet_t> (libjst::detail::delta_kind_insertion<alphabet_t> const & i)
                        {
                            // aaaa--aaaaa
                            // aaaabbaaaaa
                            first_sequence.insert(event_position, insertion_size, '-');
                            second_sequence.insert(second_sequence.begin() + event_position,
                                                   i.value().begin(), i.value().end());
                            virtual_offset += insertion_size;
                        },
                        [&] (libjst::detail::delta_kind_deletion const &)
                        {
                            // aaaaaaaaaaaa
                            // aaaaa----aaa
                            second_sequence.replace(event_position, deletion_size, deletion_size, '-');
                        },
                        [] (...)
                        {
                            FAIL();
                        }
                    } (event_kind);
                }, event.delta_variant());
            });

            // Store the generated alignment.
            alignments[i] = std::pair{libjst::test::make_gapped(first_sequence),
                                      libjst::test::make_gapped(second_sequence)};

            std::erase(second_sequence, '-');
            sequences[i] = second_sequence;
        }
    }
};

} // namespace libjst::test
