// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <seqan3/std/ranges>
#include <sstream>

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/utility/tag_invoke.hpp>
#include <libio/stream_token.hpp>
#include <libio/format/fasta/fasta_format.hpp>
#include <libio/format/fasta/fasta_token.hpp>
#include <libio/format/fastq/fastq_format.hpp>
#include <libio/format/fastq/fastq_token.hpp>
#include <libio/format/sequence/sequence_format.hpp>
#include <libio/format/sequence/sequence_token.hpp>
#include <libio/format/sequence/sequence_record.hpp>
#include <libio/record/fasta_record.hpp>
#include <libio/record/fastq_record.hpp>

// Now we just want to tokenize it!

inline constexpr std::string_view fa_input = R"fa(>SEQ_ID 1
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
>SEQ_ID 1
GGTTAAGGTTCCCCAAGGTTAC
)fa";

inline constexpr std::string_view fq_input = R"fq(@SEQ_ID 1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID 2
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
)fq";

TEST(mini_test, fasta_token)
{
    std::stringstream str_stream{fa_input.data()};
    libio::fasta_token token{str_stream};

    libio::fasta_record fa_record{};
    libio::detokenize_to(token, fa_record); // what happens now?

    std::cout << "First token:\nID: " << fa_record.id() << "\nSEQ: " << fa_record.seq() << "\n\n";
    // we may clear the record here?
    libio::fasta_record fa_record2{};
    libio::fasta_token token2{str_stream};
    libio::detokenize_to(token2, fa_record2); // what happens now?
    std::cout << "Second token:\nID: " << fa_record2.id() << "\nSEQ: " << fa_record2.seq() << "\n\n";
}

TEST(mini_test, fastq_token)
{
    std::stringstream str_stream{fq_input.data()};
    libio::fastq_token token1{str_stream};

    libio::fastq_record fq_record1{};
    libio::detokenize_to(token1, fq_record1);

    std::cout << "First token:\nID: " << fq_record1.id() << "\nSEQ: " << fq_record1.seq() << "\nQUAL: " << fq_record1.qual() << "\n\n";
    // we may clear the record here?
    libio::fastq_record fq_record2{};
    libio::fastq_token token2{str_stream};
    libio::detokenize_to(token2, fq_record2);
    std::cout << "Second token:\nID: " << fq_record2.id() << "\nSEQ: " << fq_record2.seq() << "\nQUAL: " << fq_record2.qual() << "\n\n";
}

TEST(mini_test, sequence_token_as_fasta)
{
    // right! they are more difficult to realize!
    using sequence_token_t = libio::sequence_token<libio::fasta_token<std::stringstream>, libio::fastq_token<std::stringstream>>;
    std::stringstream str_stream{fa_input.data()};
    sequence_token_t token1{libio::fasta_token{str_stream}};

    libio::sequence_record seq_record1{};
    libio::detokenize_to(token1, seq_record1);

    std::cout << "First token:\nID: " << seq_record1.id() << "\nSEQ: " << seq_record1.seq() << "\nQUAL: " << seq_record1.qual() << "\n\n";
    // we may clear the record here?
    libio::sequence_record seq_record2{};
    sequence_token_t token2{libio::fasta_token{str_stream}};
    libio::detokenize_to(token2, seq_record2);
    std::cout << "Second token:\nID: " << seq_record2.id() << "\nSEQ: " << seq_record2.seq() << "\nQUAL: " << seq_record2.qual() << "\n\n";
}

TEST(mini_test, sequence_token_as_fastq)
{
    // right! they are more difficult to realize!
    using sequence_token_t = libio::sequence_token<libio::fasta_token<std::stringstream>, libio::fastq_token<std::stringstream>>;
    std::stringstream str_stream{fq_input.data()};
    sequence_token_t token1{libio::fastq_token{str_stream}};

    libio::sequence_record seq_record1{};
    libio::detokenize_to(token1, seq_record1);

    std::cout << "First token:\nID: " << seq_record1.id() << "\nSEQ: " << seq_record1.seq() << "\nQUAL: " << seq_record1.qual() << "\n\n";
    // we may clear the record here?
    libio::sequence_record seq_record2{};
    sequence_token_t token2{libio::fastq_token{str_stream}};
    libio::detokenize_to(token2, seq_record2);
    std::cout << "Second token:\nID: " << seq_record2.id() << "\nSEQ: " << seq_record2.seq() << "\nQUAL: " << seq_record2.qual() << "\n\n";
}

TEST(mini_test, fasta_format)
{
    std::stringstream str_stream{fa_input.data()};
    libio::fasta_format fmt{};
    auto token = libio::format_token(fmt, str_stream);

    libio::fasta_record fa_record{};
    libio::detokenize_to(token, fa_record); // what happens now?

    std::cout << "First token:\nID: " << fa_record.id() << "\nSEQ: " << fa_record.seq() << "\n\n";
    // we may clear the record here?
    libio::fasta_record fa_record2{};
    auto token2 = libio::format_token(fmt, str_stream);
    libio::detokenize_to(token2, fa_record2); // what happens now?
    std::cout << "Second token:\nID: " << fa_record2.id() << "\nSEQ: " << fa_record2.seq() << "\n\n";
}
