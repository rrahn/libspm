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
#include <libio/format/fasta/fasta_token.hpp>
#include <libio/format/fastq/fastq_token.hpp>
#include <libio/record/fasta_record.hpp>
#include <libio/record/fastq_record.hpp>

// Now we just want to tokenize it!

TEST(mini_test, tokenization)
{
    // std::stringstream str_stream{">my 1. id\nAGACTGAGCTACGAGCTAGCGACT\nAGACTGAGCTACGAGCTAGCGACT\nAGACTGAGCTACGAGCTAGCGACT\nAGACTGAGCTACGAGCTAGCGACT\n>my 2. id\nGGTTAAGGTTCCCCAAGGTTAC\n"};

    // libio::stream_token token{str_stream.rdbuf(), seqan3::is_char<'>'>};

    // std::string buffer;
    // buffer.reserve(1000);
    // libio::detokenize_to(token, buffer); // what happens now?

    // std::cout << "First token:\n" << buffer << "\n\n";
    // buffer.clear();
    // libio::stream_token token2{str_stream.rdbuf(), seqan3::is_char<'>'>};
    // libio::detokenize_to(token2, buffer); // what happens now?
    // std::cout << "Second token:\n" << buffer;
}

TEST(mini_test, fasta_token)
{
    std::string_view fa_input = R"fa(>SEQ_ID 1
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
>SEQ_ID 1
GGTTAAGGTTCCCCAAGGTTAC
)fa";
    std::stringstream str_stream{fa_input.data()};
    // adapt the str_stream to inherit from it and move it into it
    // I can move the iterator to the next element
    // I can also get from the tokenized buffer

    libio::fasta_token token{str_stream.rdbuf()};

    libio::fasta_record fa_record{};
    libio::detokenize_to(token, fa_record); // what happens now?

    std::cout << "First token:\nID: " << fa_record.id() << "\nSEQ: " << fa_record.seq() << "\n\n";
    // we may clear the record here?
    libio::fasta_record fa_record2{};
    libio::fasta_token token2{str_stream.rdbuf()};
    libio::detokenize_to(token2, fa_record2); // what happens now?
    std::cout << "Second token:\nID: " << fa_record2.id() << "\nSEQ: " << fa_record2.seq() << "\n\n";
}

TEST(mini_test, fastq_token)
{
    // this is, we set the readbuffer type which provides its own fast iterator scope.

    // @SEQ_ID
    // GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    // +
    // !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
    std::string_view fq_input = R"fq(@SEQ_ID 1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID 2
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65)fq";
    std::stringstream str_stream{fq_input.data()};
    // adapt the str_stream to inherit from it and move it into it
    // I can move the iterator to the next element
    // I can also get from the tokenized buffer

    libio::fastq_token token1{str_stream.rdbuf()};

    libio::fastq_record fq_record1{};
    libio::detokenize_to(token1, fq_record1);

    std::cout << "First token:\nID: " << fq_record1.id() << "\nSEQ: " << fq_record1.seq() << "\nQUAL: " << fq_record1.qual() << "\n\n";
    // we may clear the record here?
    libio::fastq_record fq_record2{};
    libio::fastq_token token2{str_stream.rdbuf()};
    libio::detokenize_to(token2, fq_record2);
    std::cout << "Second token:\nID: " << fq_record2.id() << "\nSEQ: " << fq_record2.seq() << "\nQUAL: " << fq_record2.qual() << "\n\n";
}
