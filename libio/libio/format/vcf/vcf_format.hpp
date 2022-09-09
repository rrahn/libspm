// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

namespace libio
{

    class vcf_format_input
    {

        // tag invoke read_header() // reads tokenised into a raw_record?
        // gets a field delimitter?
        void read_header()
        {
        }
    };

    class vcf_format_output
    {
    };

    class vcf_format
    {
    };
} // namespace libio
