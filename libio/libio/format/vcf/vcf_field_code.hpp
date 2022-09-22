// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides adaption of the field code for the fastq format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libio/file/field_code.hpp>

namespace libio
{
    enum class vcf_meta_field
    {
        // -- 0 -- forbidden
        // required fields
        version = 1,
        // optional fields
        meta,
        sample_names
    };

    // define tags for the different fastq fields.
    enum class vcf_field
    {
        // -- 0 -- forbidden
        // required fields
        chrom = 1,
        pos = 2,
        id = 3,
        ref = 4,
        alt = 5,
        qual = 6,
        filter = 7,
        info = 8,
        // optional fields
        genotype_format = 9,
        genotypes = 10
    };

} // namespace libio
