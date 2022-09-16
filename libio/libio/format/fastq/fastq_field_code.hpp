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
    // define tags for the different fastq fields.
    enum class fastq_field
    {
        // -- 0 -- forbidden
        id = 1,
        seq,
        qual
    };

    namespace detail {
        // Do we need this?
        // class fastq_field_category final : public field_category
        // {
        // public:

        //     virtual std::string_view name() const noexcept override {
        //         return "fastq";
        //     }

        // private:
        //     // we can overload the equivalence value.
        //     // in this case we can convert what?
        //     // how do we prevent specific overload issues?
        // };

    } // namespace detail

    // now global variable only visible in current namespace?
    // inline constexpr detail::fastq_field_category fastq_field_category{};

    // to translate this into a field code?
    // overload which needs to be found by dependent lookup
    // but only if field code fulfils some requirements!
    // inline constexpr detail::fastq_field_category const & select_field_category(fastq_field f) noexcept
    // {
    //     return fastq_field_category;
    // }
} // namespace libio
