// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides restorable pattern wrapper.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace spm
{
    template <typename tag_t>
    struct Restorable : public tag_t {

        // constexpr operator tag_t() const noexcept {
        //     return static_cast<tag_t>(*this);
        // }
    };

}  // namespace spm
