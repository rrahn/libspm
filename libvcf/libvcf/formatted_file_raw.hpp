// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <libvcf/formatted_file_base.hpp>

namespace libvcf {

// we split the format into a function to extract the record and a tokenisation function with fields.

class formatted_file_raw : protected formatted_file_base {

public:
    template <typename format_t>
    formatted_file_raw(std::filesystem::path && file_path, format_t && format) :
        formatted_file_base{std::move(file_path), std::forward<format_t>(format)} {
        // Now we create a specfic format but how do we get the actual type?
        // can we specify the number of fields?
        // how can we specify this as std::any?
        // is there a way to store the function object type?
    }

    auto read_record() {

    }
};

} // namespace libvcf
