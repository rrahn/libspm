// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides cached branch state.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
    template <typename branch_state_t>
    class chached_branch_state : public branch_state_t
    {
    private:
        bool _is_cached{false};
        bool _cached_result{false};
    public:

        chached_branch_state() = default;
        explicit chached_branch_state(branch_state_t && base) : branch_state_t{std::move(base)}
        {}

        void reset_coverage(coverage_t const & coverage)
        {
            branch_state_t::reset_coverage(coverage);
            _is_cached = false;
        }

        template <covered_sequence_variant variant_t>
        void set_branch(variant_t const & variant)
        {
            branch_state_t::set_branch(variant);
            _is_cached = false;
        }

        constexpr bool has_value() const noexcept
        {
            if (!_is_cached) {
                _cached_result = branch_state_t::has_value();
                _is_cached = true;
            }

            return _cached_result;
        }
    };

}  // namespace libjst
