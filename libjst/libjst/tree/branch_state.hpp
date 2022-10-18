// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides value for jst node.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
    template <typename journal_t, typename coverage_t>
    class jst_branch_state
    {
    private:
        journal_t _journal{};
        coverage_t _coverage{};
        int32_t _offset{};
    protected:

        using sequence_t = typename journal_t::mapped_type;

    public:

        jst_branch_state() = default;
        template <typename base_sequence_t>
            requires (!std::same_as<base_sequence_t, jst_branch_state> &&
                      std::constructible_from<journal_t, base_sequence_t const &>)
        explicit jst_branch_state(base_sequence_t const & base, coverage_t coverage) :
            _journal{base},
            _coverage{std::move(coverage)}
        {}

        auto sequence() const noexcept
        {
            return _journal.sequence();
        }

        coverage_t const & coverage() const noexcept
        {
            return _coverage;
        }

        coverage_t & coverage() noexcept
        {
            return _coverage;
        }

        constexpr bool has_value() const noexcept
        {
            return _coverage.any();
        }

        constexpr operator bool() const noexcept
        {
            return has_value();
        }

        template <covered_sequence_variant variant_t>
        void set_branch(variant_t const & variant)
        {
            record_sequence_variant(variant);
            _offset += std::ranges::ssize(libjst::insertion(variant)) - libjst::deletion(variant);
            _coverage &= libjst::coverage(variant);
        }

        void unset(coverage_t const & coverage)
        {
            _coverage.and_not(coverage);
        }
    private:

        template <typename sequence_variant_t>
        void record_sequence_variant(sequence_variant_t const &variant)
        {
            auto position = libjst::position(variant) + _offset;
            if (libjst::is_insertion(variant)) {
                _journal.record_insertion(position, libjst::insertion(variant));
            } else if (libjst::is_deletion(variant)) {
                _journal.record_deletion(position, libjst::deletion(variant));
            } else {
                assert(libjst::is_replacement(variant));
                _journal.record_substitution(position, libjst::insertion(variant));
            }
        }

    };
}  // namespace libjst
