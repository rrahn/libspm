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
    // implements: journaled_sequence_tree_node_label
    template <typename sequence_strategy_t, typename coverage_t>
    class node_label
    {
    private:
        sequence_strategy_t _sequence_strategy{};
        coverage_t _coverage{};
    protected:

        // using sequence_t = typename journal_t::mapped_type;

    public:

        node_label() = default;
        template <typename _sequence_strategy_t>
            requires (!std::same_as<std::remove_cvref_t<_sequence_strategy_t>, node_label> &&
                      std::constructible_from<sequence_strategy_t, _sequence_strategy_t>)
        explicit node_label(_sequence_strategy_t && sequence_strategy, coverage_t coverage) :
            _sequence_strategy{(_sequence_strategy_t &&) sequence_strategy},
            _coverage{std::move(coverage)}
        {}

        auto sequence() const noexcept -> decltype(_sequence_strategy.sequence())
        {
            return _sequence_strategy.sequence();
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
        void reset(variant_t const & variant, variant_position_t<variant_t> const size)
        {
            _sequence_strategy.record(variant, size);
            _coverage = libjst::coverage(variant);
        }
    };
}  // namespace libjst
