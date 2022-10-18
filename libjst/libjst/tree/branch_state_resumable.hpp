// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides independent node state.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

    template <typename branch_state_t>
    class resumable_branch_state : public branch_state_t
    {
    private:

        using base_t = branch_state_t;

        using position_t = int32_t;

        position_t _head_position{};
        position_t _next_position{};
        position_t _end_position{};

        using sequence_iterator_t = std::ranges::iterator_t<typename base_t::sequence_t const>;
    protected:
        using sequence_t = std::ranges::subrange<sequence_iterator_t, sequence_iterator_t, ranges::subrange_kind::sized>;

    public:

        using base_type = base_t;

        resumable_branch_state() =  default;
        resumable_branch_state(base_t && base) :
            base_t{std::move(base)},
            _end_position{std::ranges::ssize(base_t::sequence())}
        {}

        // we can return the base class.
        base_type const & base() const
        {
            return static_cast<base_t const &>(*this);
        }

        sequence_t sequence() const
        {
            assert(_head_position <= _next_position);
            assert(_next_position <= _end_position);

            auto && seq = base_t::sequence();
            auto it = std::ranges::begin(seq);
            it += _head_position;

            position_t subrange_size = std::min(std::ranges::ssize(seq), _next_position) - _head_position;
            return std::ranges::subrange{it, it + subrange_size, subrange_size}; // borrowed range
        }

        template <covered_sequence_variant variant_t>
        void set_branch(variant_t const & variant, position_t const new_branch_size)
        {
            base_t::set_branch(variant);
            _head_position = libjst::position(variant);
            _next_position = _head_position;
            _end_position = new_branch_size; // next variant we set the state too.
        }

        void set_range(position_t const first, position_t const next)
        {
            _head_position = first;
            _next_position = std::min(next, _end_position);
        }
    };
}  // namespace libjst
