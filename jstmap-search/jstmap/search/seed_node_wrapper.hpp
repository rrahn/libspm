// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

namespace jstmap
{
    template <typename extender_t>
    class extension_state_manager {
    private:

        using state_t = libjst::matcher_state_t<extender_t>;
        using state_stack_t = std::stack<state_t>;

        extender_t & _extender;
        state_stack_t _states{};

    public:

        constexpr explicit extension_state_manager(extender_t & extender) noexcept :
            _extender{extender},
            _states{}
        {}

        constexpr void notify_push() {
            _states.push(_extender.capture());
        }

        constexpr void notify_pop() {
            assert(!_states.empty());
            _extender.restore(_states.top());
            _states.pop();
        }

        // constexpr void print_state() const noexcept {
        //     seqan3::debug_stream << "[";
        //     std::ranges::for_each(_extender.capture(), [] (uint64_t word) {
        //         seqan3::debug_stream << std::bitset<64>{word} << ", ";
        //     });
        //     seqan3::debug_stream << "]\n";
        // }
    };
}  // namespace jstmap
