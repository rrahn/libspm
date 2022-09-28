// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the basic adapter for the seqan::Pattern objects.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan/find.h>

#include <libcontrib/seqan/concept.hpp>
#include <libcontrib/seqan/container_adapter.hpp>

#include <libjst/concept.hpp>
#include <libjst/journaled_sequence_tree_forward.hpp>

namespace jst::contrib
{
    template <typename derived_operation_t, typename pattern_t>
    class pattern_operation
    {
    private:

        friend derived_operation_t;

        pattern_t *_pattern{};

        pattern_operation() = default;
        explicit pattern_operation(pattern_t *pattern) noexcept : _pattern(pattern)
        {}
        pattern_operation(pattern_operation const &) = default;
        pattern_operation(pattern_operation &&) = default;
        pattern_operation & operator=(pattern_operation const &) = default;
        pattern_operation & operator=(pattern_operation &&) = default;
        ~pattern_operation() = default;
    public:

        constexpr size_t window_size() const noexcept
        {
            assert(_pattern != nullptr);
            return seqan::length(seqan::needle(*_pattern));
        }

        // we will get a subrange, but is starts from the beginning!
        template <std::ranges::viewable_range haystack_t, typename callback_t>
        constexpr auto operator()(haystack_t &&haystack, callback_t &&callback)
        { // we need to construct this haystack here!
            assert(_pattern != nullptr);
            auto seqan_haystack = make_seqan_container((haystack_t &&)haystack);
            using seqan_haystack_t = decltype(seqan_haystack);

            seqan::Finder<seqan_haystack_t> finder{seqan_haystack};
            jst::contrib::set_up(derived(), finder);

            while (seqan::find(finder, *_pattern))
            {
                callback(finder);
            }

            jst::contrib::tear_down(derived(), finder);
        }

    protected:

    private:
        derived_operation_t &derived() noexcept
        {
            return static_cast<derived_operation_t &>(*this);
        }

        derived_operation_t const &derived() const noexcept
        {
            return static_cast<derived_operation_t const &>(*this);
        }

        constexpr friend size_t tag_invoke(std::tag_t<libjst::window_size>, pattern_operation const & me) noexcept
        {
            return me.window_size();
        }

        constexpr friend bool tag_invoke(std::tag_t<libjst::is_resumable>, any_instance_of_t<pattern_operation> const &)
            noexcept
        {
            return false;
        }

        template <typename cpo_t, typename ...args_t>
        constexpr friend void tag_invoke(cpo_t, pattern_operation const &, args_t &&...) noexcept
        {
            // noop
        }
    };

}  // namespace jst::contrib
