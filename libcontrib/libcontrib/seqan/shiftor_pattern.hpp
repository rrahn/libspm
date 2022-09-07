// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an adapter to make the horspool algorithm work with the JST.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan/sequence.h>
#include <seqan/find.h>

#include <libcontrib/seqan/concept.hpp>
#include <libcontrib/seqan/pattern_operation.hpp>
#include <libcontrib/std/tag_invoke.hpp>
#include <libcontrib/type_traits.hpp>

namespace jst::contrib
{
    template <typename needle_t>
    class shiftor_pattern
    {
    private:

        using pattern_t = seqan::Pattern<needle_t, seqan::ShiftOr>;
        struct operation;

        pattern_t _pattern{};

    public:

        shiftor_pattern() = delete;
        template <typename _needle_t>
            requires std::constructible_from<pattern_t, _needle_t>
        explicit shiftor_pattern(_needle_t && needle) : _pattern{(_needle_t &&) (needle)}
        {
        }

        constexpr operation search_operation() noexcept
        {
            return operation{_pattern};
        }
    };

    template <typename needle_t>
    shiftor_pattern(needle_t &&) -> shiftor_pattern<std::remove_reference_t<needle_t>>;

    template <typename needle_t>
    class shiftor_pattern<needle_t>::operation : public pattern_operation<operation, pattern_t> // run as CRTP?
    {
    private:
        // extract special state from the actual pattern?
        using base_t = pattern_operation<operation, pattern_t>;
        using match_t = std::remove_cvref_t<decltype(std::declval<pattern_t>().prefSufMatch)>;

        friend base_t;

        match_t _match{};
        size_t _cached_position{};
        bool _find_first{true};
    public:

        operation() = default;
        explicit operation(pattern_t & pattern) : base_t{std::addressof(pattern)}
        {}

    private:

        template <typename finder_t>
        constexpr friend void tag_invoke(std::tag_t<jst::contrib::set_up>, operation & me, finder_t & finder)
            noexcept
        {
            if (!me._find_first) // we overwrite the behaviour of the seqan::find method implemented for shiftor.
            {
                using std::swap;
                swap(me.pattern().prefSufMatch, me._match); // swap with last state.
                seqan::_setFinderLength(finder, static_cast<base_t const &>(me).window_size());
                seqan::_setFinderEnd(finder, me._cached_position);
                seqan::setPosition(finder, seqan::beginPosition(finder));
                seqan::_finderSetNonEmpty(finder);
            }
            me._find_first = false;
        }

        template <typename finder_t>
        constexpr friend void tag_invoke(std::tag_t<jst::contrib::tear_down>, operation & me, finder_t & finder)
            noexcept
        {
            using std::swap;
            swap(me._match, me.pattern().prefSufMatch); // copy state?
            me._cached_position = seqan::length(seqan::haystack(finder));
        }

        constexpr friend bool tag_invoke(std::tag_t<libjst::is_resumable>, any_instance_of_t<operation> const &) noexcept
        {
            return true;
        }

        // member type trait
        template <typename cpo_t, typename operation_t, typename ...args_t>
            requires std::same_as<std::remove_cvref_t<operation_t>, operation>
        constexpr friend auto tag_invoke(cpo_t cpo, operation_t && me, args_t &&...args)
            noexcept(std::is_nothrow_invocable_v<cpo_t, member_type_t<operation_t, base_t>, args_t...>)
            -> std::invoke_result_t<cpo_t, member_type_t<operation_t, base_t>, args_t...>
        {
            using fwd_base_t = member_type_t<operation_t, base_t>;
            return cpo(static_cast<fwd_base_t>(me), (args_t &&) args...);
        }
    };

}  // namespace jst::contrib
