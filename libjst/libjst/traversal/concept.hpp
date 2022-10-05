// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides CPOs for the search interface.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Operation CPOs
    // ----------------------------------------------------------------------------

    // start()
    namespace _start {
        inline constexpr struct _cpo  {

            template <typename operation_t>
                requires std::tag_invocable<_cpo, operation_t &>
            constexpr void operator()(operation_t &operation) const noexcept
            {
                std::tag_invoke(_cpo{}, operation);
            }
        } start;
    }
    using _start::start;

    // ----------------------------------------------------------------------------
    // Searcher CPOs
    // ----------------------------------------------------------------------------

    // connect(publisher)
    namespace _connect {
        inline constexpr struct _cpo  {

            template <typename searcher_t, typename publisher_t>
                requires std::tag_invocable<_cpo, searcher_t, publisher_t>
            constexpr auto operator()(searcher_t &&searcher, publisher_t &&publisher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t, publisher_t>)
                -> std::tag_invoke_result_t<_cpo, searcher_t, publisher_t> // must be a operation
            {
                return std::tag_invoke(_cpo{}, (searcher_t &&)searcher, (publisher_t &&)publisher);
            }
        } connect;
    }
    using _connect::connect;

    template <typename searcher_t, typename publisher_t>
    using operation_t = std::remove_cvref_t<std::tag_invoke_result_t<_connect::_cpo, searcher_t, publisher_t>>;

    // search_state() -> getter
    // search_state(state) -> setter
    namespace _search_state {
        inline constexpr struct _cpo  {
            // get state
            template <typename searcher_t>
                requires std::tag_invocable<_cpo, searcher_t>
            constexpr auto operator()(searcher_t &&searcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t>)
                -> std::tag_invoke_result_t<_cpo, searcher_t>
            {
                return std::tag_invoke(_cpo{}, (searcher_t &&) searcher);
            }

            // set state
            template <typename searcher_t, typename state_t>
                requires std::tag_invocable<_cpo, searcher_t, state_t>
            constexpr void operator()(searcher_t &&searcher, state_t && state) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t, state_t>)
            {
                return std::tag_invoke(_cpo{}, (searcher_t &&) searcher, (state_t &&)state);
            }
        } search_state;
    }
    using _search_state::search_state;

    template <typename searcher_t>
    using search_state_t = std::remove_cvref_t<std::tag_invoke_result_t<_search_state::_cpo, searcher_t>>;

    // ----------------------------------------------------------------------------
    // Publisher CPOs
    // ----------------------------------------------------------------------------

    // set_next(search_result) // notify next result
    namespace _set_next {
        inline constexpr struct _cpo  {

            template <typename publisher_t, typename search_result_t>
                requires std::tag_invocable<_cpo, publisher_t, search_result_t>
            constexpr void operator()(publisher_t &&publisher, search_result_t &&result) const noexcept
            {
                std::tag_invoke(_cpo{}, (publisher_t &&)publisher, (search_result_t &&)result);
            }
        } set_next;
    }
    using _set_next::set_next;

    // set_value() // notify finished
    namespace _set_value {
        inline constexpr struct _cpo  {

            template <typename publisher_t>
                requires std::tag_invocable<_cpo, publisher_t>
            constexpr void operator()(publisher_t &&publisher) const noexcept
            {
                std::tag_invoke(_cpo{}, (publisher_t &&)publisher);
            }
        } set_value;
    }
    using _set_value::set_value;

    // set_done() // notify abort
    namespace _set_done {
        inline constexpr struct _cpo  {

            template <typename publisher_t>
                requires std::tag_invocable<_cpo, publisher_t>
            constexpr void operator()(publisher_t &&publisher) const noexcept
            {
                std::tag_invoke(_cpo{}, (publisher_t &&)publisher);
            }
        } set_done;
    }
    using _set_done::set_done;

    // set_error() // notify exception
    namespace _set_error {
        inline constexpr struct _cpo  {

            template <typename publisher_t, typename error_t>
                requires std::tag_invocable<_cpo, publisher_t, error_t>
            constexpr void operator()(publisher_t &&publisher, error_t && error) const noexcept
            {
                std::tag_invoke(_cpo{}, (publisher_t &&)publisher, (error_t &&)error);
            }
        } set_error;
    }
    using _set_error::set_error;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename operation_t>
    concept search_operation = std::destructible<std::remove_reference_t<operation_t>> && requires
    (std::remove_reference_t<operation_t> & op)
    {
        { libjst::start(op) };
        requires noexcept(libjst::start(op)); // must not throw
    };

    template <typename publisher_t, typename search_result_t>
    concept search_result_publisher = requires
    (publisher_t && p, search_result_t && r)
    {
        { libjst::set_next((publisher_t &&)p, (search_result_t &&)r) };
        requires noexcept(libjst::set_next((publisher_t &&)p, (search_result_t &&)r)); // must not throw

        { libjst::set_value((publisher_t &&)p) };
        requires noexcept(libjst::set_value((publisher_t &&)p)); // must not throw
    };

    template <typename publisher_t>
    concept search_done_publisher = requires
    (publisher_t && p)
    {
        { libjst::set_done((publisher_t &&)p) };
        requires noexcept(libjst::set_done((publisher_t &&)p)); // must not throw
    };

    template <typename publisher_t, typename error_t>
    concept search_error_publisher = requires
    (publisher_t && p, error_t && e)
    {
        { libjst::set_error((publisher_t &&)p, (error_t &&)e) };
        requires noexcept(libjst::set_error((publisher_t &&)p, (error_t &&)e)); // must not throw
    };

    template <typename searcher_t, typename publisher_t>
    concept sender_to = requires
    (searcher_t && s, publisher_t && p)
    {
        { libjst::connect((searcher_t &&)s, (publisher_t &&)p) } -> search_operation;
    };
}  // namespace libjst
