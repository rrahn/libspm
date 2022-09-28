// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides common CPOs and concepts for the jst core.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <libcontrib/std/tag_invoke.hpp>
#include <libcontrib/type_traits.hpp>

namespace libjst
{

    // ----------------------------------------------------------------------------
    // Concept for searcher object.
    // ----------------------------------------------------------------------------

    // window_size(searcher);
    namespace _window_size
    {
        template <typename obj_t>
        concept has_member = requires (obj_t const & o)
        {
            { o.window_size() } -> std::integral;
        };

        inline constexpr struct _cpo
        {
            template <typename searcher_t>
                requires std::tag_invocable<_cpo, searcher_t>
            constexpr auto operator()(searcher_t && searcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t>)
                -> std::tag_invoke_result_t<_cpo, searcher_t>
            {
                return std::tag_invoke(_cpo{}, (searcher_t &&) searcher);
            }

            template <typename searcher_t>
                requires (!std::tag_invocable<_cpo, searcher_t const &> && has_member<searcher_t>)
            constexpr auto operator()(searcher_t const & searcher) const
                noexcept(noexcept(searcher.window_size()))
                -> decltype(searcher.window_size())
            {
                return  searcher.window_size();
            }

            // maybe alternative using free function?
        } window_size{};
    } // namespace _window_size
    using _window_size::window_size;

    // search_operation(searcher);
    namespace _search_operation
    {
        template <typename obj_t>
        concept has_member = requires (obj_t && o)
        {
            ((obj_t &&)o).search_operation(); // check if this makes sense!
        };

        inline constexpr struct _cpo
        {
            template <typename searcher_t>
                requires std::tag_invocable<_cpo, searcher_t>
            constexpr auto operator()(searcher_t && searcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t>)
                -> std::tag_invoke_result_t<_cpo, searcher_t>
            {
                return std::tag_invoke(_cpo{}, (searcher_t &&) searcher);
            }

            // Defaults
            template <typename searcher_t>
                requires (!std::tag_invocable<_cpo, searcher_t>)
            constexpr decltype(auto) operator()(searcher_t && searcher) const
                // noexcept()
            {
                if constexpr (has_member<searcher_t>)
                    return  ((searcher_t &&) searcher).search_operation();
                else
                    return  (searcher_t &&) searcher;
                // maybe alternative using free function?
            }
        } search_operation{};
    } // namespace _search_operation
    using _search_operation::search_operation;

    template <typename searcher_t>
    using search_operation_t = std::remove_cvref_t<std::invoke_result_t<_search_operation::_cpo, searcher_t>>;

    // is_resumable property
    namespace _is_resumable
    {
        inline constexpr struct _cpo
        {
            template <typename searcher_t>
                requires std::tag_invocable<_cpo, searcher_t const &>
            constexpr auto operator()(searcher_t const &searcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, searcher_t const &>)
                -> std::tag_invoke_result_t<_cpo, searcher_t const &>
            {
                return std::tag_invoke(_cpo{}, searcher);
            }
        } is_resumable{};
    } // namespace _is_resumable
    using _is_resumable::is_resumable;

    template <typename searcher_t>
        requires requires { is_resumable(jst::contrib::any_instance_of_v<searcher_t>); }
    inline constexpr bool is_resumable_v = is_resumable(jst::contrib::any_instance_of_v<searcher_t>);

    // ----------------------------------------------------------------------------
    // Concept for jst object.
    // ----------------------------------------------------------------------------

    // search(jst, searcher) -> sender!
    namespace _search
    {
        template <typename obj_t, typename searcher_t>
        concept has_member = requires (obj_t const & o, searcher_t && p)
        {
            o.search((searcher_t &&) p); // check if this makes sense!
        };

        inline constexpr struct _cpo
        {
            // add possible implementation for member function and maybe free function
            template <typename object_t, typename searcher_t>
                requires std::tag_invocable<_cpo, object_t, searcher_t>
            constexpr auto operator()(object_t && object, searcher_t && searcher) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, object_t, searcher_t>)
                -> std::tag_invoke_result_t<_cpo, object_t, searcher_t>
            {
                return std::tag_invoke(_cpo{}, (object_t &&) object, (searcher_t &&) searcher);
            }

            template <typename object_t, typename searcher_t>
                requires (!std::tag_invocable<_cpo, object_t, searcher_t> && has_member<object_t, searcher_t>)
            constexpr auto operator()(object_t &&obj,searcher_t && searcher) const
                noexcept(noexcept(obj.search((searcher_t &&) searcher)))
                -> decltype(obj.search((searcher_t &&) searcher))
            {
                return obj.search((searcher_t &&) searcher);
            }
        } search{};
    } // namespace _search
    using _search::search;


    // ----------------------------------------------------------------------------
    // Concept for sender object.
    // ----------------------------------------------------------------------------
    // TODO: possible interface from future design of unified executor proposal.
    // connect(sender, receiver) -> operation
    namespace _connect
    {
        inline constexpr struct _cpo
        {
            // add possible implementation for member function and maybe free function
            template <typename sender_t, typename receiver_t>
                requires std::tag_invocable<_cpo, sender_t, receiver_t>
            constexpr auto operator()(sender_t && sender, receiver_t && receiver) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, sender_t, receiver_t>)
                -> std::tag_invoke_result_t<_cpo, sender_t, receiver_t>
            {
                return std::tag_invoke(_cpo{}, (sender_t &&) sender, (receiver_t &&) receiver);
            }
        } connect{};
    } // namespace _connect
    using _connect::connect;

    // ----------------------------------------------------------------------------
    // Concept for operation object.
    // ----------------------------------------------------------------------------
    // TODO: possible interface from future design of unified executor proposal.
    // start(operation)
    namespace _start
    {
        inline constexpr struct _cpo
        {
            // add possible implementation for member function and maybe free function
            template <typename operation_t>
                requires std::tag_invocable<_cpo, operation_t>
            constexpr auto operator()(operation_t && operation) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, operation_t>)
                -> std::tag_invoke_result_t<_cpo, operation_t>
            {
                return std::tag_invoke(_cpo{}, (operation_t &&) operation);
            }
        } start{};
    } // namespace _start
    using _start::start;

}  // namespace libjst
