// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

#include <seqan3/utility/views/type_reduce.hpp>

#include <libjst/container/concept_jst.hpp>
#include <libjst/container/concept_serialiser.hpp>
#include <libjst/container/journaled_path.hpp>
#include <libjst/container/serialiser_delegate.hpp>
#include <libjst/variant/variant_store_sorted.hpp>
// #include <libjst/traversal/source_branch_node.hpp>
// #include <libjst/traversal/node_label.hpp>

namespace libjst
{

    namespace detail
    {
        template <typename jst_t>
        using maybe_reference_wrapper_t = std::conditional_t<std::is_lvalue_reference_v<jst_t>,
                                                             std::reference_wrapper<std::remove_reference_t<jst_t>>,
                                                             jst_t>;
        template <typename jst_t>
        class jst_box : public jst::contrib::copyable_box<maybe_reference_wrapper_t<jst_t>>
        {
            using base_t = jst::contrib::copyable_box<maybe_reference_wrapper_t<jst_t>>;
            // TODO: add serialiser functionality
        public:

            using base_t::base_t;

            constexpr jst_t & operator*() noexcept(noexcept(*std::declval<base_t &>()))
            {
                return *static_cast<base_t &>(*this);
            }

            constexpr jst_t const & operator*() const noexcept(noexcept(*std::declval<base_t const &>()))
            {
                return *static_cast<base_t const &>(*this);
            }
        };

        namespace _forward_jst
        {
            inline constexpr struct _fn
            {
                template <typename node_factory_t, std::ranges::view source_t, std::ranges::view store_t>
                constexpr auto operator()(node_factory_t && make_node, source_t &&source, store_t &&store) const
                {
                    return make_node((source_t &&)source, (store_t &&)store);
                }

                template <std::ranges::view source_t, std::ranges::view store_t>
                constexpr auto operator()(source_t &&source, store_t &&store) const
                    noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _fn, source_t, store_t>)
                    -> jst::contrib::closure_result_t<_fn, source_t, store_t>
                {
                    return jst::contrib::make_closure(_fn{}, (source_t &&)source, (store_t &&)store);
                }
            } root{};

        } // _forward_jst

    } // namespace detail

    template <typename base_set_t>
    class jst_forward : private traversable_jst_base
    {
    private:
        using store_t = variant_store_t<base_set_t>;
        using sorted_store_t = variant_store_sorted<store_t const>;

        class ordered_variant;

        // we need to put it into a moveable-box!
        // then we do not need to store a reference object here.
        [[no_unique_address]] detail::jst_box<base_set_t> _wrappee{};
        sorted_store_t _store{};

    public:

        jst_forward() = default;
        jst_forward(jst_forward const &) = default;
        jst_forward(jst_forward &&) = default;
        jst_forward & operator=(jst_forward const &) = default;
        jst_forward & operator=(jst_forward &&) = default;

        template <typename _base_set_t>
            requires (!std::same_as<std::remove_cvref_t<_base_set_t>, jst_forward> &&
                       std::same_as<std::remove_cvref_t<_base_set_t>, std::remove_cvref_t<base_set_t>>)
        explicit jst_forward(_base_set_t &&base_set) :
            _wrappee{(_base_set_t &&) base_set},
            _store{libjst::variant_store(*_wrappee)}
        {}

    private:
        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load>, jst_forward & me, archive_t & archive)
        {
            libjst::load_extern(archive, *me._wrappee);
            archive(me._store);
        }

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save>, jst_forward const & me, archive_t & archive)
        {
            libjst::save_extern(archive, *me._wrappee);
            archive(me._store);
        }

        template <typename cpo_t>
            requires std::invocable<cpo_t, base_set_t const &>
        constexpr friend auto tag_invoke(cpo_t cpo, jst_forward const &me)
            noexcept(std::is_nothrow_invocable_v<cpo_t, base_set_t const &>)
            -> std::invoke_result_t<cpo_t, base_set_t const &>
        {
            return cpo(*me._wrappee);
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::variant_store>, jst_forward const &me) noexcept
            -> sorted_store_t const &
        {
            return me._store;
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::path>, jst_forward const &me) noexcept
        {
            return journaled_path{libjst::base_sequence(me) | seqan3::views::type_reduce, me._store};
        }

    };

    template <typename base_set_t>
    jst_forward(base_set_t &&) -> jst_forward<base_set_t>;

    namespace _forward_jst
    {
        inline constexpr struct _fn
        {
            template <journaled_sequence_tree_c base_set_t>
            constexpr auto operator()(base_set_t &&wrappee) const
                -> jst_forward<base_set_t>
            {
                return jst_forward<base_set_t>{(base_set_t &&)wrappee};
            }

            template <typename ...args_t>
                requires (sizeof...(args_t) == 0)
            constexpr auto operator()(args_t &&...) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _fn, args_t...>)
                -> jst::contrib::closure_result_t<_fn, args_t...>
            {
                return jst::contrib::make_closure(_fn{});
            }

        } forward_jst;
    } // namespace _forward_jst

    using _forward_jst::forward_jst;
}  // namespace libjst
