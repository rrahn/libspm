// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides statistics for jst.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

namespace libjst
{

    struct tree_stats {
        size_t node_count{};
        size_t subtree_count{};
        size_t leaf_count{};
        size_t symbol_count{};
        size_t max_subtree_depth{};
        std::vector<size_t> subtree_depths{};

        constexpr void notify_push() {
        }

        constexpr void notify_pop() {
            ++leaf_count;
        }
    };

    struct node_properties {
        std::ptrdiff_t subtree_depth{};
    };

    template <typename base_tree_t>
        // requires covered tree
    class stats_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;
        class label_impl;

        base_tree_t _wrappee{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr stats_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<std::remove_cvref_t<wrapped_tree_t>, stats_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr stats_tree_impl(wrapped_tree_t && wrappee) noexcept : _wrappee{(wrapped_tree_t &&)wrappee}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee)};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }
   };

    template <typename base_tree_t>
    class stats_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend stats_tree_impl;

        node_properties _properties{.subtree_depth = 0};

        explicit constexpr node_impl(base_node_type && base_node) noexcept : base_node_type{std::move(base_node)}
        {}

        explicit constexpr node_impl(base_node_type && base_node, node_properties properties) noexcept :
            base_node_type{std::move(base_node)},
            _properties{std::move(properties)}
        {}

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit<false>(base_node_type::next_ref());
        }

        constexpr label_impl operator*() const noexcept {
            return label_impl{*static_cast<base_node_type const &>(*this), _properties};
        }

    private:

        template <bool is_alt, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                node_impl new_child{std::move(*maybe_child), _properties};
                if (base_node_type::on_alternate_path()) {
                    ++new_child._properties.subtree_depth;
                } else {
                    if constexpr (is_alt) {
                        ++new_child._properties.subtree_depth;
                    }
                }
                return new_child;
            }

            return std::nullopt;
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class stats_tree_impl<base_tree_t>::sink_impl {
    private:
        friend stats_tree_impl;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    public:
        sink_impl() = default;
    };

    template <typename base_tree_t>
    class stats_tree_impl<base_tree_t>::label_impl : public libjst::tree_label_t<base_tree_t> {
    private:
        friend stats_tree_impl;

        using base_label_type = libjst::tree_label_t<base_tree_t>;

        node_properties const & _properties;

        constexpr explicit label_impl(base_label_type base_label, node_properties const & properties) :
            base_label_type{std::move(base_label)},
            _properties{properties}
        {}

    public:
        label_impl() = delete;

        constexpr size_t subtree_depth() const noexcept {
            return _properties.subtree_depth;
        }

        constexpr bool is_subtree_root() const noexcept {
            return _properties.subtree_depth == 1;
        }
    };

    namespace _tree_adaptor {
        inline constexpr struct _stats
        {
            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, [[maybe_unused]] args_t &&...args) const
            {
                tree_stats stats{.node_count = 0,
                                 .subtree_count = 0,
                                 .leaf_count = 1,
                                 .symbol_count = 0};
                stats_tree_impl<tree_t> stats_tree{(tree_t &&)tree};

                tree_traverser_base path{stats_tree};
                path.subscribe(stats);
                for (auto it = path.begin(); it != path.end(); ++it) {
                    auto node_properties = *it;
                    ++stats.node_count;
                    stats.symbol_count += std::ranges::size(node_properties.sequence());
                    if (node_properties.is_subtree_root()) {
                        ++stats.subtree_count;
                        stats.subtree_depths.push_back(1);
                    } else if (node_properties.subtree_depth() > 0) {
                        stats.subtree_depths.back() = std::max(stats.subtree_depths.back(),
                                                               node_properties.subtree_depth());
                    }
                }

                std::ranges::for_each(stats.subtree_depths, [&] (size_t const & depth) {
                    stats.max_subtree_depth = std::max(stats.max_subtree_depth, depth);
                });
                return stats;
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t && ...args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_stats, args_t...>
            {
                return jst::contrib::make_closure(_stats{}, (args_t &&)args...);
            }
        } stats{};
    } // namespace _tree_adaptor

    using _tree_adaptor::stats;
}  // namespace libjst
