// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also availabel at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/extendable_tree.hpp>
#include <libjst/sequence_tree/journaled_sequence_label.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>

namespace libjst
{

    enum struct sequence_label_kind {
        root_path,
        node_only
    };

    template <typename derived_t, typename extended_node_t, typename label_kind_t>
    class node_label_extension {
    private:

        friend derived_t;

        using sequence_type = typename rcs_node_traits<extended_node_t>::sequence_type;
        using breakpoint_type = typename rcs_node_traits<extended_node_t>::breakpoint_type::value_type;

        using label_strategy_type = journaled_sequence_label<breakpoint_type, sequence_type const &>;

        template <typename base_label_t>
        class label_impl;

        label_strategy_type _label_strategy{}; // the label information

    private:

        node_label_extension() = default;

        constexpr void initialise() {
            _label_strategy = label_strategy_type(as_derived(*this).rcs_store().source());
            _label_strategy.reset_positions(as_derived(*this).left_breakpoint(), as_derived(*this).right_breakpoint());
        };

        constexpr node_label_extension notify(extended_node_t const & child_node) const {
            node_label_extension child_extension{*this}; // copy current extension into child.
            if (as_derived(child_node).is_alt_node()) {
                child_extension._label_strategy.record_variant(as_derived(child_node).left_variant());
            } else {
                // now we could be in a multibranch state!
                // thus we need to set the correct left position when we go into state A!

                auto right_bp = as_derived(child_node).right_breakpoint().value();
                auto left_bp = as_derived(child_node).left_breakpoint().value();
                child_extension._label_strategy.reset_positions(std::min(left_bp, right_bp), right_bp);
            }
            return child_extension;
        }

        static constexpr derived_t & as_derived(node_label_extension & me) noexcept {
            return static_cast<derived_t &>(me);
        }

        static constexpr derived_t & as_derived(extended_node_t & me) noexcept {
            return static_cast<derived_t &>(me);
        }

        static constexpr derived_t const & as_derived(node_label_extension const & me) noexcept {
            return static_cast<derived_t const &>(me);
        }

        static constexpr derived_t const & as_derived(extended_node_t const & me) noexcept {
            return static_cast<derived_t const &>(me);
        }

    public:

        template <typename base_label_t>
        constexpr auto label(base_label_t && base_label) const noexcept {
            return label_impl<std::remove_cvref_t<base_label_t>>{(base_label_t &&) base_label,
                                                                 std::addressof(_label_strategy)};
        }
    };

    template <typename derived_t, typename extended_node_t, typename label_kind_t>
    template <typename base_label_t>
    class node_label_extension<derived_t, extended_node_t, label_kind_t>::label_impl : public base_label_t {
    private:

        [[no_unique_address]] label_strategy_type const * _label_strategy{};

        friend node_label_extension;

        constexpr explicit label_impl(base_label_t base_label, label_strategy_type const * strategy) :
            base_label_t{std::move(base_label)},
            _label_strategy{strategy}
        {}
    public:

        using size_type = typename label_strategy_type::size_type;

        static constexpr size_type npos{label_strategy_type::npos};

        label_impl() = default;

        constexpr auto sequence(size_type const first = 0, size_type const last = label_strategy_type::npos) const noexcept {
            assert(_label_strategy != nullptr);
            assert(first <= last);
            return _label_strategy->sequence(first, last);
        }
    };

    template <typename base_tree_t, sequence_label_kind label_kind>
    using labelled_tree_impl = extendable_tree<base_tree_t,
                                               node_label_extension,
                                               std::integral_constant<sequence_label_kind, label_kind>>;

    namespace _tree_adaptor {
        template <sequence_label_kind label_kind>
        struct _cpo {

            template <typename tree_t, typename ...args_t>
            constexpr auto operator()(tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            labelled_tree_impl<std::remove_reference_t<tree_t>, label_kind, args_t...>>)
                -> labelled_tree_impl<std::remove_reference_t<tree_t>, label_kind, args_t...>
            {
                using adapted_tree_t = labelled_tree_impl<std::remove_reference_t<tree_t>, label_kind, args_t...>;
                return adapted_tree_t{(tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_cpo, args_t...>
            {
                return jst::contrib::make_closure(_cpo{}, (args_t &&)args...);
            }
        };

        template <sequence_label_kind label_kind>
        inline constexpr _cpo<label_kind> labelled{};
    } // namespace _tree_adaptor

    using _tree_adaptor::labelled;
}  // namespace libjst
