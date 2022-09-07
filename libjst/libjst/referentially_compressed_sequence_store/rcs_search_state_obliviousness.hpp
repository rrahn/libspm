// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides standard jst search.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
    template <typename rcs_store_t>
    class rcs_search_state_oblivious // realises burrowed_range
    {
    private:
        using jst_traversal_t = jst_traversal_t<rcs_store_t>;

        std::shared_ptr<jst_traversal_t> _jst_traversal{};

        class iterator;

        explicit constexpr rcs_search_state_oblivious(jst_traversal_t * other_traversal) noexcept :
            _jst_traversal{other_traversal}
        {}

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr rcs_search_state_oblivious() = default; //!< Default.

        template <typename pattern_t>
        explicit constexpr rcs_search_state_oblivious(rcs_store_t const & rcs_store,
                                                      pattern_t const & pattern) noexcept :
            _jst_traversal{std::make_shared<jst_traversal_t>(rcs_store, libjst::window_size(pattern))}
        {}
        //!\}

        constexpr bool empty() noexcept {
            return iterator{_jst_traversal.release()};
        }

        constexpr iterator begin() noexcept {
            return iterator{std::move(_jst_traversal)};
        }

        constexpr std::default_sentinel_t end() noexcept {
            return std::default_sentinel;
        }
    private:

        template <typename pattern_t>
            requires std::predicate<pattern_t, typename context_traits<traversal_context_t>::sequence_type>
        constexpr friend rcs_search_state_oblivious tag_invoke(std::tag_t<search>,
                                                               rcs_search_state_oblivious && me,
                                                               pattern_t && pattern) {
            for (auto && ctxt : me) {
                if (std::invoke(pattern, ctxt.sequence()))
                    break;
            }
            return rcs_search_state_oblivious{me._jst_traversal};
        }
    };

    template <typename rcs_store_t>
    class rcs_search_state_oblivious<rcs_store_t>::iterator {
    private:

        std::shared_ptr<jst_traversal_t> _jst_traversal{};

        friend rcs_store_t;

        explicit iterator(jst_traversal_t && jst_traversal) noexcept :
            _jst_traversal{std::make_shared<jst_traversal_t>(std::move(jst_traversal))}
        {}

    public:
        using value_type = typename jst_traversal_t::context_type;
        using reference = value_type const &;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type const *;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;

        constexpr reference operator*() const noexcept {
            return _jst_traversal->context();
        }

        constexpr pointer operator->() const noexcept {
            return std::addressof(_jst_traversal->context());
        }

        constexpr iterator & operator++() noexcept {

            return *this;
        }

        constexpr void operator++() const noexcept {

            return *this;
        }

        constexpr friend bool operator==(iterator const & lhs, std::default_sentinel_t const &) const noexcept {
            return !*lhs._jst_traversal; // reached sink
        }
    };
}  // namespace libjst
