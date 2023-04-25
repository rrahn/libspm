// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides base implementation of a referentially compressed multisequence (rcms).
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cmath>
#include <ranges>
#include <unordered_map>

#include <seqan3/utility/detail/multi_invocable.hpp>
#include <seqan3/core/concept/cereal.hpp>

#include <libcontrib/type_traits.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/coverage/concept.hpp>
#include <libjst/rcms/contiguous_multimap.hpp>
#include <libjst/rcms/delta_sequence_variant.hpp>
#include <libjst/rcms/generic_delta.hpp>
#include <libjst/rcms/indel_variant.hpp>
#include <libjst/rcms/packed_breakend_key.hpp>

namespace libjst
{
    namespace detail {
        enum struct delta_kind {
            snv = 0b00,
            insertion = 0b01,
            deletion = 0b10,
            indel = insertion | deletion
        };
    } // namespace detail

    // requires is_object_v<source_t>
    // requires breakpoint_coverage<coverage_t>
    template <typename source_t, typename coverage_t>
    class compressed_multisequence { // TODO: breakend multimap

        using value_t = std::ranges::range_value_t<source_t>;
        static constexpr std::array<value_t, 4> _snv_table{{value_t{0}, value_t{1}, value_t{2}, value_t{3}}};

        // so what is the value type?
        // how can we change the different implementations for static and dynamic?
        using position_type = uint32_t;
        using breakend_key_type = packed_breakend_key<position_type>;
        using breakend_map_type = contiguous_multimap<breakend_key_type, coverage_t>;

        using coverage_value_type = std::ranges::range_value_t<coverage_t>;
        using indel_key_type = std::pair<breakend_key_type, coverage_value_type>;

        struct indel_key_hash {
            constexpr std::size_t operator()(indel_key_type const & key) const noexcept {
                return std::hash<breakend_key_type>{}(key.first) ^ std::hash<coverage_value_type>{}(key.second);
            }
        };

        using deletion_type = deletion_element<std::ranges::iterator_t<breakend_map_type>>;
        using insertion_type = insertion_element<source_t>;
        using indel_type = indel_variant<deletion_type, insertion_type>;

        using indel_map_type = std::unordered_map<indel_key_type, indel_type, indel_key_hash>;
        using coverage_domain_type = libjst::coverage_domain_t<coverage_t>;

        template <bool>
        class iterator_impl;

        template <typename>
        class delta_proxy;

        source_t _source{}; // only links to the source.
        breakend_map_type _breakend_map{};
        indel_map_type _indel_map{};
        coverage_domain_type _coverage_domain{};

    public:

        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;
        using value_type = std::iter_value_t<iterator>;
        using source_type = std::span<std::ranges::range_value_t<source_t> const>;
        using size_type = typename breakend_map_type::size_type;

        compressed_multisequence() = default;
        explicit compressed_multisequence(source_t source, coverage_domain_type coverage_domain) :
            _source{std::move(source)},
            _coverage_domain{std::move(coverage_domain)}
        {
            using position_t = libjst::breakend_t<value_type>;
            position_t source_size = std::ranges::size(_source);
            coverage_t def_coverage{std::views::iota(_coverage_domain.min(), _coverage_domain.max()), _coverage_domain};
            _breakend_map.emplace(breakend_key_type{indel_breakend_kind::nil, 0}, def_coverage);
            _breakend_map.emplace(breakend_key_type{indel_breakend_kind::nil, source_size}, def_coverage);
        }

        iterator insert(value_type value) { // low_breakend, alt_sequence, coverage
            if (libjst::get_domain(libjst::coverage(value)) != _coverage_domain)
                throw std::domain_error{"Trying to insert an element from a different coverage domain!"};

            switch (select_delta_kind(value)) {
                case detail::delta_kind::snv: return insert_snv_impl(std::move(value));
                case detail::delta_kind::insertion: return insert_insertion_impl(std::move(value));
                case detail::delta_kind::deletion: return insert_deletion_impl(std::move(value));
                case detail::delta_kind::indel: {
                    insert_insertion_impl(value);
                    return insert_deletion_impl(std::move(value));
                }
                default: throw std::runtime_error{"Unkown delta kind."};
            }
        }

        // Check if given value is conflicting with some of the variants!
        constexpr bool has_conflicts(value_type const & value) const  noexcept {
            // find equal range in position
            auto candidates = std::ranges::equal_range(std::ranges::next(_breakend_map.begin()),
                                                       std::ranges::prev(_breakend_map.end()),
                                                       libjst::low_breakend(value),
                                                       std::ranges::less{},
                                                       [] (auto && breakend) {
                                                          return breakend.first.position();
                                                       });
            for (auto && breakend : candidates) {
                if (!libjst::coverage_intersection(libjst::coverage(value), breakend.second).empty())
                    return true;
            }
            return false;
        }

        // void erase(const_iterator first = end(), const_iterator last = end()) {

        // }

        // iterator find(key_type const key) {

        // }

        // constexpr compressed_multisequence extract(coverage_t) const noexcept {
        //     return leaves this without the coverages here and returns a new compressed mul
        // }

        // constexpr source_t const & merge(compressed_multisequence()) const noexcept {
        //      collision: like same coverage element in both multisequences.
        //     return _source;
        // }

        constexpr source_type source() const noexcept {
            return _source;
        }

        constexpr size_t size() const noexcept {
            return std::ranges::size(_breakend_map);
        }

        constexpr void reserve(size_type const new_capacity) {
            return _breakend_map.reserve(new_capacity);
        }

        constexpr coverage_domain_type const & coverage_domain() const noexcept {
            return _coverage_domain;
        }

        constexpr iterator begin() noexcept {
            return get_iterator(_breakend_map.begin());
        }

        constexpr const_iterator begin() const noexcept {
            return get_iterator(_breakend_map.begin());
        }

        constexpr iterator end() noexcept {
            return get_iterator(_breakend_map.end());
        }

        constexpr const_iterator end() const noexcept {
            return get_iterator(_breakend_map.end());
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_source, _breakend_map, /*_indel_map,*/ _coverage_domain);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_source, _breakend_map, /*_indel_map,*/ _coverage_domain);
        }

    private:

        constexpr detail::delta_kind select_delta_kind(value_type const & value) {
            detail::delta_kind kind{detail::delta_kind::snv};
            if (libjst::breakpoint_span(value) != 1 || std::ranges::size(libjst::alt_sequence(value)) != 1) {
                using underlying_kind_t = std::underlying_type_t<detail::delta_kind>;
                if (libjst::breakpoint_span(value) == 0) {
                    kind = detail::delta_kind::insertion;
                }
                if (std::ranges::empty(libjst::alt_sequence(value))) {
                    kind = static_cast<detail::delta_kind>(
                                static_cast<underlying_kind_t>(kind) |
                                static_cast<underlying_kind_t>(detail::delta_kind::deletion));
                }
            }

            return kind;
        }

        constexpr iterator get_iterator(std::ranges::iterator_t<breakend_map_type> it) noexcept {
            return iterator{std::move(it), std::addressof(_indel_map)};
        }

        constexpr const_iterator get_iterator(std::ranges::iterator_t<breakend_map_type const> it) const noexcept {
            return const_iterator{std::move(it), std::addressof(_indel_map)};
        }

        template <typename code_t, typename fwd_value_t>
        std::ranges::iterator_t<breakend_map_type> insert_breakend(code_t code, fwd_value_t && value) {
            using underlying_type = typename breakend_key_type::underlying_type;

            underlying_type position{libjst::low_breakend((fwd_value_t &&)value)};
            if constexpr (std::same_as<code_t, indel_breakend_kind>) {
                if (code == indel_breakend_kind::deletion_high) {
                    position = libjst::high_breakend((fwd_value_t &&)value);
                }
            }
            return _breakend_map.emplace_hint(std::ranges::prev(_breakend_map.end()),
                                              breakend_key_type{code, position},
                                              libjst::coverage((fwd_value_t &&) value));
        }

        iterator insert_snv_impl(value_type value) {
            auto breakend_it = insert_breakend(libjst::alt_sequence(value)[0], std::move(value));
            return get_iterator(std::move(breakend_it));
        }

        iterator insert_deletion_impl(value_type value) {
            // coverages are linked
            auto low_it = insert_breakend(indel_breakend_kind::deletion_low, value);
            auto high_it = insert_breakend(indel_breakend_kind::deletion_high, std::move(value));

            deletion_type low_deletion{high_it};
            deletion_type high_deletion{low_it};
            _indel_map.emplace(indel_key_type{low_it->first, low_it->second.front()}, low_deletion);
            _indel_map.emplace(indel_key_type{high_it->first, high_it->second.front()}, high_deletion);

            // TODO: fix resolving the updates of the deletion variants.

            return get_iterator(std::move(low_it));
        }

        iterator insert_insertion_impl(value_type value) {
            insertion_type insertion{std::move(libjst::alt_sequence(value))};
            auto breakend_it = insert_breakend(indel_breakend_kind::insertion_low, std::move(value));
            _indel_map.emplace(indel_key_type{breakend_it->first, breakend_it->second.front()}, std::move(insertion));

            return get_iterator(std::move(breakend_it));
        }

    };

    template <typename source_t, typename coverage_t>
    template <bool is_const>
    class compressed_multisequence<source_t, coverage_t>::iterator_impl {

        friend compressed_multisequence;

        template <bool>
        friend class iterator_impl;

        using maybe_const_map_type = std::conditional_t<is_const, breakend_map_type const, breakend_map_type>;
        using breakend_iterator = std::ranges::iterator_t<maybe_const_map_type>;

        breakend_iterator _breakend_it{};
        indel_map_type const * _indel_map{};

        explicit constexpr iterator_impl(breakend_iterator breakend_it, indel_map_type const * indel_map) noexcept :
            _breakend_it{std::move(breakend_it)},
            _indel_map{indel_map}
        {}

    public:

        using value_type = generic_delta<source_t, coverage_t>;
        using reference = delta_proxy<breakend_iterator>;
        using difference_type = std::iter_difference_t<breakend_iterator>;
        using pointer = void;
        using iterator_category = std::random_access_iterator_tag;

        constexpr iterator_impl() = default;
        constexpr iterator_impl(iterator_impl<!is_const> other) noexcept requires is_const :
            _breakend_it{std::move(other._breakend_it)},
            _indel_map{other._indel_map}
        {}

        constexpr reference operator*() const noexcept {
            return reference{*_breakend_it, *_indel_map};
        }

        constexpr reference operator[](difference_type const step) const noexcept {
            return *(this->operator+(step));
        }

        constexpr iterator_impl & operator++() noexcept {
            ++_breakend_it;
            return *this;
        }

        constexpr iterator_impl operator++(int) noexcept {
            iterator_impl tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator_impl & operator+=(difference_type const step) noexcept {
            _breakend_it += step;
            return *this;
        }

        constexpr iterator_impl & operator--() noexcept {
            --_breakend_it;
            return *this;
        }

        constexpr iterator_impl operator--(int) noexcept {
            iterator_impl tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator_impl & operator-=(difference_type const step) noexcept {
            _breakend_it -= step;
            return *this;
        }

    private:

        constexpr friend iterator_impl operator+(iterator_impl lhs, difference_type const step) noexcept {
            return lhs += step;
        }

        constexpr friend iterator_impl operator+(difference_type const step, iterator_impl rhs) noexcept {
            return rhs + step;
        }

        constexpr friend iterator_impl operator-(iterator_impl lhs, difference_type const step) noexcept {
            return lhs -= step;
        }

        constexpr friend difference_type operator-(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._breakend_it - rhs._breakend_it;
        }

        constexpr friend bool operator==(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._breakend_it == rhs._breakend_it;
        }

        constexpr friend std::strong_ordering operator<=>(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._breakend_it <=> rhs._breakend_it;
        }
    };

    template <typename source_t, typename coverage_t>
    template <typename breakend_iterator>
    class compressed_multisequence<source_t, coverage_t>::delta_proxy {

        friend compressed_multisequence;

        using breakend_reference_t = std::iter_reference_t<breakend_iterator>;
        using reference = breakend_reference_t;
        using sequence_reference = source_type;

        breakend_reference_t _breakend_reference;
        indel_map_type const & _indel_map;

        explicit constexpr delta_proxy(breakend_reference_t breakend_reference, indel_map_type const & indel_map) noexcept :
            _breakend_reference{breakend_reference},
            _indel_map{indel_map}
        {}

    public:

        delta_proxy() = delete;

        constexpr operator value_type() const noexcept {
            return value_type{libjst::breakpoint(*this), libjst::alt_sequence(*this), libjst::coverage(*this)};
        }

        constexpr breakend_key_type get_key() const noexcept {
            return _breakend_reference.first;
        }

        constexpr breakpoint_end get_breakpoint_end() const noexcept {
            return _breakend_reference.first.visit(seqan3::detail::multi_invocable{
                [pos = _breakend_reference.first.position()] (indel_breakend_kind code) {
                    switch (code) {
                        case indel_breakend_kind::deletion_low: [[fallthrough]];
                        case indel_breakend_kind::insertion_low: return breakpoint_end::low;
                        case indel_breakend_kind::nil: return (pos == 0) ? breakpoint_end::low : breakpoint_end::high;
                        default: return breakpoint_end::high;
                    }
                },
                [] (...) {
                    return breakpoint_end::low;
                }
            });
        }

        constexpr std::optional<iterator_impl<true>> jump_to_mate() const noexcept {
            using mate_iterator = iterator_impl<true>;
            using optional_mate = std::optional<mate_iterator>;
            return _breakend_reference.first.visit(seqan3::detail::multi_invocable{
                [&] (indel_breakend_kind code) {
                    switch (code) {
                        case indel_breakend_kind::deletion_low: [[fallthrough]];
                        case indel_breakend_kind::deletion_high: {
                            indel_key_type key{_breakend_reference.first, _breakend_reference.second.front()};
                            return _indel_map.find(key)->second.visit(seqan3::detail::multi_invocable{
                                [&] (deletion_type const & deletion) -> optional_mate {
                                    return iterator_impl<true>{deletion.value(), std::addressof(_indel_map)};
                                },
                                [] (insertion_type const &) -> optional_mate { return std::nullopt; }
                            });
                        }
                        default: return optional_mate{std::nullopt};
                    }
                },
                [] (...) {
                    return std::nullopt;
                }
            });
        }

    private:

        constexpr breakend_reference_t get_breakend_mate(indel_key_type key) const noexcept {
            using optional_it_t = std::optional<breakend_iterator>;
            assert(_indel_map.contains(key));
            optional_it_t mate = _indel_map.find(key)->second.visit(seqan3::detail::multi_invocable{
                [] (deletion_type const & deletion) -> optional_it_t { return deletion.value(); },
                [&] (insertion_type const &) -> optional_it_t { return std::nullopt; }
            });
            if (mate.has_value()) {
                return *(mate.value());
            } else {
                return _breakend_reference;
            }
        }

        constexpr breakpoint get_deletion_breakpoint(indel_breakend_kind const deletion_kind) const noexcept {
            assert(deletion_kind == indel_breakend_kind::deletion_low ||
                   deletion_kind == indel_breakend_kind::deletion_high);

            using signed_position_t = std::make_signed_t<typename breakend_key_type::underlying_type>;
            using breakend_t = typename breakpoint::value_type;

            indel_key_type indel_key{_breakend_reference.first, _breakend_reference.second.front()};
            breakend_reference_t breakend_mate = get_breakend_mate(std::move(indel_key));

            breakend_t low_breakend = (deletion_kind == indel_breakend_kind::deletion_low) ?
                                        _breakend_reference.first.position() :
                                        breakend_mate.first.position();


            size_t deletion_size = std::abs(static_cast<signed_position_t>(breakend_mate.first.position()) -
                                            static_cast<signed_position_t>(_breakend_reference.first.position()));
            return breakpoint{low_breakend, deletion_size};
        }

        constexpr breakpoint extract_breakpoint() const noexcept {
            // what can we have:
            using position_t = breakend_key_type::underlying_type;
            breakend_key_type key = _breakend_reference.first;
            position_t const position = key.position();
            return key.visit(seqan3::detail::multi_invocable{
                [&] (indel_breakend_kind indel_kind) {
                    if (indel_kind == indel_breakend_kind::deletion_low ||
                        indel_kind == indel_breakend_kind::deletion_high) {
                        return get_deletion_breakpoint(indel_kind);
                    } else {
                        return breakpoint{position, 0};
                    }
                },
                [&] (auto const &) { return breakpoint{position, 1}; }
            });
        }

        constexpr sequence_reference extract_alt_sequence() const noexcept {
            // what can we have:
            return _breakend_reference.first.visit(seqan3::detail::multi_invocable{
                [&] (indel_breakend_kind breakend_kind) {
                    indel_key_type indel_key{_breakend_reference.first, _breakend_reference.second.front()};
                    if (breakend_kind != indel_breakend_kind::nil) {
                        assert(_indel_map.contains(indel_key));
                        return _indel_map.find(indel_key)->second.visit(seqan3::detail::multi_invocable{
                            [] (insertion_type const & insertion) { return sequence_reference{insertion.value()}; },
                            [] (deletion_type const &) { return sequence_reference{}; }
                        });
                    } else {
                        return sequence_reference{};
                    }
                },
                [&] (...) {
                    return sequence_reference{std::addressof(_snv_table[(int)_breakend_reference.first.snv_value()]), 1};
                }
            });
        }

        template <typename result_t>
        using make_result_t = std::conditional_t<std::is_rvalue_reference_v<result_t>,
                                                 std::remove_reference_t<result_t>,
                                                 result_t>;

        template <typename cpo_t>
            requires std::tag_invocable<cpo_t, breakpoint>
        friend constexpr auto tag_invoke(cpo_t cpo, delta_proxy me)
            noexcept(std::is_nothrow_tag_invocable_v<cpo_t, breakpoint>)
            -> make_result_t<std::tag_invoke_result_t<cpo_t, breakpoint>>
        {
            return std::tag_invoke(cpo, libjst::get_breakpoint(std::move(me)));
        }

        friend constexpr breakpoint tag_invoke(std::tag_t<libjst::get_breakpoint>, delta_proxy me) noexcept
        {
            return me.extract_breakpoint();
        }

        friend constexpr sequence_reference tag_invoke(std::tag_t<libjst::alt_sequence>, delta_proxy me) noexcept
        {
            return me.extract_alt_sequence();
        }

        friend constexpr coverage_t const & tag_invoke(std::tag_t<libjst::coverage>, delta_proxy me) noexcept
        {
            return me._breakend_reference.second;
        }

        friend constexpr position_type tag_invoke(std::tag_t<libjst::position>, delta_proxy me) noexcept
        {
            return me._breakend_reference.first.position();
        }

        friend constexpr alternate_sequence_kind tag_invoke(std::tag_t<libjst::alt_kind>, delta_proxy me) noexcept
        {
            return me._breakend_reference.first.visit(seqan3::detail::multi_invocable{
                [&] (indel_breakend_kind indel_kind) {
                    switch (indel_kind) {
                        case indel_breakend_kind::insertion_low: return alternate_sequence_kind::insertion;
                        case indel_breakend_kind::deletion_low: [[fallthrough]];
                        case indel_breakend_kind::deletion_high: return alternate_sequence_kind::deletion;
                        default: return alternate_sequence_kind::unknown;
                    }
                },
                [] (...) { return alternate_sequence_kind::replacement; }
            });
        }

        friend constexpr std::ptrdiff_t tag_invoke(std::tag_t<libjst::effective_size>, delta_proxy me) noexcept
        {
            return std::ranges::ssize(libjst::alt_sequence(me)) -
                   static_cast<std::ptrdiff_t>(libjst::breakpoint_span(libjst::get_breakpoint(me)));
        }
    };

}  // namespace libjst
