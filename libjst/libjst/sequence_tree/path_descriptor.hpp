// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides path descriptor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <bit>
#include <bitset>
#include <cstdint>
#include <iterator>

#include <seqan3/core/concept/cereal.hpp>

namespace libjst
{

    class alternate_path_descriptor {
    protected:

        using word_type = uint_fast64_t;
        using path_type = std::bitset<sizeof(word_type)>;

        static constexpr bool ref_mask{0};
        static constexpr bool alt_mask{1};

        class iterator {
        private:

            friend alternate_path_descriptor;

            path_type _path{};
            size_t _active_bit{};

            explicit constexpr iterator(word_type word, size_t const active_bit) noexcept :
                _path{std::move(word)},
                _active_bit{active_bit}
            {}

        public:
            using value_type = bool;
            using reference = value_type;
            using pointer = void;
            using difference_type = std::ptrdiff_t;
            using iterator_category = std::random_access_iterator_tag;

            constexpr iterator() = default;

            constexpr reference operator*() const noexcept {
                return _path[_active_bit - 1];
            }

            constexpr reference operator[](difference_type const offset) const noexcept {
                return *(*this + offset);
            }

            constexpr iterator& operator++() noexcept {
                --_active_bit;
                return *this;
            }

            constexpr iterator operator++(int) noexcept {
                iterator tmp{*this};
                --_active_bit;
                return tmp;
            }

            constexpr iterator& operator+=(difference_type const offset) noexcept {
                _active_bit -= offset;
                return *this;
            }

            constexpr iterator& operator--() noexcept {
                ++_active_bit;
                return *this;
            }

            constexpr iterator operator--(int) noexcept {
                iterator tmp{*this};
                ++_active_bit;
                return tmp;
            }

            constexpr iterator& operator-=(difference_type const offset) noexcept {
                _active_bit += offset;
                return *this;
            }
        private:

            friend constexpr iterator operator+(iterator lhs, difference_type const offset) noexcept {
                return lhs -= offset;
            }

            friend constexpr iterator operator+(difference_type const offset, iterator rhs) noexcept {
                return rhs -= offset;
            }

            friend constexpr difference_type operator-(iterator const & lhs, iterator const & rhs) noexcept {
                return rhs._active_bit - lhs._active_bit;
            }

            friend constexpr bool operator==(iterator const & lhs, iterator const & rhs) noexcept {
                return lhs._active_bit == rhs._active_bit;
            }

            friend constexpr std::strong_ordering operator<=>(iterator const & lhs, iterator const & rhs) noexcept {
                return lhs._active_bit <=> rhs._active_bit;
            }
        };

        word_type _word{1};

    public:

        constexpr alternate_path_descriptor() = default;

        constexpr void next() noexcept {
            assert((get_word() & 1ull << (max_size() - 1)) == 0);
            get_word() <<= 1;
        }

        constexpr size_t size() const noexcept {
            return std::bit_width(get_word());
        }

        constexpr size_t max_size() const noexcept {
            return sizeof(word_type) << 3;
        }

        constexpr void set_alt() noexcept {
            get_word() |= alt_mask;
        }

        constexpr void set_ref() noexcept {
            get_word() |= ref_mask;
        }

        constexpr iterator begin() const noexcept {
            return iterator{get_word(), size()};
        }

        constexpr iterator end() const noexcept {
            return iterator{get_word(), 0};
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t const & oarchive) const
        {
            oarchive(get_word());
        }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t &iarchive) noexcept
        {
            iarchive(get_word());
        }

    private:

        constexpr word_type & get_word() noexcept {
            return _word;
        }

        constexpr word_type const & get_word() const noexcept {
            return _word;
        }
    };

}  // namespace libjst
