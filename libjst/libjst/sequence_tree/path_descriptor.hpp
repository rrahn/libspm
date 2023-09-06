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

#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <iosfwd>
#include <iterator>

#include <cereal/types/array.hpp>

#include <seqan3/core/concept/cereal.hpp>

namespace libjst
{

    template <typename word_t = uint_fast64_t>
    class extended_word {
    private:

        static constexpr size_t word_count = 4;
        static constexpr size_t bits_per_word = sizeof(word_t) << 3;
        static constexpr word_t carry_msb_mask = static_cast<word_t>(1) << (bits_per_word - 1);

        template <size_t idx, typename value_t>
        struct indexed_word {
            value_t value;

            static constexpr size_t offset = bits_per_word * idx;

            constexpr size_t bit_width() const noexcept {
                return std::bit_width(value) + (value != 0) * offset;
            }
        };

        using data_type = std::array<word_t, word_count>;
        data_type _data{};

    public:

        using value_type = word_t;

        constexpr extended_word() = default;
        constexpr explicit extended_word(word_t init_value) : _data{init_value, 0, 0, 0}
        {}

        constexpr bool operator[](std::ptrdiff_t const index) const noexcept {
            assert(index < static_cast<std::ptrdiff_t>(size()));
            auto [word_idx, word_offset] = to_local_index(index);
            return std::bitset<bits_per_word>{_data[word_idx]}[word_offset];
        }

        constexpr extended_word & operator<<=(size_t const shift) noexcept {
            assert(shift < bits_per_word);
            word_t const shift_mask = (1 << shift) - 1;
            word_t carry_bits{};
            for (word_t & word : _data) {
                word_t next_carry = std::rotl(word, shift) & shift_mask;
                word = (word << shift) | carry_bits;
                carry_bits = next_carry;
            }
            return *this;
        }

        constexpr extended_word & operator|=(word_t const & rhs) noexcept {
            _data[0] |= rhs;
            return *this;
        }

        constexpr size_t size() const noexcept {
            return unpack_apply(_data, [] (auto ...words) {
                return std::max({words.bit_width()...});
            });
        }

        static constexpr size_t max_size() noexcept {
            return bits_per_word * word_count;
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t const & oarchive) const
        {
            oarchive(_data);
        }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t &iarchive) noexcept
        {
            iarchive(_data);
        }

        word_t const* data() const noexcept {
            return _data.data();
        }
    private:

        static constexpr std::pair<size_t, size_t> to_local_index(size_t const index) noexcept {
            return {index >> 6, index & (bits_per_word - 1)};
        }

        template <typename data_t, typename fn_t>
        static constexpr auto unpack_apply(data_t && data, fn_t && apply) noexcept {
            return []
            <typename _data_t, typename _fn_t, size_t ...idx>
            (_data_t && data, _fn_t && apply, std::index_sequence<idx...> const &) constexpr {
                return apply(indexed_word<idx, std::ranges::range_reference_t<_data_t>>{data[idx]}...);
            } ((data_t &&)data, (fn_t &&)apply, std::make_index_sequence<word_count>());
        }

        friend constexpr bool operator==(extended_word const &, extended_word const &) noexcept = default;
        friend constexpr std::strong_ordering operator<=>(extended_word const &, extended_word const &) noexcept = default;
    };
    class alternate_path_descriptor {
    protected:

        using data_type = extended_word<uint_fast64_t>;

        static constexpr bool ref_mask{0};
        static constexpr bool alt_mask{1};

        class iterator {
        private:

            friend alternate_path_descriptor;

            data_type const * _path{};
            size_t _active_bit{};

            explicit constexpr iterator(data_type const * word, size_t const active_bit) noexcept :
                _path{word},
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
                assert(_path != nullptr);
                return (*_path)[_active_bit - 1];
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

        data_type _word{1};

    public:

        using underlying_type = data_type;

        constexpr alternate_path_descriptor() = default;

        constexpr void next() noexcept {
            assert(get_word().size() < get_word().max_size());
            get_word() <<= 1;
        }

        constexpr size_t size() const noexcept {
            return get_word().size();
        }

        static constexpr size_t max_size() noexcept {
            return data_type::max_size();
        }

        constexpr void set_alt() noexcept {
            get_word() |= alt_mask;
        }

        constexpr void set_ref() noexcept {
            get_word() |= ref_mask;
        }

        constexpr iterator begin() const noexcept {
            return iterator{std::addressof(get_word()), size()};
        }

        constexpr iterator end() const noexcept {
            return iterator{std::addressof(get_word()), 0};
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

        typename data_type::value_type const* data() const noexcept {
            return get_word().data();
        }

    private:

        constexpr data_type & get_word() noexcept {
            return _word;
        }

        constexpr data_type const & get_word() const noexcept {
            return _word;
        }

        friend constexpr bool operator==(alternate_path_descriptor const &,
                                         alternate_path_descriptor const &) noexcept = default;

        friend constexpr std::strong_ordering operator<=>(alternate_path_descriptor const &,
                                                          alternate_path_descriptor const &) noexcept = default;
    };

    template <typename char_t, typename char_traits_t, typename alt_path_descriptor_t>
        requires std::same_as<std::remove_cvref_t<alt_path_descriptor_t>, alternate_path_descriptor>
    inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                                  alt_path_descriptor_t && descriptor)
    {
        stream << "[";
        if (descriptor.size() != 0) {
            auto it = descriptor.begin();
            for (auto it = descriptor.begin(); it != std::ranges::prev(descriptor.end(), 1, descriptor.begin()); ++it) {
                stream << *it << ", ";
            }
            stream << *it;
        }
        stream << "]";
        return stream;
    }

}  // namespace libjst
