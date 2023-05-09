// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an adapter to make the pigeonhole online pattern matching algorithm work with the JST.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan/index.h>

#include <libcontrib/seqan/container_adapter.hpp>

#include <libjst/matcher/seqan_pattern_base.hpp>

namespace seqan {

    template <typename TRange, typename TFinderSpec, typename TIteratorSpec>
    struct Iterator<Finder<jst::contrib::seqan_container_adapter<TRange>, TFinderSpec> const, TIteratorSpec>
    {
        using Type = typename Iterator<jst::contrib::seqan_container_adapter<TRange> const, TIteratorSpec>::Type;
    };

    using PigeonholeSeedOnlyTag = seqan::Tag<struct SeedOnly>;
    template <>
    struct Pigeonhole<PigeonholeSeedOnlyTag>
    {
        enum { ONE_PER_DIAGONAL = 0 };    // disable for heuristic.
        enum { HAMMING_ONLY = 0 };        // 1..no indels; 0..allow indels
    };

    using PigeonholeSeedOnly = Pigeonhole<PigeonholeSeedOnlyTag>;

    struct PigeonholeSeedOnlyPosition {
        std::ptrdiff_t index{};
        std::ptrdiff_t offset{};
        std::ptrdiff_t count{};

    private:

        constexpr friend bool operator==(PigeonholeSeedOnlyPosition const &, PigeonholeSeedOnlyPosition const &) noexcept = default;

        template <typename stream_t, typename me_t>
            requires std::same_as<std::remove_cvref_t<me_t>, PigeonholeSeedOnlyPosition>
        friend stream_t & operator<<(stream_t & stream, me_t && me) {
            stream << "<" << me.index << ", " << me.offset << ", " << me.count << ">";
            return stream;
        }
    };

    template <typename TIndex>
    struct Position<Pattern<TIndex, PigeonholeSeedOnly>> {
        using Type = PigeonholeSeedOnlyPosition;
    };

    template <typename THaystack, typename TPattern>
    struct FindResult<Finder<THaystack, PigeonholeSeedOnly>, TPattern>
    {
        using Type = SwiftHit_<int64_t>;
    };

    template <typename TFinder, typename TIndex, typename THValue>
    inline bool _pigeonholeProcessQGram(TFinder &finder,
                                        Pattern<TIndex, PigeonholeSeedOnly> &pattern,
                                        THValue hash)
    {
        typedef typename Fibre<TIndex, QGramSA>::Type const TSA;
        typedef typename Iterator<TSA, Standard>::Type      TSAIter;
        typedef typename TFinder::TPigeonholeHit            THit;

        TIndex const &index = host(pattern);

        // all previous matches reported -> search new ones
        clear(finder.hits);

        TSAIter saBegin = begin(indexSA(index), Standard());
        Pair<unsigned> ndlPos;
        THit hit;

        unsigned bktNo = getBucket(index.bucketMap, hash);
        TSAIter occ = saBegin + indexDir(index)[bktNo];
        TSAIter occEnd = saBegin + indexDir(index)[bktNo + 1];

        for(; occ != occEnd; ++occ)
        {
            // Haystack occurrence site
            hit.hstkPos = finder.curPos;
            hit.bucketWidth = length(getFibre(index, QGramShape{}));
            // Needle occurrence site
            posLocalize(ndlPos, *occ, stringSetLimits(index));
            hit.ndlSeqNo = getSeqNo(ndlPos);
            hit.ndlPos = getSeqOffset(ndlPos);
            hit.hitLengthNeedle = hit.bucketWidth;
            appendValue(finder.hits, hit);
        }

        finder.curHit = begin(finder.hits, Standard());
        finder.endHit = end(finder.hits, Standard());

        return !empty(finder.hits);
    }

    template <typename TFinder, typename TIndex>
    inline void _copyPigeonholeHit(TFinder &finder, Pattern<TIndex, PigeonholeSeedOnly> &pattern)
    {
        auto const & hit = *finder.curHit;
        pattern.curSeqNo = hit.ndlSeqNo;
        pattern.curBeginPos = hit.ndlPos;
        pattern.curEndPos = hit.ndlPos + hit.hitLengthNeedle;
    }

    template <typename TIndex>
    inline PigeonholeSeedOnlyPosition position(Pattern<TIndex, PigeonholeSeedOnly> const & pattern) noexcept
    {
        return {.index = static_cast<std::ptrdiff_t>(pattern.curSeqNo),
                .offset = static_cast<std::ptrdiff_t>(pattern.curBeginPos),
                .count = static_cast<std::ptrdiff_t>(pattern.curEndPos - pattern.curBeginPos)};
    }
} // namespace seqan

namespace libjst
{

    template <std::ranges::random_access_range needle_t>
    class pigeonhole_matcher : public seqan_pattern_base<pigeonhole_matcher<needle_t>>
    {
    private:

        using base_t = seqan_pattern_base<pigeonhole_matcher<needle_t>>;

        friend base_t;

        using compatible_needle_type = jst::contrib::seqan_container_t<needle_t>;
        using multi_needle_type = seqan::StringSet<compatible_needle_type>;
        using qgram_shape_type = seqan::Shape<std::ranges::range_value_t<compatible_needle_type>, seqan::SimpleShape>;
        using finder_spec_type = seqan::PigeonholeSeedOnly;
        using index_type = seqan::Index<multi_needle_type, seqan::IndexQGram<qgram_shape_type, seqan::OpenAddressing>>;
        using pattern_type = seqan::Pattern<index_type, finder_spec_type>;

        multi_needle_type _multi_needle{};
        index_type _needle_index{_multi_needle};
        pattern_type _pattern{_needle_index};
        double _error_rate{};

    public:

        pigeonhole_matcher() = delete;
        template <std::ranges::viewable_range _needle_t>
            requires (!std::same_as<_needle_t, pigeonhole_matcher> &&
                       std::constructible_from<compatible_needle_type, _needle_t>)
        explicit pigeonhole_matcher(_needle_t && needle, double error_rate = 0.0) :
            _error_rate{error_rate}
        {
            appendValue(getFibre(_needle_index, seqan::QGramText{}),
                        jst::contrib::make_seqan_container(std::views::all((_needle_t &&) needle)));
            _patternInit(_pattern, _error_rate);
        }

        template <std::ranges::viewable_range _multi_needle_t>
            requires (!std::same_as<_multi_needle_t, pigeonhole_matcher> &&
                       std::constructible_from<compatible_needle_type,
                                               std::ranges::range_reference_t<_multi_needle_t>>)
        explicit pigeonhole_matcher(_multi_needle_t && multi_needle, double error_rate = 0.0) :
            _error_rate{error_rate}
        {
            for (auto && needle : multi_needle)
                appendValue(getFibre(_needle_index, seqan::QGramText{}),
                            jst::contrib::make_seqan_container(std::views::all((decltype(needle) &&) needle)));

            _patternInit(_pattern, _error_rate);
        }

        constexpr auto position() const noexcept {
            return seqan::position(_pattern);
        }
    private:

        template <typename haystack_t>
        constexpr auto make_finder(haystack_t & haystack) const noexcept
        {
            // TODO: configure repeat length and finder.
            if (std::ranges::size(haystack) > 1000)
                return seqan::Finder<haystack_t, finder_spec_type>{haystack, 1000, 1};
            else
                return seqan::Finder<haystack_t, finder_spec_type>{haystack};
        }

        constexpr pigeonhole_matcher & get_pattern() noexcept {
            return *this;
        }

        constexpr auto custom_find_arguments() const noexcept {
            return std::tuple{_error_rate};
        }

        constexpr friend std::size_t tag_invoke(std::tag_t<window_size>, pigeonhole_matcher const & me) noexcept {
            return length(getFibre(needle(me._pattern), seqan::QGramShape{}));
        }

        template <typename haystack_t>
        constexpr bool initialise(seqan::Finder<haystack_t, finder_spec_type> & finder,
                                  pattern_type & pattern)
        {
            pattern.finderLength = std::ranges::size(haystack(finder));
            pattern.finderPosOffset = 0;
            pattern.finderPosNextOffset = pattern.maxSeqLen + pattern.finderLength;

            _finderSetNonEmpty(finder);
            finder.dotPos = 100000;
            finder.dotPos2 = 10 * finder.dotPos;

            if (!_firstNonRepeatRange(finder, pattern))
                return false;

            clear(finder.hits);
            if (_pigeonholeProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
                _copyPigeonholeHit(finder, pattern);

            return true;
        }

        template <typename haystack_t, typename ...args_t>
        friend bool find(seqan::Finder<haystack_t, finder_spec_type> & finder,
                         pigeonhole_matcher & matcher,
                         args_t && ...args)
        {
            pattern_type & pattern = matcher._pattern;
            if (empty(finder)) {
                if (!matcher.initialise(finder, pattern)) {
                    return false;
                } else if (finder.curHit != finder.endHit) {
                    return true;
                }
            }

            return find(finder, pattern, (args_t &&)args...);
        }
    };

    template <std::ranges::viewable_range needle_t>
    pigeonhole_matcher(needle_t &&) -> pigeonhole_matcher<std::views::all_t<needle_t>>;

    template <std::ranges::viewable_range needle_t>
    pigeonhole_matcher(needle_t &&, double) -> pigeonhole_matcher<std::views::all_t<needle_t>>;

    template <std::ranges::viewable_range multi_needle_t>
        requires std::ranges::random_access_range<std::ranges::range_reference_t<multi_needle_t>>
    pigeonhole_matcher(multi_needle_t &&) -> pigeonhole_matcher<std::views::all_t<std::ranges::range_reference_t<multi_needle_t>>>;

    template <std::ranges::viewable_range multi_needle_t>
        requires std::ranges::random_access_range<std::ranges::range_reference_t<multi_needle_t>>
    pigeonhole_matcher(multi_needle_t &&, double) -> pigeonhole_matcher<std::views::all_t<std::ranges::range_reference_t<multi_needle_t>>>;

}  // namespace libjst
