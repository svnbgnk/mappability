#ifndef COMMON_H_
#define COMMON_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;

// reduce space consumption.
// requires genome not to have more than ~ 4 gigabases
// multi-sequence fasta file must contain less than ~64k sequences of at most ~ 4 gigabases in total
namespace seqan {

template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};
template <typename TSpec = void, typename TLengthSum = size_t, unsigned LEVELS = 2, unsigned WORDS_PER_BLOCK = 1>
struct FastFMIndexConfigS1
{
     typedef TLengthSum                                                                         LengthSum;
     typedef Levels<void, LevelsPrefixRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >   Bwt;
     typedef Levels<void, LevelsRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >         Sentinels;
     static const unsigned SAMPLING = 1;
};

};

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}



using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;
/*
using TMyFastConfig = seqan::FastFMIndexConfigS1<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;*/

/*
typedef FastFMIndexConfigS1<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;*/


template <typename TText>
using TIndex = seqan::Index<TText, TIndexConfig>;

// template <typename TString, typename TPos, typename TLength>
// auto extractNeedle(TString const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     return infix(text, pos, pos + len);
// }
//
// template <typename TString, typename TOwner, typename TPos, typename TLength>
// auto extractNeedle(seqan::StringSet<TString, TOwner> const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     auto const & limits = stringSetLimits(text);
//     Pair<unsigned, unsigned> res;
//     posLocalize(res, pos, limits);
//     if (res.i2 + len > length(text[res.i1]))
//     {
//         TString r(suffix(text[res.i1], res.i2));
//         r += prefix(text[res.i1 + 2], length(text[res.i1]) - res.i2);
//         return r;
//     }
//     else
//     {
//         return infix(text[res.i1], res.i2, res.i2 + len);
//     }
// }

// template <typename T>
// struct is_int_vector
// {
//     static constexpr bool value = false;
// };
//
// template <uint8_t width_t>
// struct is_int_vector<sdsl::int_vector<width_t> >
// {
//     static constexpr bool value = true;
// };
//
// template <typename T>
// constexpr bool is_int_vector<T>::value;
//
// template <uint8_t width_t>
// constexpr bool is_int_vector<sdsl::int_vector<width_t> >::value;
#endif

