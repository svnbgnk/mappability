#include <seqan/index.h>

using namespace std;
using namespace seqan;

namespace seqan {

template <typename TChar, typename TOwner>
struct SAValue<StringSet<String<TChar>, TOwner > >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TOwner>
struct SAValue<String<TChar, TOwner > >
{
    typedef uint32_t Type;
};

};

typedef StringSet<DnaString, Owner<ConcatDirect<> > > TText;
typedef Concatenator<TText>::Type TConcat;
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;
typedef Index<StringSet<String<Dna, MMap<> >, Owner<ConcatDirect<> > >, TIndexConfig> TIndex;
typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

auto mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

int main() {
    TIndex index;
    open(index, toCString("/srv/public/cpockrandt/schemes-index-hg38-without-N/index"), OPEN_RDONLY);
    TIter iter(index);
    cout << mytime() << "Index loaded\n";

    constexpr unsigned errors = 2;
    unsigned m = 36;

    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, m);

    auto & genome = indexText(index);
    TConcat concatGenome = concat(genome);
    uint64_t genome_length = length(concatGenome);
    cout << mytime() << "Genome length: " << genome_length << "\n";
    cout << "Solving: (" << m << ", " << errors << ")\n"; 
    vector<uint16_t> count;
    count.resize(genome_length - m + 1);
    cout << mytime() << "Vector built.\n";
    #pragma omp parallel for schedule(dynamic, 1000000)
    for (uint64_t i = 0; i < genome_length - m + 1; ++i)
    {
        unsigned hits = 0;
        auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
            hits += countOccurrences(it);
        };
        
        auto const & needle = infix(concatGenome, i, i + m);
        goRoot(iter);
        _optimalSearchScheme(delegate, iter, needle, scheme, HammingDistance());
        count[i] = hits;
    }
    cout << mytime() << "Done.\n";

    std::ofstream outfile("/group/ag_abi/cpockrandt/mappability_results_2_36"/* + errors*/, std::ios::out | std::ofstream::binary);
    std::copy(count.begin(), count.end(), std::ostream_iterator<int>(outfile));
    outfile.close();
    cout << mytime() << "Saved to disk.\n";
}

