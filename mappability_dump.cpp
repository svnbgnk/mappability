#include <vector>

#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;

template <typename value_t>
void dump(CharString const & inputPath, CharString const & outputPath, uint64_t const runBound)
{
    vector<value_t> v;
    ifstream file(toCString(inputPath), std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        //std::cout << "File size: " << fileSize << '\n';
        //std::cout << "Sizeof: " << sizeof(value_t) << '\n';
        v.resize(fileSize / sizeof(value_t));
        file.seekg(0, std::ios_base::beg);
        file.read((char*) & v[0], fileSize);
        file.close();

        cout << "Load successful\n";

        //cout << "v.size(): " << v.size() << '\n';

        // TODO: progress bar
        ofstream outfile(toCString(outputPath), std::ios::out | std::ofstream::binary);
        copy(v.begin(), v.end(), (std::ostream_iterator<value_t>(outfile), std::ostream_iterator<int>(outfile, " ")));

        if (runBound > 0)
        {
            uint64_t i = 0;
            cout << "Number of runs:\n";
            for (uint64_t j = 0; j < v.size(); ++j)
            {
                while (v[i + j] == 1 && v[i + j] < v.size())
                    ++i;

                if (runBound > 100)
                    cout <<"Pos: " << j << ", run's length: " << i << '\n';

                j += i;
                i = 0;
            }
        }

        cout << "Done.\n";

        return;
    }

    cout << "Something went wrong ...\n";
}

int main(int argc, char *argv[])
{
    // Argument Parser
    // TODO: allow more than 4 gigabases
    ArgumentParser parser("Mappability Dumper");
    addDescription(parser, "Transforms the output of the mappability program into human readable format. Please be aware that the file size will increase significantly and most likely multiply.");

    addOption(parser, ArgParseOption("I", "input", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "input", "gmapp8 gmapp16");
	setRequired(parser, "input");

    addOption(parser, ArgParseOption("O", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("R", "runs", "Output number of runs with mappability value of at least R", ArgParseArgument::INTEGER, "INT"));


    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    uint64_t runBound;
    CharString inputPath, outputPath;
    getOptionValue(inputPath, parser, "input");
    getOptionValue(outputPath, parser, "output");

    if (isSet(parser, "runs"))
        getOptionValue(runBound, parser, "runs");
    else
        runBound = 0;

    string _inputPath = toCString(inputPath);

    if (_inputPath.substr(_inputPath.find_last_of(".") + 1) == "gmapp8")
        dump<uint8_t>(inputPath, outputPath, runBound);
    else
        dump<uint16_t>(inputPath, outputPath, runBound);

    return 0;
}
