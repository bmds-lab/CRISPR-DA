#ifndef FMT_HEADER_ONLY
#   define FMT_HEADER_ONLY
#endif
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <boost/iostreams/device/mapped_file.hpp>

#include "fmt/format.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::filesystem::path;
using std::filesystem::exists;
using std::filesystem::remove;
using std::filesystem::file_size;
using namespace boost::iostreams;

// Char to binary encoding
const uint8_t  nucleotideIndexOffset = 65;
const vector<uint8_t> nucleotideIndex{ 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
// Binary to char encoding
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

string getLineEnding(path file, uint64_t seqLength)
{
    string lineEnding;
    char currentChar;
    ifstream inFile;
    inFile.open(file);
    inFile.ignore(seqLength);
    for (;;)
    {
        inFile.read(&currentChar, 1);
        if (currentChar == '\r' || currentChar == '\n')
        {
            lineEnding.push_back(currentChar);
        }
        else
        {
            break;
        }  
    }
    inFile.close();
    return lineEnding;
}

uint64_t encode2bit(const char* ptr, uint64_t seqLength)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLength; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[*ptr - nucleotideIndexOffset]) << (j * 2);
        ptr++;
    }
    return signature;
}

string decode2bit(uint64_t signature, uint64_t seqLength)
{
    string sequence = string(seqLength, ' ');
    for (uint64_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

int main(int argc, char** argv)
{
    // Check number of args
    if (argc != 5)
    {
        cerr << fmt::format("Usage: {} <offtarget-sites> <slice-config> <sequence-length> <output-file>\n", argv[0]) << endl;
        exit(1);
    }


    // Check seq length 
    uint64_t seqLength = atoi(argv[3]);
    if (seqLength > 21)
    {
        cerr << "Sequence length is greater than 21, which is the maximum supported currently\n" << endl;
        exit(1);
    }


    // Check offtarget sites exists
    path otFile(argv[1]);
    if (!exists(otFile))
    {
        cerr << fmt::format("Could not find the specified offtarget sites file: {}", otFile.string()) << endl;
        exit(1);
    }


    // Chek offtarget sites file size
    uintmax_t otFileSize = file_size(otFile);
    uint64_t seqLineLength = seqLength + getLineEnding(otFile, seqLength).size();
    if (otFileSize % seqLineLength != 0)
    {
        cerr << fmt::format("fileSize: {}\n", otFileSize);
        cerr << fmt::format("Error: offtargetSites.txt file does is not a multiple of the expected line length ({})\n", seqLineLength);
        cerr << "The sequence length may be incorrect; alternatively, the line endings\n";
        cerr << fmt::format("may be something other than {}, or there may be junk at the end of the file.", getLineEnding(otFile, seqLength)) << endl;
        exit(1);
    }

    // Slice config file exists
    path scFile(argv[2]);
    if (!exists(scFile))
    {
        cerr << fmt::format("Could not find the specified slice config file: {}", scFile.string()) << endl;
        exit(1);
    }


    // Check slice config file size
    uintmax_t scFileSize = file_size(scFile);
    uint64_t scLineLength = seqLength + getLineEnding(scFile, seqLength).size();
    if (scFileSize % scLineLength != 0) {
        cerr << fmt::format("fileSize: {}\n", scFileSize);
        cerr << fmt::format("Error: sliceconfig.txt file does is not a multiple of the expected line length ({})\n", scLineLength);
        cerr << "The sequence length may be incorrect; alternatively, the line endings\n";
        cerr << fmt::format("may be something other than {}, or there may be junk at the end of the file.", getLineEnding(scFile, seqLength)) << endl;
        exit(1);
    }

    // Read in and genereate slice masks
    ifstream scInFile;
    scInFile.open(argv[2], std::ios::in | std::ios::binary);
    vector<vector<uint64_t>> sliceMasks;
    vector<uint64_t> sliceMasksBinary;
    for (string line; std::getline(scInFile, line);)
    {
        vector<uint64_t> mask;
        uint64_t maskBinary = 0ULL;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (line[j] == '1')
            {
                maskBinary |= 1ULL << j;
                mask.push_back(j);
            }   
        }
        sliceMasks.push_back(mask);
        sliceMasksBinary.push_back(maskBinary);
    }
    scInFile.close();
    size_t sliceCount = sliceMasks.size();


    // Begin counting off targets
    uint64_t seqCount = otFileSize / seqLineLength;
    cout << fmt::format("Number of sequences: {}", seqCount) << endl;

    uint64_t globalCount = 0;
    uint64_t offtargetsCount = 0;

    boost::iostreams::mapped_file_source entireDataSet;
    entireDataSet.open(argv[1]);

    ofstream tempSeqSignatures;
    ofstream tempSeqSignaturesOccurrences;
    tempSeqSignatures.open(fmt::format("{}.tmp.1", argv[4]), std::ios::out | std::ios::binary);
    tempSeqSignaturesOccurrences.open(fmt::format("{}.tmp.2", argv[4]), std::ios::out | std::ios::binary);

    uint64_t progressCount = 0;
    uint64_t offtargetId = 0;
    cout << "Counting occurrences..." << endl;
    while (progressCount < seqCount) {
        const char* ptr =  entireDataSet.data() + (progressCount * seqLineLength);
        uint64_t signature = encode2bit(ptr, seqLength);
        // check how many times the off-target appears
        // (assumed the list is sorted)
        uint32_t occurrences = 1;
        while (memcmp(ptr, ptr + (seqLineLength * occurrences), seqLength) == 0) {
            occurrences++;
            if ((seqCount - progressCount - occurrences) < 100)
                cout << fmt::format("{}/{} : {}", (progressCount + occurrences), seqCount, offtargetsCount) << endl;
        }

        tempSeqSignatures.write(reinterpret_cast<const char *>(&signature), sizeof(signature));
        tempSeqSignaturesOccurrences.write(reinterpret_cast<const char *>(&occurrences), sizeof(occurrences));
        offtargetsCount++;
        if (progressCount % 10000 == 0)
            cout << fmt::format("{}/{} : {}", progressCount, seqCount, offtargetsCount) << endl;
        progressCount += occurrences;
    }
    entireDataSet.close();
    tempSeqSignatures.close();
    tempSeqSignaturesOccurrences.close();
    cout << "Finished!" << endl;

    boost::iostreams::mapped_file_source seqSignatures;
    boost::iostreams::mapped_file_source seqSignaturesOccurrences;
    seqSignatures.open(fmt::format("{}.tmp.1", argv[4]));
    seqSignaturesOccurrences.open(fmt::format("{}.tmp.2", argv[4]));
    uint64_t seqSignaturesCount = std::distance(reinterpret_cast<const uint64_t*>(seqSignatures.begin()), reinterpret_cast<const uint64_t*>(seqSignatures.end()));

    cout << "Writing index header to file..." << endl;
    ofstream isslIndex;
    isslIndex.open(argv[4], std::ios::out | std::ios::binary);
    isslIndex.write(reinterpret_cast<char*>(&offtargetsCount), sizeof(uint64_t));
    isslIndex.write(reinterpret_cast<char*>(&seqLength), sizeof(uint64_t));
    isslIndex.write(reinterpret_cast<char*>(&sliceCount), sizeof(size_t));
    cout << "Finished!" << endl;

    cout << "Writing offtargets to file..." << endl;
    isslIndex.write(seqSignatures.data(), seqSignatures.size());
    cout << "Finished!" << endl;

    cout << "Writing slice masks to file..." << endl;
    for (uint64_t& maskBinary : sliceMasksBinary)
    {
        isslIndex.write(reinterpret_cast<char*>(&maskBinary), sizeof(uint64_t));
    }
    cout << "Finished!" << endl;

    isslIndex.close();

    cout << "Constructing index..." << endl;
    for (size_t i = 0; i < sliceMasks.size(); i++)
    {
        cout << fmt::format("\tBuilding slice list {}", i+1) << endl;
        size_t sliceListSize = 1ULL << (sliceMasks[i].size() * 2);
        vector<vector<uint64_t>> sliceList(sliceListSize);
        for (uint32_t signatureId = 0; signatureId < seqSignaturesCount; signatureId++) {
            const uint64_t* signature = reinterpret_cast<const uint64_t*>(seqSignatures.data()) + signatureId;
            const uint32_t* occurrences = reinterpret_cast<const uint32_t*>(seqSignaturesOccurrences.data()) + signatureId;
            uint32_t sliceVal = 0ULL;
            for (size_t j = 0; j < sliceMasks[i].size(); j++)
            {
                sliceVal |= ((*signature >> (sliceMasks[i][j] * 2)) & 3ULL) << (j * 2);
            }
            // seqSigIdVal represnets the sequence signature ID and number of occurrences of the associated sequence.
            // (((uint64_t)occurrences) << 32), the most significant 32 bits is the count of the occurrences.
            // (uint64_t)signatureId, the index of the sequence in `seqSignatures`
            uint64_t seqSigIdVal = (static_cast<uint64_t>(*occurrences) << 32) | static_cast<uint64_t>(signatureId);
            sliceList[sliceVal].push_back(seqSigIdVal);
        }
        cout << "\tFinished!" << endl;

        cout << fmt::format("\tWriting slice list {} to file...", i+1) << endl;
        isslIndex.open(argv[4], std::ios::out | std::ios::binary | std::ios::app);
        // Write slice list lengths
        for (size_t j = 0; j < sliceListSize; j++) { // Slice limit given slice width
            size_t sz = sliceList[j].size();
            isslIndex.write(reinterpret_cast<char*>(&sz), sizeof(size_t));
        }
        // write slice list data
        for (size_t j = 0; j < sliceListSize; j++) { // Slice limit given slice width
            isslIndex.write(reinterpret_cast<char*>(sliceList[j].data()), sizeof(uint64_t) * sliceList[j].size());
        }
        isslIndex.close();
        cout << "\tFinished!" << endl;
    }
    seqSignatures.close();
    seqSignaturesOccurrences.close();
    remove(fmt::format("{}.tmp.1", argv[4]));
    remove(fmt::format("{}.tmp.2", argv[4]));

    cout << "Finished!" << endl;
    cout << "Done" << endl;
    return 0;
}
