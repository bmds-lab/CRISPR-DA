#ifndef FMT_HEADER_ONLY
#   define FMT_HEADER_ONLY
#endif
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <boost/iostreams/device/mapped_file.hpp>

#include "libpopcnt/libpopcnt.h"
#include "penalties.hpp"

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

enum otScoreMethod { mit = 0, cfd = 1, mitAndCfd = 2, mitOrCfd = 3, avgMitCfd = 4 };

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

uint64_t encode2bit(string seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[seq[j] - nucleotideIndexOffset]) << (j * 2);
    }
    return signature;
}

string decode2bit(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method]\n", argv[0]);
        exit(1);
    }

    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);

    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);

    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    string argScoreMethod = argv[5];
    otScoreMethod scoreMethod;
    bool calcCfd = false;
    bool calcMit = false;
    if (!argScoreMethod.compare("and")) {
        scoreMethod = otScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("or")) {
        scoreMethod = otScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("avg")) {
        scoreMethod = otScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("mit")) {
        scoreMethod = otScoreMethod::mit;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("cfd")) {
        scoreMethod = otScoreMethod::cfd;
        calcCfd = true;
    }
    else
    {
        fprintf(stderr, "Invalid scoring method. Acceptable options are: 'and', 'or', 'avg', 'mit', 'cfd'");
        exit(1);
    }

    /** Begin reading the binary encoded ISSL, structured as:
     *  - The header (3 items)
     *  - All binary-encoded off-target sites
     *  - Slice masks
     *  - Size of slice 1 lists
     *  - Contents of slice 1 lists
     *  ...
     *  - Size of slice N lists (N being the number of slices)
     *  - Contents of slice N lists
     */
    mapped_file_source isslFp;
    isslFp.open(argv[1]);

    if (!isslFp.is_open())
    {
        throw std::runtime_error("Error reading index: could not open file\n");
    }
    const void* inFileFp = static_cast<const void*>(isslFp.data());

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the number of slices
     */
    const size_t* headerPtr = static_cast<const size_t*>(inFileFp);
    size_t offtargetsCount = *headerPtr++;
    size_t seqLength = *headerPtr++;
    size_t sliceCount = *headerPtr++;

    /** Load in all of the off-target sites */
    const uint64_t* offtargetsPtr = static_cast<const uint64_t*>(headerPtr);

    /** Read the slice masks and generate 2 bit masks */
    const uint64_t* sliceMasksPtr = static_cast<const uint64_t*>(offtargetsPtr + offtargetsCount);
    vector<vector<uint64_t>> sliceMasks;
    for (size_t i = 0; i < sliceCount; i++)
    {
        vector<uint64_t> mask;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (*sliceMasksPtr & (1ULL << j))
            {
                mask.push_back(j);
            }
        }
        sliceMasks.push_back(mask);
        sliceMasksPtr++;
    }

    /** The contents of the slices. Stored by slice
    * Contains:
    *   - Size of each list within the slice stored contiguously
    *   - The contents of all the lists stored contiguously
    */
    vector<const size_t*> allSlicelistSizes(sliceCount);
    vector<const uint64_t*> allSliceSignatures(sliceCount);
    const size_t* listSizePtr = static_cast<const size_t*>(sliceMasksPtr);
    const uint64_t* signaturePtr = static_cast<const uint64_t*>(sliceMasksPtr);
    for (size_t i = 0; i < sliceCount; i++)
    {
        allSlicelistSizes[i] = listSizePtr;
        signaturePtr = static_cast<const uint64_t*>(listSizePtr + (1ULL << (sliceMasks[i].size() * 2)));
        allSliceSignatures[i] = signaturePtr;
        listSizePtr = static_cast<const size_t*>(signaturePtr + offtargetsCount);
    }

    /** Prevent assessing an off-target site for multiple slices
     *
     *      Create enough 1-bit "seen" flags for the off-targets
     *      We only want to score a candidate guide against an off-target once.
     *      The least-significant bit represents the first off-target
     *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
     *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
     */
    uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

    /** Start constructing index in memory
     *
     *      To begin, reverse the contiguous storage of the slices,
     *         into the following:
     *
     *         + Slice 0 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |----  ...
     *         |
     *         + Slice 1 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |---- ...
     *         | ...
     */

    vector<vector<const uint64_t*>> sliceLists(sliceCount);
    // Assign sliceLists size based on each slice length
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<const uint64_t*>(1ULL << (sliceMasks[i].size() * 2));
    }

    for (size_t i = 0; i < sliceCount; i++) {
        const uint64_t* sliceList = allSliceSignatures[i];
        size_t sliceLimit = 1ULL << (sliceMasks[i].size() * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            sliceLists[i][j] = sliceList;
            sliceList += allSlicelistSizes[i][j];
        }
    }

    //TODO: rewrite
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    std::filesystem::path queryFile(argv[2]);
    size_t fileSize = std::filesystem::file_size(queryFile);
    size_t seqLineLength = seqLength + getLineEnding(queryFile, seqLength).size();
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    FILE* fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    vector<vector<uint64_t>> queryOfftargets(queryCount);
    vector<vector<uint64_t>> queryCounts(queryCount);
    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

    /** Binary encode query sequences */
    #pragma omp parallel
    {
    #pragma omp for
        for (int i = 0; i < queryCount; i++) {
            char* ptr = &queryDataSet[i * seqLineLength];
            uint64_t signature = encode2bit(ptr, 20);
            querySignatures[i] = signature;
        }
    }

    /** Begin scoring */
    #pragma omp parallel
    {
        vector<uint64_t> offtargetToggles(numOfftargetToggles);
        uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;
        /** For each candidate guide */
        // TODO: update to openMP > v2 (Use clang compiler)
        #pragma omp for
        for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            auto searchSignature = querySignatures[searchIdx];
            size_t sliceLimitOffset = 0;
            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                vector<uint64_t>& sliceMask = sliceMasks[i];
                auto& sliceList = sliceLists[i];

                uint64_t searchSlice = 0ULL;
                for (int j = 0; j < sliceMask.size(); j++)
                {
                    searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                }

                size_t idx = sliceLimitOffset + searchSlice;

                size_t signaturesInSlice = allSlicelistSizes[i][searchSlice];
                const uint64_t* sliceOffset = sliceList[searchSlice];

                /** For each off-target signature in slice */
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    auto signatureWithOccurrencesAndId = sliceOffset[j];
                    auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
                    uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

                    /** Prevent assessing the same off-target for multiple slices */
                    uint64_t seenOfftargetAlready = 0;
                    uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
                    seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;

                    if (!seenOfftargetAlready) {
                        *ptrOfftargetFlag |= (1ULL << (signatureId % 64));

                        /** Find the positions of mismatches
                            *
                            *  Search signature (SS):    A  A  T  T    G  C  A  T
                            *                           00 00 11 11   10 01 00 11
                            *
                            *        Off-target (OT):    A  T  A  T    C  G  A  T
                            *                           00 11 00 11   01 10 00 11
                            *
                            *                SS ^ OT:   00 00 11 11   10 01 00 11
                            *                         ^ 00 11 00 11   01 10 00 11
                            *                  (XORd) = 00 11 11 00   11 11 00 00
                            *
                            *        XORd & evenBits:   00 11 11 00   11 11 00 00
                            *                         & 10 10 10 10   10 10 10 10
                            *                   (eX)  = 00 10 10 00   10 10 00 00
                            *
                            *         XORd & oddBits:   00 11 11 00   11 11 00 00
                            *                         & 01 01 01 01   01 01 01 01
                            *                   (oX)  = 00 01 01 00   01 01 00 00
                            *
                            *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
                            *                         | 00 01 01 00   01 01 00 00
                            *            mismatches   = 00 01 01 00   01 01 00 00
                            *
                            *   popcount(mismatches):   4
                            */
                        uint64_t xoredSignatures = searchSignature ^ offtargetsPtr[signatureId];
                        uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                        uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                        uint64_t mismatches = (evenBits >> 1) | oddBits;
                        uint64_t dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= 4) {
                            queryOfftargets[searchIdx].push_back(offtargetsPtr[signatureId]);
                            queryCounts[searchIdx].push_back(occurrences);
							*ptrOfftargetFlag |= (1ULL << (signatureId % 64));
                        }
                    }
                }
            }
            memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
        }
    }

    /** Print global scores to stdout */
    for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
        auto querySequence = decode2bit(querySignatures[searchIdx], seqLength);
        printf("%s:", querySequence.c_str());
        for (uint64_t offTargetidx = 0; offTargetidx < queryOfftargets[searchIdx].size(); ++offTargetidx)
        {
            printf("%ld-%s,", queryCounts[searchIdx][offTargetidx], decode2bit(queryOfftargets[searchIdx][offTargetidx], seqLength).c_str());
        }
        printf("\n");
    }

    return 0;
}
