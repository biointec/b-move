/******************************************************************************
 *  b-move: bidirectional move structure                                      *
 *  Copyright (C) 2020-2024 - Lore Depuydt <lore.depuydt@ugent.be> and        *
 *                            Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

/******************************************************************************
 *  This file includes code that was inspired by the                          *
 *  code from the repository: https://github.com/U-Ar/full-br-index           *
 *  by Yuma Arakawa.                                                          *
 *                                                                            *
 *  Citation of the br-index:                                                                 *
 *  Arakawa, Y., Navarro, G., & Sadakane, K. (2022). Bi-Directional r-Indexes.*
 *  In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022).    *
 *  Schloss Dagstuhl-Leibniz-Zentrum für Informatik.                          *
 ******************************************************************************/

using namespace std;

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "../external/Big-BWT/utils.h"
#include "alphabet.h"
#include "buildhelpers.h"
#include "move.h"
#include "moveElement.h"
#include "plcp.h"
#include "sparseBitvec.h"
#include "util.h"

using namespace std;

void generateBWT(const string& T, const vector<length_t>& SA, string& BWT) {
    // build the BWT
    std::cout << "Generating BWT..." << std::endl;
    BWT.resize(T.size());
    for (size_t i = 0; i < SA.size(); i++)
        if (SA[i] > 0)
            BWT[i] = T[SA[i] - 1];
        else
            BWT[i] = T.back();
}

/**
 * Get run index of a given position.
 * @param position The position for which to find the run index.
 * @param size The amount of input-output intervals.
 * @param dPair Vector with pairs of input-output interval start positions.
 * @returns The run index of the given position.
 */
length_t getRunIndex(const length_t position, const length_t size,
                     const vector<MoveElement>& rows) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the posible range smaller by binary search, untill only
    // 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the
        // test value.
        if (rows[testIndex].inputStart <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex - 1;
        }
    }

    assert(position >= rows[leftBoundary].inputStart);
    assert(leftBoundary == size - 1 ||
           position < rows[leftBoundary + 1].inputStart);

    return leftBoundary;
}

/**
 * Fill dIndex array.
 * @param tIn Balanced tree: stores pairs of I_out with increasing input
 * interval index.
 * @param rows The vector to fill, rows[i] = input-output interval i +
 * index j of input interval containing q_i.
 * @param size The amount of input-output intervals.
 * @param bwtSize The BWT size.
 */
void fillRows(const map<length_t, pair<uint8_t, length_t>>& tIn,
              vector<MoveElement>& rows, const length_t size,
              const length_t bwtSize) {
    rows.resize(size + 1);

    map<length_t, pair<uint8_t, length_t>>::const_iterator tInIt = tIn.begin();
    length_t i = 0;
    while (tInIt != tIn.end()) {
        length_t inputStart = tInIt->first;
        uint8_t c = tInIt->second.first;
        length_t outputStart = tInIt->second.second;

        rows[i] = MoveElement(c, inputStart, outputStart, 0);

        i++;
        tInIt++;
    }

    tInIt = tIn.begin();
    i = 0;
    while (tInIt != tIn.end()) {
        length_t outputRunIndex = getRunIndex(tInIt->second.second, size, rows);

        rows[i].mappingIndex = outputRunIndex;

        i++;
        tInIt++;
    }

    rows[size] = MoveElement(0, bwtSize, bwtSize, size);
}

/**
 * Creates the Move structure.
 * @param baseFN the base filename to write the move structure to.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 * @param charCounts Accumulated number of characters in lex order.
 * @param sigma The alphabet.
 */
length_t createAndWriteMove(const string& baseFN, const string& BWT,
                            const vector<length_t>& charCounts,
                            const Alphabet<ALPHABET>& sigma) {

    // Get the alphabet size
    size_t S = sigma.size();

    // Set bwt size
    length_t bwtSize = BWT.size();

    // balanced tree: stores pairs of I_out with increasing input interval index
    map<length_t, pair<uint8_t, length_t>> tIn;

    // Create accumulated charCounts
    vector<length_t> charCountsAcc(S, 0);
    length_t total = 0;
    for (size_t i = 0; i < S; i++) {
        char c = sigma.i2c(i);
        charCountsAcc[i] = total;
        total += charCounts[c];
    }

    // fill tIn, tOut
    vector<length_t> charsVisisted(S, 0);
    length_t prevC = S;
    length_t zeroCharPos;
    for (length_t i = 0; i < bwtSize; i++) {

        length_t c = sigma.c2i(BWT[i]);
        if (c == 0) {
            zeroCharPos = i;
        }
        if (prevC != c) {
            length_t lf = charCountsAcc[c] + charsVisisted[c];
            tIn[i] = make_pair(c, lf);
        }

        charsVisisted[c]++;
        prevC = c;
    }

    // Set dPair and dIndex size
    length_t arraySize = tIn.size();
    length_t size = arraySize;
    std::cout << "There are " + to_string(size) + " runs." << std::endl;
    // Fill dPair and dIndex
    vector<MoveElement> rows;
    fillRows(tIn, rows, size, BWT.size());

    // Write Move structures to files
    string fileName = baseFN;
    fileName += ".move";
    ofstream ofs(fileName, ios::binary);

    // Write bwtSize, amount of input-output intervals and the alphabet size to
    // the output stream.
    ofs.write((char*)&bwtSize, sizeof(bwtSize));
    ofs.write((char*)&size, sizeof(size));
    ofs.write((char*)&zeroCharPos, sizeof(zeroCharPos));

    // Write rows to the output stream.
    for (length_t i = 0; i < size + 1; i++) {
        rows[i].serialize(ofs);
    }

    ofs.close();
    std::cout << "Wrote file " + fileName << std::endl;
    return size;
}

/**
 * Print usage to stdout
 */
void showUsage() {
    cout << "Usage: ./bmove-build [--pfp] [--preprocess] <fasta file>\n\n";
    cout << "Required arguments:\n";
    cout << "\t fasta file: Path to the FASTA file containing the reference "
            "sequences.\n\n";

    cout << "Optional arguments:\n";
    cout << "\t--preprocess:                   Preprocess the fasta file only. "
            "Do not build the "
            "index.\n";
    cout << "\t--pfp:                          Start the index construction "
            "after the prefix-free "
            "parsing step.\n\n";

    cout << "Report bugs to lore.depuydt@ugent.be" << endl;
}

bool parseArguments(int argc, char* argv[], string& baseFN,
                    bool& preprocessOnly, bool& pfp) {
    if (argc < 2)
        return false;

    // Set defaults
    preprocessOnly = false;
    pfp = false;

    // Process flags and options
    for (int i = 1; i < argc - 1; ++i) {
        string arg(argv[i]);
        if (arg == "--preprocess") {
            preprocessOnly = true;
        } else if (arg == "--pfp") {
            pfp = true;
        } else {
            std::cerr << "Unknown argument or missing value for build argument "
                         "flag: " +
                             arg
                      << std::endl;
            return false;
        }
    }

    // The last argument is the base filename
    baseFN = argv[argc - 1];

    return true;
}

void writeIntVectorBinary(const std::string& filename,
                          const std::vector<length_t>& array) {
    // convert to int_vector
    uint8_t width =
        (uint8_t)ceil(log2(*max_element(array.begin(), array.end())));
    sdsl::int_vector<> intVector(array.size(), 0, width);
    for (size_t i = 0; i < array.size(); i++) {
        intVector[i] = array[i];
    }
    std::ofstream ofs(filename, std::ios::binary);
    intVector.serialize(ofs);
    ofs.close();
}

/**
 * Builds samplesFirst and samplesLast.
 * @param samplesFirst Vector to fill. First sa samples of each input interval.
 * @param samplesLast Vector to fill. Last sa samples of each input interval.
 * @param SA The suffix array.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 */
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast,
                  const vector<length_t>& SA, const string& BWT) {

    samplesFirst.push_back(SA[0] > 0 ? SA[0] - 1 : BWT.size() - 1);
    for (size_t pos = 0; pos < BWT.size() - 1; pos++) {
        if (BWT[pos] != BWT[pos + 1]) {
            samplesLast.push_back(SA[pos] > 0 ? SA[pos] - 1 : BWT.size() - 1);
            samplesFirst.push_back(SA[pos + 1] > 0 ? SA[pos + 1] - 1
                                                   : BWT.size() - 1);
        }
    }
    samplesLast.push_back(SA[BWT.size() - 1] > 0 ? SA[BWT.size() - 1] - 1
                                                 : BWT.size() - 1);
}

/**
 * Generate predecessor structures
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 * @param verbose print verbose output to stdout
 */
void generatePredecessors(const vector<length_t>& samplesFirst,
                          const vector<length_t>& samplesLast,
                          SparseBitvec& predFirst, SparseBitvec& predLast,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength, false);
    vector<bool> predLastBV(textLength, false);

    for (length_t SApos : samplesFirst) {
        predFirstBV[SApos] = true;
    }

    for (length_t SApos : samplesLast) {
        predLastBV[SApos] = true;
    }

    predFirst = SparseBitvec(predFirstBV);
    predLast = SparseBitvec(predLastBV);
}

/**
 * Generate predToRun arrays
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] firstToRun mapping between rank of ones in predFirst bitvector
 * and run indices
 * @param [out] lastToRun mapping between rank of ones in predLast bitvector and
 * run indices
 */
void generatePredTorun(const vector<length_t>& samplesFirst,
                       const vector<length_t>& samplesLast,
                       vector<length_t>& firstToRun,
                       vector<length_t>& lastToRun) {

    firstToRun.resize(samplesFirst.size());
    iota(firstToRun.begin(), firstToRun.end(), 0);
    sort(firstToRun.begin(), firstToRun.end(),
         [&samplesFirst](length_t a, length_t b) {
             return samplesFirst[a] < samplesFirst[b];
         });

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(),
         [&samplesLast](length_t a, length_t b) {
             return samplesLast[a] < samplesLast[b];
         });
}

void readSuffixArrayFile(const std::string& baseFN,
                         const std::string& extension,
                         vector<length_t>& samples, vector<length_t>& toRun,
                         length_t size, length_t nrOfRuns, bool reverse) {
    std::cout << "Reading suffix array samples from " + baseFN + extension +
                     "..."
              << std::endl;

    string fileName = baseFN + extension;

    // Open the file
    FILE* file = fopen(fileName.c_str(), "rb");

    samples = vector<length_t>(nrOfRuns, 0);
    toRun = vector<length_t>(nrOfRuns, 0);
    std::vector<std::pair<length_t, length_t>> pos_run_pairs(nrOfRuns);

    length_t* buf = new length_t[1];

    for (length_t i = 0; i < nrOfRuns; ++i) {
        // Read the first SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Read the next SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Calculate sa_val from buf: take the first byte, ensure it's within
        // range, and adjust as needed based on the data size
        length_t sa_val = buf[0] % (1UL << (8 * SABYTES));
        if (!reverse) {
            // Only for the regular suffix array, sa_val must be corrected.
            // For the reverse case, the value is already correct since Big-BWT
            // puts the sentinel character at the end of the reverse text as
            // well.
            sa_val = (sa_val > 0) ? (sa_val - 1) : (size - 1);
        }

        // Store sa_val in samples array at index i
        samples[i] = sa_val;

        // Store {sa_val, i} pair in pos_run_pairs array
        pos_run_pairs[i] = {sa_val, i};
    }

    delete[] buf;

    if (!reverse) {

        std::cout
            << "Creating the mapping between the predecessor bits and the "
               "runs..."
            << std::endl;

        std::sort(pos_run_pairs.begin(), pos_run_pairs.end());

        std::vector<length_t> positions;
        for (length_t i = 0; i < nrOfRuns; ++i) {
            positions.push_back(pos_run_pairs[i].first);
            toRun[i] = pos_run_pairs[i].second;
        }
    }

    fclose(file);
}

void constructRunLengthEncodedPLCP(const SparseBitvec& first,
                                   const vector<length_t>& first_to_run,
                                   const std::vector<length_t>& charCounts,
                                   length_t size, length_t r, PLCP& plcp,
                                   const std::string& baseFN,
                                   const Alphabet<ALPHABET>& sigma) {
    std::cout << "Constructing run-length encoded PLCP..." << std::endl;
    // Construct cumulative charcounts
    vector<length_t> charCountsCumulative(256, 0);
    for (size_t i = 1; i < 256; i++) {
        auto temp = charCountsCumulative[i - 1] + charCounts[i];
        charCountsCumulative[i] = temp;
    }

    std::cout << "Reading " + baseFN + ".move..." << std::endl;

    Move<ALPHABET> rows;
    rows.load(baseFN, false);

    std::cout << "Creating bit vectors to support select for each "
                 "character in the BWT..."
              << std::endl;

    sdsl::bit_vector char_bv[ALPHABET];
    sdsl::bit_vector::select_1_type select[ALPHABET];

    {

        for (length_t i = 0; i < ALPHABET; i++) {
            char_bv[i] = sdsl::bit_vector(size, 0);
        }

        for (length_t i = 0; i < r; i++) {
            size_t c = rows.getRunHead(i);
            for (length_t j = rows.getInputStartPos(i);
                 j < rows.getInputStartPos(i + 1); j++) {
                char_bv[c][j] = true;
            }
        }

        for (length_t i = 0; i < ALPHABET; i++) {
            select[i] = sdsl::bit_vector::select_1_type(&char_bv[i]);
        }
    }

    std::cout << "Creating PLCP..." << std::endl;

    std::vector<length_t> ones, zeros;

    char c, c0;
    length_t p, p0, pos, l, gap;
    length_t prev_pos = 0;
    length_t prev_l = 0;
    length_t acc0 = 0, acc1 = 0;

    try {

        for (length_t i = 0; i < r; ++i) {
            // Log a progress message that takes r into account (e.g 100 times
            // from 0 to r-1)
            if (i % (r / 100) == 0) {
                cout << "Progress: " + std::to_string(i) + " / " +
                            std::to_string(r) + " (" +
                            std::to_string((i * 100) / r) + "%)" << "\r";
                cout.flush();
            }

            try {
                pos = first.select(i);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << '\n';
                std::cerr << "While trying to perform select on first, with i: "
                          << i << '\n';
                exit(1);
            }

            gap = i == 0 ? pos : (pos - prev_pos - 1);

            l = 0;

            length_t tempIdx = first_to_run[i];
            p = rows.getInputStartPos(tempIdx);
            c = sigma.i2c(rows.getRunHead(tempIdx));
            rows.findLF(p, tempIdx);
            if (p != 0) {
                p0 = p - 1;
                c0 = (std::upper_bound(charCountsCumulative.begin(),
                                       charCountsCumulative.end(), p0) -
                      charCountsCumulative.begin());
                while (c == c0) {
                    l++;
                    length_t i = p - charCountsCumulative[c - 1];
                    try {
                        p = select[sigma.c2i(c)](i + 1);
                    } catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << '\n';
                        std::cerr << "While trying to perform select on "
                                     "sigma.c2i(c), with i: "
                                  << i << " and c: " << c << '\n';
                        exit(1);
                    }
                    c = (std::upper_bound(charCountsCumulative.begin(),
                                          charCountsCumulative.end(), p) -
                         charCountsCumulative.begin());
                    length_t i0 = p0 - charCountsCumulative[c0 - 1];
                    try {
                        p0 = select[sigma.c2i(c0)](i0 + 1);
                    } catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << '\n';
                        std::cerr << "While trying to perform select on "
                                     "sigma.c2i(c0), with i0: "
                                  << i0 << " and c0: " << c0 << '\n';
                        exit(1);
                    }
                    c0 = (std::upper_bound(charCountsCumulative.begin(),
                                           charCountsCumulative.end(), p0) -
                          charCountsCumulative.begin());
                }
            }

            if (i == 0) {
                if (l + gap > 0) {
                    ones.push_back(acc1);
                    acc0 += l + gap - 1;
                    zeros.push_back(acc0);
                    acc1 += gap + 1;
                } else {
                    acc1 += gap + 1;
                }
            } else if (l + gap + 1 - prev_l) {
                ones.push_back(acc1);
                acc0 += l + gap + 1 - prev_l;
                zeros.push_back(acc0);
                acc1 += gap + 1;
            } else {
                acc1 += gap + 1;
            }

            prev_pos = pos;
            prev_l = l;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "While trying to construct PLCP" << '\n';
        exit(1);
    }

    ones.push_back(acc1);

    std::cout << "Copying the PCLP into its correct object..." << std::endl;

    plcp = PLCP(size, ones, zeros);

    std::cout << "Constructed run-length encoded PLCP." << std::endl;
}

void createIndex(const string& baseFN, const string& T, bool fromFasta) {
    // count the frequency of each characters in T
    vector<length_t> charCounts;
    countChars(T, charCounts);
    // write the character counts table
    writeCharCounts(baseFN, charCounts);

    // Create the alphabet
    Alphabet<ALPHABET> sigma;
    createAlphabet(T, charCounts, sigma);

    {
        // read the suffix array
        vector<length_t> SA;
        if (fromFasta) {
            // create suffix array from concatenated using radixSA64
            createSuffixArray(T, SA);
        } else {
            std::cout << "Reading " + baseFN + ".sa..." << std::endl;
            readSA(baseFN + ".sa", SA, T.size());
        }

        // perform a sanity check on the suffix array
        std::cout << "Performing sanity checks..." << std::endl;
        sanityCheck(T, SA);
        std::cout << "\tSanity checks OK" << std::endl;
        // build the BWT
        string BWT;
        generateBWT(T, SA, BWT);

        { // generate and write PLCP array
            std::cout
                << "Generating the permuted longest common prefix array..."
                << std::endl;
            PLCP plcp(T, SA);
            plcp.write(baseFN + ".plcp");
            std::cout << "Wrote file " + baseFN + ".plcp" << std::endl;
        }

        // Create samplesLast and samplesFirst
        std::cout << "Sampling the suffix array values at run boundaries..."
                  << std::endl;
        vector<length_t> samplesFirst;
        vector<length_t> samplesLast;
        buildSamples(samplesFirst, samplesLast, SA, BWT);

        // Clear the suffix array
        SA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        std::cout << "Wrote file " + baseFN + ".smpf" << std::endl;
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        std::cout << "Wrote file " + baseFN + ".smpl" << std::endl;
        SparseBitvec predFirst;
        SparseBitvec predLast;
        vector<length_t> firstToRun;
        vector<length_t> lastToRun;

        // Generate the predecessor bitvectors
        std::cout << "Generating the predecessor bitvectors for the samples..."
                  << std::endl;
        generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                             BWT.size());
        // Generate the predToRun arrays
        std::cout
            << "Mapping the predecessor bits to their corresponding runs..."
            << std::endl;
        generatePredTorun(samplesFirst, samplesLast, firstToRun, lastToRun);

        samplesFirst.clear();
        samplesLast.clear();

        // Write predecessor structures
        predFirst.write(baseFN + ".prdf");
        std::cout << "Wrote file " + baseFN + ".prdf" << std::endl;
        predLast.write(baseFN + ".prdl");
        std::cout << "Wrote file " + baseFN + ".prdl" << std::endl;
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        std::cout << "Wrote file " + baseFN + ".ftr" << std::endl;
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        std::cout << "Wrote file " + baseFN + ".ltr" << std::endl;
        firstToRun.clear();
        lastToRun.clear();

        // Create the Move structure
        std::cout << "Creating the move table..." << std::endl;
        createAndWriteMove(baseFN, BWT, charCounts, sigma);
        BWT.clear();
    }

    std::cout << "Switching to reversed text..." << std::endl;
    { // read or create the reverse suffix array
        vector<length_t> revSA;
        if (fromFasta) {
            std::string revT = T;
            std::reverse(revT.begin(), revT.end());

            // create the suffix array of the reverse text
            createSuffixArray(revT, revSA);
            revT.clear();

        } else {
            std::cout << "Reading " + baseFN + ".rev.sa..." << std::endl;
            readSA(baseFN + ".rev.sa", revSA, T.size());
        }

        // perform a sanity check on the suffix array
        std::cout << "Performing sanity checks..." << std::endl;
        sanityCheck(T, revSA);
        std::cout << "\tSanity checks OK" << std::endl;
        // build the reverse BWT
        string revBWT;
        createRevBWT(baseFN, T, revSA, sigma, revBWT);

        // Create samplesLast and samplesFirst
        std::cout
            << "Sampling the reverse suffix array values at run boundaries..."
            << std::endl;
        vector<length_t> revSamplesFirst;
        vector<length_t> revSamplesLast;
        buildSamples(revSamplesFirst, revSamplesLast, revSA, revBWT);

        // Clear the reverse suffix array
        revSA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        std::cout << "Wrote file " + baseFN + ".rev.smpf" << std::endl;
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        std::cout << "Wrote file " + baseFN + ".rev.smpl" << std::endl;
        revSamplesFirst.clear();
        revSamplesLast.clear();

        // Create the Move structure
        std::cout << "Creating the reverse move table..." << std::endl;
        createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();
    }

    // Write the BMOVE_BUILD_INDEX_TAG to a tag file
    ofstream tagFile(baseFN + ".tag");
    tagFile << BMOVE_BUILD_INDEX_TAG;
    tagFile.close();
    // Write the 64 or 32-bit compiled info to a file
    ofstream compiledInfoFile(baseFN + ".comp");
    compiledInfoFile << sizeof(length_t);
    compiledInfoFile.close();
}

void createIndexPFP(const string& baseFN) {

    // build the BWT
    string BWT;
    // Read the BWT from disk
    std::cout << "Reading " + baseFN + ".bwt..." << std::endl;
    readText(baseFN + ".bwt", BWT);

    // Replace '\0' with '$'
    std::replace(BWT.begin(), BWT.end(), '\0', '$');

    length_t bwtSize = BWT.size();

    // count the frequency of each characters in T
    vector<length_t> charCounts;
    countChars(BWT, charCounts);

    // write the character counts table
    writeCharCounts(baseFN, charCounts);

    // Create the alphabet
    Alphabet<ALPHABET> sigma;
    createAlphabet(BWT, charCounts, sigma);

    { // Create the Move structure
        std::cout << "Creating the move table..." << std::endl;
        length_t nrOfRuns = createAndWriteMove(baseFN, BWT, charCounts, sigma);

        BWT.clear();

        vector<length_t> samplesFirst;
        vector<length_t> firstToRun;
        readSuffixArrayFile(baseFN, ".ssa", samplesFirst, firstToRun, bwtSize,
                            nrOfRuns, false);

        vector<length_t> samplesLast;
        vector<length_t> lastToRun;
        readSuffixArrayFile(baseFN, ".esa", samplesLast, lastToRun, bwtSize,
                            nrOfRuns, false);

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        std::cout << "Wrote file " + baseFN + ".smpf" << std::endl;
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        std::cout << "Wrote file " + baseFN + ".smpl" << std::endl;
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        std::cout << "Wrote file " + baseFN + ".ftr" << std::endl;
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        std::cout << "Wrote file " + baseFN + ".ltr" << std::endl;
        lastToRun.clear();

        SparseBitvec predFirst;
        SparseBitvec predLast;

        // Generate the predecessor bitvectors
        std::cout << "Generating the predecessor bitvectors for the samples..."
                  << std::endl;
        generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                             bwtSize);

        samplesFirst.clear();
        samplesLast.clear();

        // Write predecessor structures
        predFirst.write(baseFN + ".prdf");
        std::cout << "Wrote file " + baseFN + ".prdf" << std::endl;
        predLast.write(baseFN + ".prdl");
        std::cout << "Wrote file " + baseFN + ".prdl" << std::endl;
        // generate and write PLCP array
        std::cout << "Generating the permuted longest common prefix array..."
                  << std::endl;
        PLCP plcp;
        constructRunLengthEncodedPLCP(predFirst, firstToRun, charCounts,
                                      bwtSize, nrOfRuns, plcp, baseFN, sigma);
        plcp.write(baseFN + ".plcp");
        std::cout << "Wrote file " + baseFN + ".plcp" << std::endl;
        firstToRun.clear();
    }

    std::cout << "Switching to reversed text..." << std::endl;
    {
        // build the BWT
        string revBWT;
        // Read the BWT from disk
        std::cout << "Reading " + baseFN + ".rev.bwt..." << std::endl;
        readText(baseFN + ".rev.bwt", revBWT);

        // Replace '\0' with '$'
        std::replace(revBWT.begin(), revBWT.end(), '\0', '$');

        // Create the Move structure
        std::cout << "Creating the move table..." << std::endl;
        length_t nrOfRuns =
            createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();

        vector<length_t> revSamplesFirst;
        vector<length_t> revFirstToRun;
        readSuffixArrayFile(baseFN, ".rev.ssa", revSamplesFirst, revFirstToRun,
                            bwtSize, nrOfRuns, true);

        revFirstToRun.clear();

        vector<length_t> revSamplesLast;
        vector<length_t> revLastToRun;
        readSuffixArrayFile(baseFN, ".rev.esa", revSamplesLast, revLastToRun,
                            bwtSize, nrOfRuns, true);
        revLastToRun.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        std::cout << "Wrote file " + baseFN + ".rev.smpf" << std::endl;
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        std::cout << "Wrote file " + baseFN + ".rev.smpl" << std::endl;
        revSamplesFirst.clear();
        revSamplesLast.clear();
    }

    // Write the BMOVE_BUILD_INDEX_TAG to a tag file
    ofstream tagFile(baseFN + ".tag");
    tagFile << BMOVE_BUILD_INDEX_TAG;
    tagFile.close();
    // Write the 64 or 32-bit compiled info to a file
    ofstream compiledInfoFile(baseFN + ".comp");
    compiledInfoFile << sizeof(length_t);
    compiledInfoFile.close();
}

void createBMove(const string& baseFN) {

    // read the text file from disk
    std::cout << "Reading " + baseFN + ".txt..." << std::endl;
    string T;
    readText(baseFN + ".txt", T);
    if (!T.empty() && T.back() != '$') {
        T += '$';
    }
    checkTextSize(T.size());

    std::cout << "Converting text to uppercase..." << std::endl;
    std::transform(T.begin(), T.end(), T.begin(), ::toupper);

    // create the sequence names vector, positions vector and the header
    vector<length_t> positions = {0, (length_t)T.size()};
    vector<string> seqNames = {baseFN};

    // create the header lines with these reference sequences
    createAndWriteHeaderInfo(baseFN, seqNames, positions);

    // Write the positions to disk
    writePositionsAndSequenceNames(baseFN, positions, seqNames);

    createIndex(baseFN, T, false);
}

void processFastaFile(const std::string& fastaFile, const std::string& baseFN) {
    std::string T; // the (concatenated) text
    preprocessFastaFile(fastaFile, baseFN, T);

    createIndex(baseFN, T, true);
}

int main(int argc, char* argv[]) {

    std::array<string, 4> allowedExtensionsFasta = {".fasta", ".fa", ".FASTA",
                                                    ".FA"};
    std::array<string, 4> requiredExtensionsTxt = {".txt", ".sa", ".rev.sa",
                                                   ".rev.txt"};
    std::array<string, 6> requiredExtensionsPFP = {
        ".bwt", ".esa", ".ssa", ".rev.bwt", ".rev.esa", ".rev.ssa"};

    string baseFN;
    bool preprocessOnly;
    bool pfp;

    if (!parseArguments(argc, argv, baseFN, preprocessOnly, pfp)) {
        showUsage();
        return EXIT_FAILURE;
    }
    std::cout << "Welcome to the bidirectional move structure construction!"
              << std::endl;
    std::cout << "Alphabet size is " + std::to_string(ALPHABET - 1) + " + 1"
              << std::endl;
    // check which mode to use
    bool txtMode = true;
    bool fastaMode = true;

    for (auto& ext : requiredExtensionsTxt) {
        string filename = baseFN + ext;

        // check if file with filename exists
        ifstream ifs(filename);
        if (!ifs) {
            // file does not exist
            txtMode = false;
            break;
        }
    }

    // to do set fastaMode to true if any of the allowed extensions
    for (auto& ext : allowedExtensionsFasta) {
        string filename = baseFN + ext;
        // check if file with filename exists
        ifstream ifs(filename);
        if (ifs) {
            // file does exist -> use fasta mode
            fastaMode = true;
            break;
        }
    }

    // check if the base file already ends with fasta extension
    if (baseFN.size() >= 6) {
        string posExt = baseFN.substr(baseFN.size() - 6, 6);
        for (char& c : posExt) {
            c = std::tolower(c);
        }
        if (posExt == ".fasta") {
            fastaMode = true;
            baseFN = baseFN.substr(0, baseFN.size() - 6);
        }
    }

    if (preprocessOnly) {
        std::cout << "Preprocessing only..." << std::endl;
        if (!fastaMode) {
            std::cerr << "Preprocessing only works with fasta files."
                      << std::endl;
            return EXIT_FAILURE;
        }
        std::string T; // the (concatenated) text
        preprocessFastaFile(baseFN + ".fasta", baseFN, T, true);

        // Remove the sentinel character
        T.pop_back();

        // Write the preprocessed text to disk
        std::cout << "Writing concatenated uppercase sequence to disk..."
                  << std::endl;
        std::ofstream ofs(baseFN);
        ofs << T;
        ofs.close();

        // Reverse the text
        std::cout << "Reversing text..." << std::endl;
        std::reverse(T.begin(), T.end());

        // Write the reversed text to disk
        std::cout << "Writing reversed text to disk..." << std::endl;
        ofs = std::ofstream(baseFN + ".rev");
        ofs << T;
        ofs.close();

        std::cout << "Fasta input preprocessed successfully!" << std::endl;
        std::cout << "Exiting... bye!" << std::endl;
        return EXIT_SUCCESS;
    }

    if (pfp) {
        std::cout
            << "Starting index construction after prefix-free parsing step..."
            << std::endl;

        for (auto& ext : requiredExtensionsPFP) {
            string filename = baseFN + ext;

            // check if file with filename exists
            ifstream ifs(filename);
            if (!ifs) {
                std::cerr << "Missing file: " + filename << std::endl;
                std::cerr << "Please run the prefix-free parsing step first."
                          << std::endl;
                return EXIT_FAILURE;
            }
        }

        createIndexPFP(baseFN);

        std::cout << "Index construction completed successfully!" << std::endl;
        std::cout << "Exiting... bye!" << std::endl;
        return EXIT_SUCCESS;
    }

    if (txtMode) {
        // txt takes precedence as it has already preprocessed suffix arrays
        fastaMode = false;
    }

    try {
        if (txtMode) {
            createBMove(baseFN);
        } else if (fastaMode) {
            processFastaFile(baseFN + ".fasta", baseFN);
        } else {
            throw runtime_error("No valid input files found");
        }
    } catch (const std::exception& e) {
        std::cerr << "Fatal: " + string(e.what()) << std::endl;

        return EXIT_FAILURE;
    }

    std::cout << "Index construction completed successfully!" << std::endl;
    std::cout << "Exiting... bye!" << std::endl;
    return EXIT_SUCCESS;
}
