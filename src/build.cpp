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
 *  This file includes code inspired by the full-br-index repository          *
 *  Original source: https://github.com/U-Ar/full-br-index                    *
 *  Author: Yuma Arakawa                                                      *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Affero General Public License  *
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include <cstdint> // needed for Big-BWT/utils.h
#include <cstdio>  // needed for Big-BWT/utils.h

#include "../external/Big-BWT/utils.h" // for SABYTES

#include "alphabet.h"     // for Alphabet
#include "buildhelpers.h" // for preproc...
#include "definitions.h"  // for length_t
#include "logger.h"       // for Logger
#include "move.h"      // for MoveLFRepr
#include "movephirepr.h"     // for MovePhiRepr
#include "moveElement.h"         // for MoveRow
#include "plcp.h"            // for PLCP
#include "sparseBitvec.h"    // for SparseB...

#include <algorithm>                   // for sort
#include <array>                       // for array
#include <assert.h>                    // for assert
#include <cmath>                       // for ceil, log2
#include <exception>                   // for exception
#include <ext/alloc_traits.h>          // for __alloc...
#include <iostream>                    // for operator<<
#include <map>                         // for map
#include <memory>                      // for allocat...
#include <numeric>                     // for iota
#include <sdsl/int_vector.hpp>         // for bit_vector
#include <sdsl/select_support_mcl.hpp> // for select_...
#include <stdexcept>                   // for runtime...
#include <stdio.h>                     // for fread
#include <stdlib.h>                    // for exit
#include <string>                      // for operator+
#include <type_traits>                 // for __strip...
#include <utility>                     // for pair
#include <vector>                      // for vector

using namespace std;

/**
 * Get run index of a given position.
 * @param position The position for which to find the run index.
 * @param size The amount of input-output intervals.
 * @param dPair Vector with pairs of input-output interval start positions.
 * @returns The run index of the given position.
 */
length_t getRunIndex(const length_t position, const length_t size,
                     const MoveLFReprBP& rows) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the possible range smaller by binary search, until only
    // 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the
        // test value.
        if (rows.getInputStartPos(testIndex) <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex - 1;
        }
    }

    assert(position >= rows.getInputStartPos(leftBoundary));
    assert(leftBoundary == size - 1 ||
           position < rows.getInputStartPos(leftBoundary + 1));

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
              MoveLFReprBP& rows, const length_t size, const length_t bwtSize) {
    rows.initialize(tIn.size(), bwtSize);

    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->second.first, it->first, it->second.second, 0);
        assert(rows.getRunHead(i) == it->second.first);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second.second);
        i++;
    }

    i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        length_t outputRunIndex = getRunIndex(it->second.second, size, rows);
        rows.setOutputStartRun(i, outputRunIndex);
        assert(rows.getOutputStartRun(i) == outputRunIndex);
        i++;
    }

    rows.setRowValues(size, 0, bwtSize, bwtSize, size);
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
    vector<length_t> charsVisited(S, 0);
    length_t prevC = S;
    length_t zeroCharPos = bwtSize;
    for (length_t i = 0; i < bwtSize; i++) {

        char _c = BWT[i];
        length_t c = sigma.c2i(_c);
        assert(c < S);
        if (c == 0) {
            zeroCharPos = i;
        }
        if (prevC != c) {
            length_t lf = charCountsAcc[c] + charsVisited[c];
            tIn[i] = make_pair(c, lf);
        }

        charsVisited[c]++;
        prevC = c;
    }

    // Set dPair and dIndex size
    length_t arraySize = tIn.size();
    length_t size = arraySize;
    logger.logInfo("There are " + to_string(size) + " runs.");

    // Fill dPair and dIndex
    MoveLFReprBP rows;
    rows.setZeroCharPos(zeroCharPos);
    fillRows(tIn, rows, size, BWT.size());

    string moveFileName = baseFN + ".LFBP";
    rows.write(moveFileName);
    logger.logInfo("\tWrote file " + moveFileName);

    // Write non-BP version
    moveFileName = baseFN + ".LFFull";
    ofstream ofs(moveFileName, ios::binary);
    ofs.write((char*)&bwtSize, sizeof(bwtSize));
    ofs.write((char*)&size, sizeof(size));
    ofs.write((char*)&zeroCharPos, sizeof(zeroCharPos));

    for (length_t i = 0; i < size + 1; i++) {
        uint8_t runChar = rows.getRunHead(i);
        length_t inputStartPos = rows.getInputStartPos(i);
        length_t outputStartPos = rows.getOutputStartPos(i);
        length_t outputStartRun = rows.getOutputStartRun(i);
        ofs.write((char*)&runChar, sizeof(runChar));
        ofs.write((char*)&inputStartPos, sizeof(inputStartPos));
        ofs.write((char*)&outputStartPos, sizeof(outputStartPos));
        ofs.write((char*)&outputStartRun, sizeof(outputStartRun));
    }

    ofs.close();
    logger.logInfo("\tWrote file " + moveFileName);

    return size;
}

/**
 * Print usage to stdout
 */
void showUsage() {
    cout << "Usage: ./bmove-build [-l "
            "<seedLength>] [--pfp] [--preprocess] <fasta file>\n\n";
    cout << "Required arguments:\n";
    cout << "\t fasta file: Path to the FASTA file containing the reference "
            "sequences.\n\n";

    cout << "Optional arguments:\n";
    cout << "\t--preprocess:                   Preprocess the fasta file only. "
            "Do not build the "
            "index.\n";
    cout << "\t--pfp:                          Start the index construction "
            "after the prefix-free "
            "parsing step.\n";
    cout << "\t-l <seedLength> (default: " << DEFAULT_SEEDLENGTH_BMOVE
         << "): Seed length for replacing non-ACGT characters"
            ". Seed length 0 means that no seed is used.\n\n";

    cout << "Report bugs to lore.depuydt@ugent.be" << endl;
}

bool parseArguments(int argc, char* argv[], string& baseFN, int& seedLength,
                    bool& preprocessOnly, bool& pfp) {
    if (argc < 2)
        return false;

    // Set defaults
    seedLength = DEFAULT_SEEDLENGTH_BMOVE;
    preprocessOnly = false;
    pfp = false;

    // Process flags and options
    for (int i = 1; i < argc - 1; ++i) {
        string arg(argv[i]);
        if (arg == "-l" && i + 1 < argc) {
            try {
                // Convert seed length argument to integer
                seedLength = stoi(argv[++i]);

                // Check that it is positive
                if (seedLength < 0) {
                    logger.logError("Seed length must be a positive integer.");
                    return false;
                }
            } catch (const invalid_argument& e) {
                logger.logError("Seed length must be an integer.");
                return false;
            }
        } else if (arg == "--preprocess") {
            preprocessOnly = true;
        } else if (arg == "--pfp") {
            pfp = true;
        } else {
            logger.logError(
                "Unknown argument or missing value for build argument "
                "flag: " +
                arg);
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
 * Builds samplesFirst and samplesLast from the suffix array (SA) and BWT.
 * Fills samplesFirst and samplesLast at character changes in BWT for length_t.
 *
 * @param samplesFirst The samplesFirst array to fill.
 * @param samplesLast The samplesLast array to fill.
 * @param SA The suffix array.
 * @param BWT The Burrows-Wheeler transform.
 */
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast,
                  const vector<length_t>& SA, const string& BWT) {

    samplesFirst.emplace_back(SA[0]);
    for (size_t pos = 0; pos < BWT.size() - 1; pos++) {
        if (BWT[pos] != BWT[pos + 1]) {
            samplesLast.emplace_back(SA[pos]);
            samplesFirst.emplace_back(SA[pos + 1]);
        }
    }
    samplesLast.emplace_back(SA[BWT.size() - 1]);
}

/**
 * Initiates the build for the MovePhi structure.
 *
 * @param tIn The tIn map.
 * @param tOut The tOut map.
 * @param SamplesFirst The samplesFirst array.
 * @param SamplesLast The samplesLast array.
 */
void startBuildForPhiMove(map<length_t, length_t>& tIn,
                          map<length_t, length_t>& tOut,
                          vector<length_t>& samplesFirst,
                          vector<length_t>& samplesLast) {

    tIn[samplesFirst.front()] = samplesLast.back();
    tOut[samplesLast.back()] = samplesFirst.front();
    // Remove the front of the samplesFirst array
    samplesFirst.erase(samplesFirst.begin());
    // Remove the back of the samplesLast array
    samplesLast.pop_back();
    while (!samplesFirst.empty()) {
        tIn[samplesFirst.back()] = samplesLast.back();
        tOut[samplesLast.back()] = samplesFirst.back();
        samplesFirst.pop_back();
        samplesLast.pop_back();
    }
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param textLength length of the text
 */
void generatePredecessors(const MovePhiReprBP& move, SparseBitvec& pred,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength + 1, false);

    for (length_t i = 0; i < move.size(); i++) {
        predFirstBV[move.getInputStartPos(i) > 0 ? move.getInputStartPos(i) - 1
                                                 : textLength - 1] = true;
    }

    pred = SparseBitvec(predFirstBV);
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param movePhiInv MovePhiInv array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 */
void generatePredecessors(const vector<length_t>& movePhi,
                          const vector<length_t>& movePhiInv,
                          SparseBitvec& predFirst, SparseBitvec& predLast,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength + 1, false);
    vector<bool> predLastBV(textLength + 1, false);

    for (const length_t& row : movePhi) {
        predFirstBV[row > 0 ? row - 1 : textLength - 1] = true;
    }

    for (const length_t& row : movePhiInv) {
        predLastBV[row > 0 ? row - 1 : textLength - 1] = true;
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
void generatePredToRun(const vector<length_t>& samplesFirst,
                       const vector<length_t>& samplesLast,
                       vector<length_t>& firstToRun,
                       vector<length_t>& lastToRun, length_t textLength) {

    firstToRun.resize(samplesFirst.size());
    iota(firstToRun.begin(), firstToRun.end(), 0);
    sort(firstToRun.begin(), firstToRun.end(),
         [&samplesFirst, &textLength](length_t a, length_t b) {
             return (samplesFirst[a] > 0 ? samplesFirst[a] - 1
                                         : textLength - 1) <
                    (samplesFirst[b] > 0 ? samplesFirst[b] - 1
                                         : textLength - 1);
         });

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(),
         [&samplesLast, &textLength](length_t a, length_t b) {
             return (samplesLast[a] > 0 ? samplesLast[a] - 1 : textLength - 1) <
                    (samplesLast[b] > 0 ? samplesLast[b] - 1 : textLength - 1);
         });
}

/**
 * @brief Fill tUnbalanced with unbalanced intervals
 * 
 * @param tUnbalanced The map to fill with unbalanced intervals
 * @param tIn The input intervals
 * @param tOut The output intervals
 * @param maxOverlap The maximum overlap allowed
 */
void fillTUnbalanced(map<length_t, length_t>& tUnbalanced,
                     const map<length_t, length_t>& tIn,
                     const map<length_t, length_t>& tOut, length_t maxOverlap) {
    auto tInIt = tIn.begin();
    auto tOutIt = tOut.begin();

    while (tOutIt != tOut.end()) {
        length_t overlapCount = 0;
        length_t outStart = tOutIt->first;
        length_t outEnd = (std::next(tOutIt) != tOut.end())
                              ? std::next(tOutIt)->first
                              : numeric_limits<length_t>::max();

        // Skip input intervals that are before the current output interval
        while (tInIt != tIn.end() && tInIt->first < outStart) {
            ++tInIt;
        }

        // Determine if there is overlap or if the interval just started
        if (tInIt == tIn.end() || outStart < tInIt->first) {
            overlapCount = 1;
        }

        // Count overlapping input intervals
        while (tInIt != tIn.end() && tInIt->first < outEnd &&
               overlapCount < maxOverlap) {
            ++overlapCount;
            ++tInIt;
        }

        // Mark as unbalanced if overlap exceeds maxOverlap
        if (overlapCount >= maxOverlap) {
            tUnbalanced[tOutIt->second] = outStart;
        }

        ++tOutIt;
    }
}

/**
 * @brief Balance the intervals
 *
 * @param tUnbalanced The unbalanced intervals
 * @param rows The phi move array (output)
 * @param tIn The input intervals
 * @param tOut The output intervals
 * @param maxOverlap The maximum overlap allowed
 * @param textLength The length of the text
 */
void balanceIntervals(map<length_t, length_t>& tUnbalanced, MovePhiReprBP& rows,
                      map<length_t, length_t> tIn, map<length_t, length_t> tOut,
                      length_t maxOverlap, length_t textLength) {
    while (!tUnbalanced.empty()) {
        // Get the first unbalanced interval
        auto unbalancedIt = tUnbalanced.begin();
        length_t outStart = unbalancedIt->second;
        length_t inputStart = unbalancedIt->first;

        // Compute the interval difference for splitting
        auto nextInputIt = tIn.upper_bound(outStart);
        ++nextInputIt; // Todolore: this is only correct for maxOverlap = 4
        length_t intervalDiff = nextInputIt->first - outStart;

        // Insert the new intervals
        tIn.emplace(inputStart + intervalDiff, outStart + intervalDiff);
        tOut.emplace(outStart + intervalDiff, inputStart + intervalDiff);

        // Remove the balanced interval from unbalanced
        tUnbalanced.erase(inputStart);

        // Check overlapping intervals for rebalancing
        length_t outputsToCheck[] = {
            tOut.lower_bound(inputStart + intervalDiff) != tOut.begin()
                ? std::prev(tOut.lower_bound(inputStart + intervalDiff))->first
                : numeric_limits<length_t>::max(),
            outStart, outStart + intervalDiff};

        for (length_t outputStart : outputsToCheck) {
            auto endIt = tOut.upper_bound(outputStart);
            length_t endInterval = (endIt != tOut.end())
                                       ? endIt->first
                                       : numeric_limits<length_t>::max();

            // Count overlapping input intervals
            length_t overlapCount = 1;
            auto overlapInputIt = tIn.upper_bound(outputStart);
            while (overlapInputIt != tIn.end() &&
                   overlapInputIt->first < endInterval &&
                   overlapCount <= maxOverlap) {
                ++overlapCount;
                ++overlapInputIt;
            }

            // Rebalance if needed
            if (overlapCount > maxOverlap) {
                tUnbalanced[tOut[outputStart]] = outputStart;
            }
        }
    }

    rows.initialize(tIn.size(), textLength);
    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->first, it->second, 0);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second);
        i++;
    }
    length_t size = tIn.size();
    rows.setRowValues(i, textLength, textLength, size);
}

void createUnbalancedMoveTable(MovePhiReprBP& rows,
                               const map<length_t, length_t>& tIn,
                               length_t textLength) {

    rows.initialize(tIn.size(), textLength);
    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->first, it->second, 0);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second);
        i++;
    }
    length_t size = tIn.size();
    rows.setRowValues(i, textLength, textLength, size);
}

/**
 * Generate the mapping from start indices in both output intervals (phi and phi
 * inverse) to their corresponding run indices
 *
 * @param move Move structure
 * @param pred predecessor bitvector of the start samples
 * @param textLength length of the text
 */
void generatePhiRunMapping(MovePhiReprBP& move, const SparseBitvec& pred,
                           length_t textLength) {
    for (length_t i = 0; i < move.size(); i++) {
        length_t textPos = move.getOutputStartPos(i);
        length_t runIndexInText = pred.rank(textPos);
        move.setOutputStartRun(i, runIndexInText);
        assert(move.getOutputStartRun(i) == runIndexInText);
    }
}

void readSuffixArrayFile(const std::string& baseFN,
                         const std::string& extension,
                         vector<length_t>& samples, vector<length_t>& toRun,
                         length_t size, length_t nrOfRuns, bool reverse) {
    logger.logInfo("Reading suffix array samples from " + baseFN + extension +
                   "...");

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
        if (reverse) {
            // Only for the reverse suffix array, sa_val must be corrected.
            // For the reverse case, the value is too small since Big-BWT
            // puts the sentinel character at the end of the reverse text as
            // well.
            sa_val = (sa_val < size - 1) ? (sa_val + 1) : 0;
        }

        // Store sa_val in samples array at index i
        assert(sa_val >= 0 && sa_val < size);
        samples[i] = sa_val;

        // Store {sa_val, i} pair in pos_run_pairs array
        pos_run_pairs[i] = {sa_val > 0 ? sa_val - 1 : size - 1, i};
    }

    delete[] buf;

    if (!reverse) {

        logger.logInfo(
            "Creating the mapping between the predecessor bits and the "
            "runs...");

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
    logger.logInfo("Constructing run-length encoded PLCP...");

    // Construct cumulative char counts more efficiently
    vector<length_t> charCountsCumulative(256, 0);
    for (size_t i = 1; i < 256; ++i) {
        charCountsCumulative[i] = charCountsCumulative[i - 1] + charCounts[i];
    }

    logger.logDeveloper("Reading " + baseFN + ".LFBP...");
    MoveLFReprBP rows;
    if (!rows.load(baseFN)) {
        logger.logError("Error: Could not load " + baseFN + ".LFBP");
        exit(1);
    }

    logger.logDeveloper("Creating bitvectors for select support...");

    // Create bitvectors and select structures
    sdsl::bit_vector char_bv[ALPHABET];
    sdsl::bit_vector::select_1_type select[ALPHABET];

    for (length_t i = 0; i < ALPHABET; ++i) {
        char_bv[i] = sdsl::bit_vector(size, 0); // Initialize bit vector
    }

    for (length_t i = 0; i < r; ++i) {
        size_t c = rows.getRunHead(i);
        for (length_t j = rows.getInputStartPos(i);
             j < rows.getInputStartPos(i + 1); ++j) {
            char_bv[c][j] = true; // Mark positions with true
        }
    }

    for (length_t i = 0; i < ALPHABET; ++i) {
        select[i] = sdsl::bit_vector::select_1_type(
            &char_bv[i]); // Initialize select support
    }

    logger.logDeveloper("Creating PLCP...");

    std::vector<length_t> ones, zeros;
    ones.reserve(r + 1);
    zeros.reserve(r + 1);

    length_t pos = 0, gap = 0, l = 0, p = 0, p0 = 0, prev_pos = 0, prev_l = 0;
    length_t acc0 = 0, acc1 = 0;

    for (length_t i = 0; i < r; ++i) {
        if (r > 100 && i % (r / 100) == 0) {
            logger.logDeveloper("Progress: " + std::to_string(i) + " / " +
                                std::to_string(r) + " (" +
                                std::to_string((i * 100) / r) + "%)");
        }

        // Perform select on the sparse bit vector
        try {
            pos = first.select(i);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << " while selecting index " << i
                      << '\n';
            exit(1);
        }

        gap = (i == 0) ? pos : (pos - prev_pos - 1);

        l = 0;
        length_t tempIdx = first_to_run[i];
        p = rows.getInputStartPos(tempIdx);
        char c = sigma.i2c(rows.getRunHead(tempIdx));
        rows.findLF(p, tempIdx);

        if (p != 0) {
            p0 = p - 1;
            char c0 = std::upper_bound(charCountsCumulative.begin(),
                                       charCountsCumulative.end(), p0) -
                      charCountsCumulative.begin();

            // Loop until characters diverge
            while (c == c0) {
                l++;
                length_t i = p - charCountsCumulative[c - 1];

                // Perform select on the current character
                try {
                    p = select[sigma.c2i(c)](i + 1);
                } catch (const std::exception& e) {
                    std::cerr << "Error: " << e.what()
                              << " while selecting character " << c << '\n';
                    exit(1);
                }

                c = std::upper_bound(charCountsCumulative.begin(),
                                     charCountsCumulative.end(), p) -
                    charCountsCumulative.begin();
                length_t i0 = p0 - charCountsCumulative[c0 - 1];

                // Perform select on the previous character
                try {
                    p0 = select[sigma.c2i(c0)](i0 + 1);
                } catch (const std::exception& e) {
                    std::cerr << "Error: " << e.what()
                              << " while selecting character " << c0 << '\n';
                    exit(1);
                }

                c0 = std::upper_bound(charCountsCumulative.begin(),
                                      charCountsCumulative.end(), p0) -
                     charCountsCumulative.begin();
            }
        }

        // Update ones and zeros based on gap and l values
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

    // Final push to ensure last values are stored
    ones.push_back(acc1);

    logger.logDeveloper("Copying the PLCP into its correct object...");
    plcp = PLCP(size, ones, zeros);

    logger.logInfo("Constructed run-length encoded PLCP.");
}

void writePhiMoveStructures(MovePhiReprBP& phiMove, const string& baseFN,
                            const string& suffix, bool inverse,
                            length_t bwtSize) {
    // Generate and write predecessor bitvector
    logger.logInfo("\tGenerating the predecessor bitvector...");
    SparseBitvec pred;
    generatePredecessors(phiMove, pred, bwtSize);

    string predFileName =
        baseFN + ".move." + suffix + "." + (inverse ? "prdl" : "prdf");
    pred.write(predFileName);
    logger.logInfo("\tWrote file " + predFileName);

    // Generate and write phi run mapping
    logger.logInfo("\tGenerating the mapping between output interval start "
                   "positions and their run identifiers...");
    generatePhiRunMapping(phiMove, pred, bwtSize);

    string moveFileName = baseFN + ".phiBP." + suffix + (inverse ? ".inv" : "");
    phiMove.write(moveFileName);
    logger.logInfo("\tWrote file " + moveFileName);

    // Write non-BP version
    moveFileName = baseFN + ".phiFull." + suffix + (inverse ? ".inv" : "");
    ofstream ofs(moveFileName, ios::binary);

    length_t size = phiMove.size();
    ofs.write((char*)&bwtSize, sizeof(bwtSize));
    ofs.write((char*)&size, sizeof(size));

    for (length_t i = 0; i < size + 1; i++) {
        length_t inputStartPos = phiMove.getInputStartPos(i);
        length_t outputStartPos = phiMove.getOutputStartPos(i);
        length_t outputStartRun = phiMove.getOutputStartRun(i);
        ofs.write((char*)&inputStartPos, sizeof(inputStartPos));
        ofs.write((char*)&outputStartPos, sizeof(outputStartPos));
        ofs.write((char*)&outputStartRun, sizeof(outputStartRun));
    }

    ofs.close();
    logger.logInfo("\tWrote file " + moveFileName);
}

void processPhiMoveStructures(const string& baseFN,
                              map<length_t, length_t>& inputIntervals,
                              map<length_t, length_t>& outputIntervals,
                              bool inverse, length_t bwtSize,
                              length_t maxOverlap = 4) {
    { // Unbalanced move table
        MovePhiReprBP phiMove;
        createUnbalancedMoveTable(phiMove, inputIntervals, bwtSize);
        writePhiMoveStructures(phiMove, baseFN, "unbalanced", inverse, bwtSize);
    }

    { // Balanced move table
        logger.logInfo("\tCollecting the unbalanced intervals...");
        map<length_t, length_t> unbalancedIntervals;
        fillTUnbalanced(unbalancedIntervals, inputIntervals, outputIntervals,
                        maxOverlap);

        logger.logInfo("\tBalancing the intervals...");
        MovePhiReprBP phiMove;
        balanceIntervals(unbalancedIntervals, phiMove, inputIntervals,
                         outputIntervals, maxOverlap, bwtSize);
        unbalancedIntervals.clear();

        writePhiMoveStructures(phiMove, baseFN, "balanced", inverse, bwtSize);
    }
}

void processPhiAndPhiInverse(const string& baseFN,
                             vector<length_t>& samplesFirst,
                             vector<length_t>& samplesLast, length_t bwtSize) {
    // Initialize move structures for phi
    map<length_t, length_t> inputIntervals;
    map<length_t, length_t> outputIntervals;

    logger.logInfo("Initiating the samples for the phi move structures...");
    startBuildForPhiMove(inputIntervals, outputIntervals, samplesFirst,
                         samplesLast);

    samplesFirst.clear();
    samplesLast.clear();

    length_t maxOverlap = 4;

    // Log creation of move table for phi
    logger.logInfo("Creating the move table for phi...");
    processPhiMoveStructures(baseFN, inputIntervals, outputIntervals, false,
                             bwtSize, maxOverlap);

    // Log creation of move table for phi inverse
    logger.logInfo("Creating the move table for phi inverse...");
    processPhiMoveStructures(baseFN, outputIntervals, inputIntervals, true,
                             bwtSize, maxOverlap);

    inputIntervals.clear();
    outputIntervals.clear();
}

void createIndex(const string& baseFN, const string& T, bool fromFasta) {
    Alphabet<ALPHABET> sigma;
    vector<length_t> charCounts;
    writeCharCountsAndCreateAlphabet(baseFN, T, sigma, charCounts);

    {
        // read the suffix array
        vector<length_t> SA;
        readOrCreateSAWithSanityCheck(baseFN, SA, T, fromFasta);

        // build the BWT
        string BWT;
        generateBWT(T, SA, BWT);

        { // generate and write PLCP array
            logger.logInfo(
                "Generating the permuted longest common prefix array...");
            PLCP plcp(T, SA);
            plcp.write(baseFN + ".plcp");
            logger.logInfo("Wrote file " + baseFN + ".plcp");
        }

        // Create samplesLast and samplesFirst
        logger.logInfo("Sampling the suffix array values at run boundaries...");
        vector<length_t> samplesFirst;
        vector<length_t> samplesLast;
        buildSamples(samplesFirst, samplesLast, SA, BWT);

        // Clear the suffix array
        SA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".smpf");
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        logger.logInfo("Wrote file " + baseFN + ".smpl");

        vector<length_t> firstToRun;
        vector<length_t> lastToRun;

        // Generate the predToRun arrays
        logger.logInfo(
            "Mapping the predecessor bits to their corresponding runs...");
        generatePredToRun(samplesFirst, samplesLast, firstToRun, lastToRun,
                          BWT.size());

        SparseBitvec predFirst;
        SparseBitvec predLast;

        logger.logInfo(
            "Generating the predecessor bitvectors for the samples...");
        generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                             BWT.size());

        // Write predecessor structures
        predFirst.write(baseFN + ".og.prdf");
        logger.logInfo("Wrote file " + baseFN + ".og.prdf");
        predLast.write(baseFN + ".og.prdl");
        logger.logInfo("Wrote file " + baseFN + ".og.prdl");

        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        logger.logInfo("Wrote file " + baseFN + ".ftr");
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        logger.logInfo("Wrote file " + baseFN + ".ltr");

        firstToRun.clear();
        lastToRun.clear();

        processPhiAndPhiInverse(baseFN, samplesFirst, samplesLast, BWT.size());

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        createAndWriteMove(baseFN, BWT, charCounts, sigma);
        BWT.clear();
    }

    logger.logInfo("Switching to reversed text...");

    { // read or create the reverse suffix array
        vector<length_t> revSA;
        readOrCreateRevSAWithSanityCheck(baseFN, revSA, T, fromFasta);

        // build the reverse BWT
        string revBWT;
        createRevBWT(baseFN, T, revSA, sigma, revBWT);

        // Create samplesLast and samplesFirst
        logger.logInfo(
            "Sampling the reverse suffix array values at run boundaries...");
        vector<length_t> revSamplesFirst;
        vector<length_t> revSamplesLast;
        buildSamples(revSamplesFirst, revSamplesLast, revSA, revBWT);

        // Clear the reverse suffix array
        revSA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpf");
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpl");

        revSamplesFirst.clear();
        revSamplesLast.clear();

        // Create the Move structure
        logger.logInfo("Creating the reverse move table...");
        createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();
    }

    writeTagAndCompiledInfo(baseFN);
}

// Optimized helper function to replace multiple characters and log their counts
void replaceSentinel(std::string& str) {
    logger.logInfo("Replacing special characters in the BWT...");

    size_t countNull = 0; // Count for '\0'
    size_t countOne = 0;  // Count for '\1'
    size_t countTwo = 0;  // Count for '\2'

    // Traverse the string once to count and replace
    for (auto& ch : str) {
        if (ch == '\0') {
            countNull++; // Count '\0'
            ch = '$';    // Replace with '$'
        } else if (ch == '\1') {
            countOne++; // Count '\1'
            ch = '$';   // Replace with '$'
        } else if (ch == '\2') {
            countTwo++; // Count '\2'
            ch = '$';   // Replace with '$'
        }
    }

    // Log results after the single pass
    logger.logDeveloper("Found " + std::to_string(countNull) +
                        " \\0 characters.");
    logger.logDeveloper("Found " + std::to_string(countOne) +
                        " \\1 characters.");
    logger.logDeveloper("Found " + std::to_string(countTwo) +
                        " \\2 characters.");

    size_t totalCount = countNull + countOne + countTwo;

    // Error check
    if (totalCount != 1) {
        logger.logError("Found " + std::to_string(totalCount) +
                        " special characters in the BWT. Expected 1.");
        exit(1);
    }
}

void createIndexPFP(const string& baseFN) {

    // build the BWT
    string BWT;
    // Read the BWT from disk
    logger.logInfo("Reading " + baseFN + ".bwt...");
    readText(baseFN + ".bwt", BWT);

    replaceSentinel(BWT);

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
        logger.logInfo("Creating the move table...");
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
        logger.logInfo("Wrote file " + baseFN + ".smpf");
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        logger.logInfo("Wrote file " + baseFN + ".smpl");
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        logger.logInfo("Wrote file " + baseFN + ".ftr");
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        logger.logInfo("Wrote file " + baseFN + ".ltr");

        lastToRun.clear();

        { // Original phi operation + PLCP
            SparseBitvec predFirst;
            SparseBitvec predLast;

            // Generate the predecessor bitvectors
            logger.logInfo(
                "Generating the predecessor bitvectors for the samples...");
            generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                                 bwtSize);

            // Write predecessor structures
            predFirst.write(baseFN + ".og.prdf");
            logger.logInfo("Wrote file " + baseFN + ".og.prdf");
            predLast.write(baseFN + ".og.prdl");
            logger.logInfo("Wrote file " + baseFN + ".og.prdl");
            // generate and write PLCP array
            logger.logInfo(
                "Generating the permuted longest common prefix array...");
            PLCP plcp;
            constructRunLengthEncodedPLCP(predFirst, firstToRun, charCounts,
                                          bwtSize, nrOfRuns, plcp, baseFN,
                                          sigma);
            plcp.write(baseFN + ".plcp");
            logger.logInfo("Wrote file " + baseFN + ".plcp");

            firstToRun.clear();
        }

        processPhiAndPhiInverse(baseFN, samplesFirst, samplesLast, bwtSize);
    }

    logger.logInfo("Switching to reversed text...");

    {
        // build the BWT
        string revBWT;
        // Read the BWT from disk
        logger.logInfo("Reading " + baseFN + ".rev.bwt...");
        readText(baseFN + ".rev.bwt", revBWT);

        replaceSentinel(revBWT);

        // Create the Move structure
        logger.logInfo("Creating the move table...");
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
        logger.logInfo("Wrote file " + baseFN + ".rev.smpf");
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpl");

        revSamplesFirst.clear();
        revSamplesLast.clear();
    }

    writeTagAndCompiledInfo(baseFN);
}

void createBMove(const string& baseFN) {
    string T; // the text
    preprocessTextMode(baseFN, T);
    createIndex(baseFN, T, false);
}

void processFastaFile(const std::string& fastaFile, const std::string& baseFN,
                      length_t seedLength) {
    std::string T; // the (concatenated) text
    preprocessFastaFile(fastaFile, baseFN, T, seedLength, true);

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
    int seedLength;
    bool preprocessOnly;
    bool pfp;

    if (!parseArguments(argc, argv, baseFN, seedLength, preprocessOnly, pfp)) {
        showUsage();
        return EXIT_FAILURE;
    }
    logger.logInfo("Welcome to b-move's index construction!");
    logger.logInfo("Alphabet size is " + std::to_string(ALPHABET - 1) + " + 1");

    bool txtMode, fastaMode;
    std::string fastaExtension;
    determineInputMode(baseFN, txtMode, fastaMode, fastaExtension);

    if (preprocessOnly) {
        logger.logInfo("Preprocessing only...");
        if (!fastaMode) {
            logger.logError("Preprocessing only works with fasta files.");
            return EXIT_FAILURE;
        }
        std::string T; // the (concatenated) text
        preprocessFastaFile(baseFN + fastaExtension, baseFN, T, seedLength,
                            true);

        // Remove the sentinel character
        T.pop_back();

        // Write the preprocessed text to disk
        logger.logInfo("Writing concatenated uppercase sequence to disk...");
        std::ofstream ofs(baseFN);
        ofs << T;
        ofs.close();

        // Reverse the text
        logger.logInfo("Reversing text...");
        std::reverse(T.begin(), T.end());

        // Write the reversed text to disk
        logger.logInfo("Writing reversed text to disk...");
        ofs = std::ofstream(baseFN + ".rev");
        ofs << T;
        ofs.close();

        logger.logInfo("Fasta input preprocessed successfully!");
        logger.logInfo("Exiting... bye!");

        return EXIT_SUCCESS;
    }

    if (pfp) {
        logger.logInfo(
            "Starting index construction after prefix-free parsing step...");

        for (auto& ext : requiredExtensionsPFP) {
            string filename = baseFN + ext;

            // check if file with filename exists
            ifstream ifs(filename);
            if (!ifs) {
                logger.logError("Missing file: " + filename);
                logger.logError(
                    "Please run the prefix-free parsing step first.");
                return EXIT_FAILURE;
            }
        }

        createIndexPFP(baseFN);

        logger.logInfo("Index construction completed successfully!");
        logger.logInfo("Exiting... bye!");

        return EXIT_SUCCESS;
    }

    try {
        if (txtMode) {
            createBMove(baseFN);
        } else if (fastaMode) {
            processFastaFile(baseFN + fastaExtension, baseFN, seedLength);
        } else {
            throw runtime_error("No valid input files found");
        }
    } catch (const std::exception& e) {
        logger.logError("Fatal: " + string(e.what()));

        return EXIT_FAILURE;
    }

    logger.logInfo("Index construction completed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}
