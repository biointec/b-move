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
#include "move.h"
#include "moveElement.h"
#include "rindexhelpers.h"
#include "wordlength.h"

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
length_t LF_call_count = 0;
length_t elapsed_LF = 0;
#endif

using namespace std;

template <size_t S>
void
Move<S>::leftExtension(const PositionRange& inputRange, PositionRange& newRange,
                               const length_t c,
                               const bool onlyPosIndex) const {
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    // Start timing with rdtscp
    unsigned int start_aux;
    uint64_t start = __builtin_ia32_rdtscp(&start_aux);
#endif

    // Select the first occurence of cP (run index) in the input range.
    Position lfStartPos;
    if (!getNextOccurence(inputRange, lfStartPos, c)) {
        newRange = PositionRange::getEmptyPositionRange();
        return;
    }

    // Select the last occurence of cP (run index) in the input range.
    Position lfEndPos;
    getPreviousOccurence(inputRange, lfEndPos, c);


    // Perform LF and return
    lf(lfStartPos, newRange.getBeginPosMutable(), onlyPosIndex);
    lf(lfEndPos, newRange.getEndPosMutable(), onlyPosIndex);

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    // End timing with rdtscp
    unsigned int end_aux;
    uint64_t end = __builtin_ia32_rdtscp(&end_aux);

    // Check if the core has switched
    if (start_aux != end_aux) {
        cerr << "Core switch detected in leftExtension after " << LF_call_count
             << " calls" << endl;

    } else {

        elapsed_LF += end - start;
        LF_call_count++;
    }
#endif

    return;
}

template <size_t S>
bool Move<S>::getPreviousOccurence(
    const PositionRange& startPosition, Position& result,
    const length_t c) const {
    result = startPosition.getEndPos();

    while (!runHeadEquals(result.runIndex, c) &&
           result.runIndex >= startPosition.getBeginPos().runIndex) {
        result.runIndex--;
        result.posIndex = mappings[result.runIndex + 1].inputStart - 1;
    }

    return true;
}

template <size_t S>
bool Move<S>::getNextOccurence(const PositionRange& startPosition,
                                       Position& result,
                                       const length_t c) const {
    result = startPosition.getBeginPos();

    while (!runHeadEquals(result.runIndex, c) &&
           result.runIndex <= startPosition.getEndPos().runIndex) {
        result.runIndex++;
        result.posIndex = mappings[result.runIndex].inputStart;
    }

    if (result.runIndex > startPosition.getEndPos().runIndex) {
        return false;
    }

    return true;
}

template <size_t S>
void Move<S>::lf(const Position& input, Position& output,
                 const bool onlyPosIndex) const {
    // Get the output position index by using mappings and adding the offset
    // from the input position to the start of its run.
    MoveElement lfElement = this->mappings[input.runIndex];
    output.posIndex =
        lfElement.outputStart + (input.posIndex - lfElement.inputStart);

    if (onlyPosIndex)
        return;

    // Find the corresponding input run by linear search, starting from the
    // first input interval overlapping the output interval.
    length_t b = lfElement.mappingIndex;
    while (this->mappings[b].inputStart <= output.posIndex) {
        b++;
    }

    output.runIndex = b - 1;

    assert(output.posIndex >= mappings[output.runIndex].inputStart);
    assert(output.runIndex == runsSize - 1 ||
           output.posIndex < mappings[output.runIndex + 1].inputStart);

    return;
}

template <size_t S>
length_t Move<S>::getCharCountInRange(const PositionRange& range, const length_t c) const {
    
    if (c == 0) {
        // zero character, check inputRange to zeroCharPos
        if (range.getBeginPos().posIndex > zeroCharPos || range.getEndPos().posIndex < zeroCharPos) {
            // zero char does not occure in input range, return empty range
            return 0;
        }
        return 1;
    }

    PositionRange charRange;
    leftExtension(range, charRange, c, true);

    return charRange.width();

}


template <size_t S>
length_t Move<S>::getSmallerCharsCountInRange(const PositionRange& range, const length_t c) const {
    length_t count = 0;
    for (length_t smallerC = 0; smallerC < c; smallerC ++) {
        count += getCharCountInRange(range, smallerC);
    }

    return count;
}


template <size_t S>
PositionRange Move<S>::fullRange() const {
    // Create the start position (pos and run 0)
    Position start(0, 0);

    // Create the end position (last bwtIndex and last run)
    Position end(this->bwtSize - 1, runsSize - 1);

    // Return the full range
    return PositionRange(start, end);
}


template <size_t S>
void Move<S>::getRunIndex(Position& position, pair<length_t,length_t> possibleRange) const {
    // Iteratively make the posible range smaller by binary search, untill only 1 interval remains.
    while (possibleRange.second-possibleRange.first >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((possibleRange.second+possibleRange.first + 1) / 2);

        // Eliminate half of the possible range by comparing the value to the test value.
        if (this->mappings[testIndex].inputStart <= position.posIndex) {
            possibleRange.first = testIndex;
        } else {
            possibleRange.second = testIndex-1;
        }
    }

    position.runIndex = possibleRange.first;

    assert(position.posIndex >= mappings[possibleRange.first].inputStart);
    assert(possibleRange.first == runsSize-1 || position.posIndex < mappings[possibleRange.first+1].inputStart);

    return;
}


template <size_t S>
PositionRange Move<S>::computeRunIndices(const PositionRange& positionRange) const {
    pair<length_t, length_t> possibleRange(positionRange.getBeginPos().runIndex, positionRange.getEndPos().runIndex);
    
    Position beginPos = positionRange.getBeginPos();
    Position endPos = positionRange.getEndPos();
    
    getRunIndex(beginPos, possibleRange);
    getRunIndex(endPos, possibleRange);
    
    return PositionRange(beginPos, endPos);
}


template <size_t S>
size_t Move<S>::printMemSize(ofstream& out) {

    size_t totalSize = sizeof(bwtSize);

    cout << endl << "Sizes in Move structure:" << endl;

    size_t mappingsSize = runsSize * sizeof(MoveElement);
    cout << "Mappings: " << mappingsSize << "\n";
    totalSize += mappingsSize;

    totalSize += sizeof(zeroCharRun) + sizeof(zeroCharPos);

    cout << "Total Move size:" << totalSize << endl << endl;

    return totalSize;
}


template <size_t S>
void Move<S>::printData() {
    // Print the contents of mappings (loop through mappings)
    cout << "Mappings:" << endl;
    for (length_t i = 0; i < runsSize; i ++) {
        cout << i << ": (" << this->mappings[i].inputStart << ", " << this->mappings[i].outputStart << ") Index output interval is " << this->mappings[i].mappingIndex << endl;
    }

    cout << endl;
}

template <size_t S> bool Move<S>::load(const string& baseFile, bool verbose) {
    string fileName = baseFile + ".move";

    ifstream ifs(fileName, ios::binary);
    if (!ifs) {
        return false;
    }

    // Load the bwtSize, amount of input intervals and the alphabet size.
    ifs.read((char*)&bwtSize, sizeof(bwtSize));
    ifs.read((char*)&runsSize, sizeof(runsSize));
    ifs.read((char*)&zeroCharPos, sizeof(zeroCharPos));

    // Load rows: allocate memory and read.
    mappings.reserve(runsSize + 1);
    for (length_t i = 0; i < runsSize + 1; i++) {
        mappings[i].load(ifs);
    }

    ifs.close();
    return true;
}

template class Move<ALPHABET>;
