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
#include "rindex.h"
#include <cassert>

using namespace std;

thread_local Direction RIndex::dir = BACKWARD;
thread_local BRExtraCharPtr RIndex::extraChar;

thread_local vector<vector<BRPosExt>> RIndex::stacks;
thread_local vector<BitParallelED> RIndex::matrices;


#ifndef ENABLE_BENCHMARK_FUNCTIONALITY

/**
 * Read a binary file and stores content in sdsl int_vector
 * @param filename File name
 * @param intVector output int_vector (contents will be overwritten)
 * @returns True if successful, false otherwise
 */
bool readIntVector(const std::string& filename, sdsl::int_vector<>& intVector) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        return false;
    }
    intVector.load(ifs);
    ifs.close();
    return true;
}

#endif

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------


void RIndex::fromFiles(const string& baseFile, bool verbose) {

    // read the compiled info
    {
        ifstream ifs(baseFile + ".comp");
        if (!ifs) {
            std::cout <<
                "The index was built with an older version of bmove-build or "
                "the compile info (" +
                baseFile +
                ".comp) is missing! "
                "Proceed with caution or rebuild the index!." << std::endl;
        } else {

            size_t build_length_t_size;
            ifs >> build_length_t_size;
            if (build_length_t_size != sizeof(length_t)) {
                size_t thisSize = sizeof(length_t) * 8;
                build_length_t_size *= 8;

                throw runtime_error(
                    "The index was built with a compiled version that uses " +
                    to_string(build_length_t_size) +
                    "-bit numbers, while the current programme was compiled "
                    "using " +
                    to_string(thisSize) +
                    "-bit numbers. Recompile the programme with the correct "
                    "THIRTY_TWO flag set.");
            }
        }
    }
    // Check the build tag of the file
    {
        length_t tag = 0;
        ifstream ifs(baseFile + ".tag");
        if (!ifs) {
            std::cout << 
                "The index was built with an older version of bmove-build or "
                "the tag is missing! "
                "Proceed with caution or update re-build the index." << std::endl;
        } else {
            ifs >> tag;
            ifs.close();
            if (tag < BMOVE_BUILD_INDEX_TAG) {
                std::cout << 
                    "The index was built with an older version of "
                    "bmove-build. "
                    "Proceed with caution or update re-build the index." << std::endl;
            } else if (tag > BMOVE_BUILD_INDEX_TAG) {
                std::cout << 
                    "The index was built with a newer version of "
                    "bmove-build. "
                    "Proceed with caution or update re-build the index." << std::endl;
            }
        }
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".cct" << "..." << std::endl;
    }

    // Read the counts table
    vector<length_t> charCounts(256, 0);
    if (!readArray(baseFile + ".cct", charCounts)) {
        throw runtime_error("Cannot open file: " + baseFile + ".cct");
    }

    length_t cumCount = 0; // Cumulative character counts
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        counts.push_back(cumCount);
        cumCount += charCounts[i];
    }
    textLength = cumCount;
    sigma = Alphabet<ALPHABET>(charCounts);

    if (verbose) {

#ifndef ENABLE_BENCHMARK_FUNCTIONALITY

        std::cout << "Reading " << baseFile << ".plcp..." << std::endl;
    }
    if (!plcp.read(baseFile + ".plcp")) {
        throw runtime_error("Cannot open file: " + baseFile + ".plcp");
    }

    if (verbose) {
#endif

        std::cout << "Reading " << baseFile << ".move..." << std::endl;
    }
    if (!move.load(baseFile, verbose)) {
        throw runtime_error("Error loading move file: " + baseFile + ".move");
    }

    if (verbose) {
#ifndef ENABLE_BENCHMARK_FUNCTIONALITY

        std::cout << "Reading " << baseFile << ".smpf..." << std::endl;
    }
    if (!readIntVector(baseFile + ".smpf", samplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpf");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".smpl..." << std::endl;
    }
    if (!readIntVector(baseFile + ".smpl", samplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpl");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".prdf..." << std::endl;
    }
    if (!predFirst.read(baseFile + ".prdf")) {
        throw runtime_error("Cannot open file: " + baseFile + ".prdf");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".prdl..." << std::endl;
    }
    if (!predLast.read(baseFile + ".prdl")) {
        throw runtime_error("Cannot open file: " + baseFile + ".prdl");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".ftr..." << std::endl;
    }
    if (!readIntVector(baseFile + ".ftr", firstToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ftr");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".ltr..." << std::endl;
    }
    if (!readIntVector(baseFile + ".ltr", lastToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ltr");
    }

    if (verbose) {
#endif


        std::cout << "Reading " << baseFile << ".rev.move..." << std::endl;
    }
    if (!moveR.load(baseFile + ".rev", verbose)) {
        throw runtime_error("Error loading reverse move file: " + baseFile + ".rev.move");
    }

    if (verbose) {
#ifndef ENABLE_BENCHMARK_FUNCTIONALITY

        std::cout << "Reading " << baseFile << ".rev.smpf..." << std::endl;
    }
    if (!readIntVector(baseFile + ".rev.smpf", revSamplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpf");
    }

    if (verbose) {
        std::cout << "Reading " << baseFile << ".rev.smpl..." << std::endl;
    }
    if (!readIntVector(baseFile + ".rev.smpl", revSamplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpl");
    }

    if (verbose) {
#endif
    }
}

void RIndex::populateTable(bool verbose) {
    if (verbose) {
        cout << "Populating FM-range table with " << wordSize << "-mers...";
    }
    cout.flush();

    Kmer::setWordSize(wordSize);

    table.resize(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize
    setDirection(FORWARD);

    string word;
    vector<BRPosExt> stack;
    Counters counters;
    extendBRPos(getInitialSample(), stack, counters);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if (curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);
            auto sample = curr.getSample();
            auto newPosRange = move.computeRunIndices(sample.getRanges().getRangeSA());
            SAPositionRangePair ranges(newPosRange, sample.getRanges().getRangeSARev());
            sample.setRanges(ranges);
            table.insert(make_pair(k, std::move(sample)));

        } else // add extra characters
            extendBRPos(curr.getSample(), stack, counters, curr.getRow());
    }
    if (verbose) {
        cout << "done." << endl;
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------

BRSample RIndex::matchStringBidirectionally(const Substring& pattern,
                                            BRSample sampleOfPrev,
                                            Counters& counters) const {

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        if (!addChar(c, sampleOfPrev, counters)) {
            // sampleOfPrev was made empty
            break;
        }
    }

    return sampleOfPrev;
}

bool RIndex::addChar(const char& c, BRSample& startSample,
                     Counters& counters) const {

    int posInAlphabet = sigma.c2i((unsigned char)c);
    if (posInAlphabet > -1) {

        assert(posInAlphabet < ALPHABET);
        if ((this->*extraChar)(posInAlphabet, startSample, startSample)) {
            // each character that we look at is a new node that is visited
            counters.incNodeCounter();
            return true;
        }
    }
    // the range is now empty
    return false;
}


length_t RIndex::phi(length_t pos) const {
    
    // find predecessor of pos
    length_t predRank = predFirst.predecessorRankCircular(pos);
    length_t pred = predFirst.select(predRank);

    // distance from predecessor
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // check if phi(SA[0]) is not called
    assert(firstToRun[predRank] > 0);

    length_t prev_sample = samplesLast[firstToRun[predRank]-1];

    return (prev_sample + delta) % textLength;
}


length_t RIndex::phiInverse(length_t pos) const {

    // find predecessor of pos
    length_t predRank = predLast.predecessorRankCircular(pos);
    length_t pred = predLast.select(predRank);

    // distance from predecessor
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // check if phiInverse(SA[n-1]) is not called
    assert(lastToRun[predRank] < samplesFirst.size()-1);

    length_t prev_sample = samplesFirst[lastToRun[predRank]+1];

    return (prev_sample + delta) % textLength;
}

length_t RIndex::getToehold(const PositionRange& range, const length_t c) const {

    length_t endRunHead = move.getRunHead(range.getEndPos());

    Position previousOcc;
    move.getPreviousOccurence(range, previousOcc, c);
    
    if (endRunHead == c) {
        return samplesFirst[previousOcc.runIndex];
    }
    return samplesLast[previousOcc.runIndex];

}


length_t RIndex::getToeholdRev(const PositionRange& range, const length_t c) const {

    length_t endRunHead = moveR.getRunHead(range.getEndPos());

    Position previousOcc;
    moveR.getPreviousOccurence(range, previousOcc, c);

    if (endRunHead == c) {
        return revSamplesFirst[previousOcc.runIndex];
    }
    return revSamplesLast[previousOcc.runIndex];

}


length_t RIndex::getInitialToehold() const {
    return samplesLast[samplesLast.size() - 1];
}


// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

bool RIndex::findRangesWithExtraCharBackward(
    length_t positionInAlphabet, const BRSample& sampleOfP,
    BRSample& childSample) const {

    // first make the trivial range by searching cP using B
    PositionRange fwdRangeOfP = sampleOfP.getRanges().getRangeSA();

    // If the run indices of revRangeOfP are not valid, compute them
    if (! fwdRangeOfP.getRunIndicesValid()) {
        fwdRangeOfP = move.computeRunIndices(fwdRangeOfP);
    }

    assert(positionInAlphabet < ALPHABET);
    PositionRange newFwdRange;
    move.leftExtension(fwdRangeOfP, newFwdRange, positionInAlphabet, false);

    if (newFwdRange.empty()) {
        SAPositionRangePair childRanges = SAPositionRangePair(newFwdRange, newFwdRange);
        // set empty range
        childSample = BRSample(childRanges, 0, 0, 0);
        return false;
    }

    // then make the less trivial range by counting the sizes of the ranges
    // of dP for all d < c using B

    // first get the start of the range we are looking for of the parent
    PositionRange revRangeOfP = sampleOfP.getRanges().getRangeSARev();

    if(newFwdRange.width() == fwdRangeOfP.width()) {
        // only increment length
        SAPositionRangePair childRanges = SAPositionRangePair(newFwdRange, revRangeOfP);
        childSample = BRSample(childRanges, sampleOfP.getTextPos(), 
                                sampleOfP.getOffset() + 1, sampleOfP.getLength() + 1);
        return true;
    }
    length_t s = revRangeOfP.getBeginPos().posIndex;

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = move.getSmallerCharsCountInRange(fwdRangeOfP, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    PositionRange newRevRange = PositionRange(
        Position::getPositionWithIndex(s + x, revRangeOfP.getBeginPos()), 
        Position::getPositionWithIndex(s + x + newFwdRange.width() - 1, revRangeOfP.getEndPos()));
    newRevRange.setRunIndicesValid(false);


    SAPositionRangePair childRanges = SAPositionRangePair(newFwdRange, newRevRange);

    // aP occurs for some a != c, sample position must be updated

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    length_t newTextpos = 0;
#else
    length_t newTextpos = getToehold(fwdRangeOfP, positionInAlphabet);
#endif

    // update sample with newTextPos, reset offset, increment length
    childSample = BRSample(childRanges, newTextpos, 0, sampleOfP.getLength() + 1);
    return true;
}

bool RIndex::findRangesWithExtraCharForward(
    length_t positionInAlphabet, const BRSample& sampleOfP,
    BRSample& childSample) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    PositionRange revRangeOfP = sampleOfP.getRanges().getRangeSARev();

    // If the run indices of revRangeOfP are not valid, compute them
    if (! revRangeOfP.getRunIndicesValid()) {
        revRangeOfP = moveR.computeRunIndices(revRangeOfP);
    }

    assert(positionInAlphabet < ALPHABET);
    PositionRange newRevRange;
    moveR.leftExtension(revRangeOfP, newRevRange, positionInAlphabet, false);

    if (newRevRange.empty()) {
        SAPositionRangePair childRanges = SAPositionRangePair(newRevRange, newRevRange);
        // set empty range
        childSample = BRSample(childRanges, 0, 0, 0);
        return false;
    }

    // first get the start of the range we are looking for of the parent
    PositionRange fwdRangeOfP = sampleOfP.getRanges().getRangeSA();

    if(newRevRange.width() == revRangeOfP.width()) {
        // only increment length
        SAPositionRangePair childRanges = SAPositionRangePair(fwdRangeOfP, newRevRange);
        childSample = BRSample(childRanges, sampleOfP.getTextPos(), 
                                sampleOfP.getOffset(), sampleOfP.getLength() + 1);
        return true;
    }

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)
    length_t s = fwdRangeOfP.getBeginPos().posIndex;

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = moveR.getSmallerCharsCountInRange(revRangeOfP, positionInAlphabet);

    // make the new range
    PositionRange newFwdRange = PositionRange(
        Position::getPositionWithIndex(s + x, fwdRangeOfP.getBeginPos()), 
        Position::getPositionWithIndex(s + x + newRevRange.width() - 1, fwdRangeOfP.getEndPos()));
    newFwdRange.setRunIndicesValid(false);

    SAPositionRangePair childRanges = SAPositionRangePair(newFwdRange, newRevRange);
    
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    length_t newTextpos = 0;
#else
    // Pa occurs for some a != c, sample position must be updated
    length_t newTextpos = textLength - 1 - getToeholdRev(revRangeOfP, positionInAlphabet);
#endif

    // update sample with newTextPos, reset offset, increment length
    childSample = BRSample(childRanges, newTextpos, 
                            sampleOfP.getLength(), sampleOfP.getLength() + 1);
    return true;
}

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
const size_t RIndex::approxMatchesNaive(const string& pattern,
                                                 length_t maxED,
                                                 Counters& counters) {
#else
const vector<TextOcc> RIndex::approxMatchesNaive(const string& pattern,
                                                 length_t maxED,
                                                 Counters& counters) {
#endif

    counters.resetCounters();
    BROccurrences occurrences;

    BitParallelED matrix;
    matrix.setSequence(pattern);
    matrix.initializeMatrix(maxED);

    setDirection(FORWARD);

    vector<BRPosExt> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendBRPos(getInitialSample(), stack, counters, 0);

    length_t lastcol = pattern.size();

    while (!stack.empty()) {
        const BRPosExt currentNode = stack.back();
        stack.pop_back();
        length_t row = currentNode.getDepth();

        if (row >= matrix.getNumberOfRows()) {
            continue;
        }

        bool valid = matrix.computeRow(row, currentNode.getCharacter());

        if (!valid) {
            // backtrack
            continue;
        }

        if (matrix.inFinalColumn(row)) {

            // full pattern was matched
            if (matrix(row, lastcol) <= maxED) {
                occurrences.addIndexOcc(currentNode, matrix(row, lastcol));
            }
        }

        extendBRPos(currentNode, stack, counters);
    }
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    // Count the total width of all occurrences:
    size_t totalWidth = 0;
    for (const auto& occ : occurrences.getIndexOccurrences()) {
        totalWidth += occ.getWidth();
    }
    return totalWidth;
#else

    return getUniqueTextOccurrences(occurrences, maxED, counters);
#endif
}

void RIndex::recApproxMatchEditNaive(const Search& s, const BROcc& startMatch,
                                      BROccurrences& occ,
                                      const vector<Substring>& parts,

                                      Counters& counters, const int& idx) {
    const Substring& p = parts[s.getPart(idx)];   // this part
    const length_t& maxED = s.getUpperBound(idx); // maxED for this part
    const length_t& minED = s.getLowerBound(idx); // minED for this part
    const length_t& pSize = p.size();
    const Direction& dir = s.getDirection(idx); // direction

    // set the direction
    setDirection(dir);

    length_t matrixID = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    auto& bpED = matrices[matrixID];

    if (!bpED.sequenceSet()) {
        bpED.setSequence(p);
    }

    bpED.initializeMatrix(maxED, {(uint)startMatch.getDistance()});

    if (bpED.getLastColumn(0) == pSize && bpED(0, pSize) <= maxED) {
        // an occurrence found by gapping entire part
        if (s.isEnd(idx)) {
            occ.addIndexOcc(startMatch.getSample(), bpED(0, pSize),
                         startMatch.getDepth());
        } else {
            // go to the next index
            recApproxMatchEditNaive(s,
                                    BROcc(startMatch.getSample(),
                                          bpED(0, pSize),
                                          startMatch.getDepth()),
                                    occ, parts, counters, idx + 1);
            // set direction correct again
            setDirection(dir);
        }
    }

    auto& stack = stacks[idx]; // stack for this partition

    extendBRPos(startMatch.getSample(), stack, counters);

    while (!stack.empty()) {
        const auto currentNode = stack.back();
        stack.pop_back();

        const length_t& row = currentNode.getRow();
        const length_t firstCol = bpED.getFirstColumn(row);
        const length_t lastCol = bpED.getLastColumn(row);

        bool valid = bpED.computeRow(row, currentNode.getCharacter());

        if (!valid) {
            // backtracking
            continue;
        }

        if (lastCol == pSize && bpED(row, pSize) <= maxED &&
            bpED(row, pSize) >= minED) {
            if (s.isEnd(idx)) {
                occ.addIndexOcc(currentNode.getSample(), bpED(row, pSize),
                             startMatch.getDepth() + currentNode.getDepth());
            } else {
                // go deeper in search
                Direction originalDir = dir;
                recApproxMatchEditNaive(
                    s,
                    BROcc(currentNode.getSample(), bpED(row, pSize),
                          startMatch.getDepth() + currentNode.getDepth()),
                    occ, parts, counters, idx + 1);
                // set direction correct again
                setDirection(originalDir);
            }
        }
        if (firstCol == pSize) {
            // final cell of matrix filled in => backtracking
            continue;
        }

        // extend to the children
        extendBRPos(currentNode, stack, counters);
    }
}

void RIndex::recApproxMatchEditOptimized(
    const Search& s, const BROcc& startMatch, BROccurrences& occ,
    const vector<Substring>& parts, Counters& counters, const int& idx,
    const vector<BRPosExt>& descPrevDir, const vector<uint>& initPrevDir,
    const vector<BRPosExt>& descNotPrevDir,
    const vector<uint>& initNotPrevDir) {

    // shortcut Variables
    const Substring& p = parts[s.getPart(idx)];      // this part
    const length_t& maxED = s.getUpperBound(idx);    // maxED for this part
    const Direction& dir = s.getDirection(idx);      // direction
    const bool& dSwitch = s.getDirectionSwitch(idx); // has direction switched?
    auto& stack = stacks[idx];                       // stack for this partition
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    BitParallelED& bpED = matrices[matrixIdx]; // matrix for this partition

    // get the correct initED and descendants based on switch
    const vector<uint>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<BRPosExt>& descendants =
        dSwitch ? descNotPrevDir : descPrevDir;
    const vector<uint>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<BRPosExt>& descOther = dSwitch ? descPrevDir : descNotPrevDir;

    // set the direction
    setDirection(dir);

    // calculate necessary increase for first column of bandmatrix
    vector<uint> initED;
    if (initEds.empty()) {
        initED = vector<uint>(1, startMatch.getDistance());
    } else {
        length_t prevED =
            (dSwitch ? *min_element(initEds.begin(), initEds.end())
                     : initEds[0]);
        length_t increase = startMatch.getDistance() - prevED;
        initED = vector<uint>(initEds.size());
        for (size_t i = 0; i < initED.size(); i++) {
            initED[i] = initEds[i] + increase;
        }
    }

    // encode the sequence of this partition in the matrix if this has not been
    // done before
    if (!bpED.sequenceSet())
        bpED.setSequence(p);

    // initialize bit-parallel matrix
    bpED.initializeMatrix(maxED, initED);

    // initialize  cluster
    BRCluster clus(bpED.getSizeOfFinalColumn(), maxED, startMatch.getDepth(),
                 startMatch.getShift());

    if (bpED.inFinalColumn(0)) {
        // the first row is part of the final column of the banded matrix
        // Update the first cell of the cluster with the startmatch and the
        // value found at the psizeth column of the initialization row the
        // character does not matter as the first cell of a cluster is never
        // a descendant of any of the other cells in the cluster, so this
        // cell will not be reused for the next part of the pattern
        clus.setValue(0, BRPosExt((char)0, startMatch.getSample(), 0),
                      bpED(0, p.size()));
    }

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        length_t maxRow = bpED.getNumberOfRows() - 1;

        for (length_t i = 0;
             i < descendants.size() && descendants[i].getDepth() <= maxRow;
             i++) {

            if (branchAndBound(
                    bpED, clus, descendants[i], s, idx, parts, occ, counters,
                    initOther, descOther,
                    {descendants.begin() + i + 1, descendants.end()})) {
                return;
            }
        }
        if (descendants.back().getDepth() == maxRow) {
            //  no more rows to possibly check
            return;
        }

        BRSample sample = descendants.back().getSample();
        if (dSwitch) {
            // after a switch the ranges of final descendant should be
            // updated
            sample = startMatch.getSample();
        }

        // push children of final descendant
        extendBRPos(sample, stack, counters, descendants.back().getDepth());

    } else { // get the initial nodes to check
        extendBRPos(startMatch.getSample(), stack, counters);
    }

    while (!stack.empty()) {
        const BRPosExt currentNode = stack.back();
        stack.pop_back();

        if (branchAndBound(bpED, clus, currentNode, s, idx, parts, occ,
                           counters, initOther, descOther)) {

            continue;
        }

        // continue the search for children of this node in-index
        extendBRPos(currentNode, stack, counters);
    }
}

bool RIndex::branchAndBound(BitParallelED& bpED, BRCluster& clus,
                             const BRPosExt& currentNode, const Search& s,
                             const length_t& idx,
                             const vector<Substring>& parts, BROccurrences& occ,
                             Counters& counters, const vector<uint>& initOther,
                             const vector<BRPosExt>& descOther,
                             const vector<BRPosExt>& remainingDesc) {

    // compute, in a bit-parallel manner, a single row of the ED matrix
    const length_t row = currentNode.getDepth();
    bool validED = bpED.computeRow(row, currentNode.getCharacter());

    // check if we have reached the final column of the matrix
    if (bpED.inFinalColumn(row)) {
        // update the cluster
        length_t clusIdx = clus.size() + row - bpED.getNumberOfRows();
        clus.setValue(clusIdx, currentNode,
                      bpED(row, bpED.getNumberOfCols() - 1));

        if (!validED || bpED.onlyVerticalGapsLeft(row)) {
            // no need to further explore this branch for this part -> go to
            // next part

            goDeeper(clus, idx + 1, s, parts, occ, counters, descOther,
                     initOther, remainingDesc);
            return true;
        }
    }

    return !validED;
}

void RIndex::goDeeper(BRCluster& cluster, const length_t& nIdx, const Search& s,
                       const vector<Substring>& parts, BROccurrences& occ,
                       Counters& counters, const vector<BRPosExt>& descOtherD,
                       const vector<uint>& initOtherD,
                       const vector<BRPosExt>& remDesc) {

    bool isEdge = s.isEdge(nIdx - 1);
    const auto& lowerBound = s.getLowerBound(nIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nIdx == parts.size()) {
            auto matches = cluster.reportCentersAtEnd();
            for (const auto& match : matches) {
                if (match.isValid() && match.getDistance() >= lowerBound) {
                    occ.addIndexOcc(match);
                }
            }
        } else {
            BROcc match = cluster.reportDeepestMinimum(this->dir);
            if (match.isValid() && match.getDistance() >= lowerBound) {
                // go deeper in search
                Direction originalDir = this->dir;
                recApproxMatchEditOptimized(s, match, occ, parts, counters,
                                            nIdx, {}, {}, descOtherD,
                                            initOtherD);
                // set direction back again
                setDirection(originalDir);
            }
        }

        return;
    }

    // one of the later stages will return to this point, so keep track of
    // the descendants and eds at this branch
    vector<BRPosExt> descendants;
    vector<uint> initEds;

    BROcc newMatch = cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (!newMatch.isValid()) {
        // no centre above lower bound found
        return;
    }

    // add the remaining descendants (copy)
    descendants.insert(descendants.end(), remDesc.begin(), remDesc.end());

    // reset the depth of all the descendants
    for (length_t i = 0; i < descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    length_t maxEDNext = s.getUpperBound(nIdx);

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.getDirectionSwitch(nIdx);

    if (switchAfter) {
        // switching direction as this is not the end of a search direction,
        // this means we'll get back here, thus range of newmatch should be
        // deepest point in branch
        if (!descendants.empty()) {
            newMatch.setSample(descendants.back().getSample());

            //  edit distance for search in other direction should be lowest
            //  value possible
            newMatch.setDistance(*min_element(initEds.begin(), initEds.end()));
        }

        Direction originalDir = this->dir;

        recApproxMatchEditOptimized(s, newMatch, occ, parts, counters, nIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatchEditOptimized(s, newMatch, occ, parts, counters, nIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);
    }
}


void RIndex::recApproxMatchHamming(const Search& s, const BROcc& startMatch,
                                    BROccurrences& occ,
                                    const vector<Substring>& parts,
                                    Counters& counters, const int& idx) {

    // shortcut variables
    const Substring& p = parts[s.getPart(idx)]; // the current part
    const length_t& pSize = p.size();           // the size of the current part
    const Direction& d = s.getDirection(idx);   // direction of current part
    const length_t& maxED = s.getUpperBound(idx); // upper bound of current part
    const length_t& minED = s.getLowerBound(idx); // lower bound of currentpart
    setDirection(d);

    // create vector for the scores
    vector<length_t> vec(p.size() + 1, 0);
    // set root element of vector to the distance of the startmatch
    vec[0] = startMatch.getDistance();
    // get stack for current part
    auto& stack = stacks[idx];

    extendBRPos(startMatch.getSample(), stack, counters);

    while (!stack.empty()) {
        const BRPosExt node = stack.back();
        stack.pop_back();

        // update the vector
        length_t row = node.getRow();
        vec[row] = vec[row - 1] + (node.getCharacter() != p[row - 1]);

        if (vec[row] > maxED) {
            // backtrack
            continue;
        }

        if (row == pSize) {
            // end of part
            if (vec[row] >= minED) {
                // valid occurrence
                BROcc match = BROcc(node.getSample(), vec[row],
                                    startMatch.getDepth() + pSize);
                if (s.isEnd(idx)) {
                    // end of search
                    occ.addIndexOcc(match);
                } else {
                    // continue search
                    recApproxMatchHamming(s, match, occ, parts, counters,
                                          idx + 1);
                    setDirection(s.getDirection(idx));
                }
            }
            continue;
        }
        extendBRPos(node, stack, counters);
    }
}


void RIndex::extendBRPos(const BRSample& parentSample,
                          vector<BRPosExt>& stack, Counters& counters,
                          length_t row) const {

    // iterate over the entire alphabet
    for (length_t i = 1; i < sigma.size(); i++) {

        BRSample sampleForNewChar;

        // check if this character occurs in the specified range
        assert(i < ALPHABET);
        if ((this->*extraChar)(i, parentSample, sampleForNewChar)) {
            // push this range and character for the next iteration
            stack.emplace_back(sigma.i2c(i), sampleForNewChar, row + 1);

            counters.incNodeCounter();
        }
    }
}

void RIndex::extendBRPos(const BRPosExt& pos, vector<BRPosExt>& stack,
                          Counters& counters) const {
    extendBRPos(pos.getSample(), stack, counters, pos.getDepth());
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<length_t> RIndex::textPositionsFromSample(const BRSample& sample) const {
    assert(sample.getTextPos() >= sample.getOffset());

    vector<length_t> textPositions;
    textPositions.reserve(sample.getRanges().width());

    if (!sample.isValid()) {
        return textPositions;
    }

    length_t startPos = sample.getTextPos() - sample.getOffset();
    length_t currentPos = startPos;

    textPositions.push_back(currentPos);

    while (plcp[currentPos] >= sample.getLength()) {
        currentPos = phi(currentPos);
        textPositions.push_back(currentPos);
    }

    currentPos = startPos;
    while (true) {
        if (currentPos == getInitialToehold()) break;
        currentPos = phiInverse(currentPos);
        if (plcp[currentPos] < sample.getLength()) break;
        textPositions.push_back(currentPos);
    }
    return textPositions;
}

vector<TextOcc> RIndex::convertToMatchesInText(const BROcc& saMatch) const {

    vector<TextOcc> textMatches;
    textMatches.reserve(saMatch.getRanges().width());

    vector<length_t> textPositions = textPositionsFromSample(saMatch.getSample());

    for (length_t pos: textPositions) {

        length_t startPos = pos + saMatch.getShift();

        length_t endPos = startPos + saMatch.getDepth();

        textMatches.emplace_back(Range(startPos, endPos),
                                 saMatch.getDistance());
    }
    return textMatches;
}

vector<TextOcc> RIndex::getUniqueTextOccurrences(BROccurrences& occ,
                                                       const length_t& maxED,
                                                       Counters& counters) {

    // increment reporte position counter
    counters.totalReportedPositions += occ.textOccSize();

    // erase equal occurrences from the in-index occurrences
    occ.eraseDoublesIndex();

    // convert the in-index occurrences to in-text occurrences
    const vector<BROcc>& broccs = occ.getIndexOccurrences();
    for (const BROcc& brocc : broccs) {

        Range saRange(brocc.getSample().getRanges().getRangeSA().getBeginPos().posIndex,
                      brocc.getSample().getRanges().getRangeSA().getEndPos().posIndex + 1);

        // increment reported positions counter
        counters.totalReportedPositions += saRange.width();

        vector<TextOcc> textoccs = convertToMatchesInText(brocc);

        for (TextOcc textocc: textoccs) {
            occ.addTextOcc(textocc);
        }
    }
    // erase equal occurrences from the in-text occurrences
    occ.eraseDoublesText();

    // find the non-redundant occurrences
    vector<TextOcc> nonRedundantOcc;
    nonRedundantOcc.reserve(occ.textOccSize());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    length_t prevDepth = numeric_limits<length_t>::max();
    length_t prevED = maxED + 1;

    const auto& textocc = occ.getTextOccurrences();
    for (const auto& o : textocc) {
        // find the difference between this and the previous
        // occurrence
        auto diff = abs_diff<length_t>(o.getRange().getBegin(), prevBegin);

        if (diff == 0) {
            // same location -> skip
            continue;
        }

        if (diff <= maxDiff) {
            // check if this later occurrence is better than the
            // previous one
            if (o.getDistance() > prevED ||
                (o.getDistance() == prevED &&
                 o.getRange().width() >= prevDepth)) {
                continue;
            }

            // prev was worse so pop_back
            nonRedundantOcc.pop_back();
        }

        prevBegin = o.getRange().getBegin();
        prevED = o.getDistance();
        prevDepth = o.getRange().width();

        nonRedundantOcc.emplace_back(o);
    }

    return nonRedundantOcc;
}
