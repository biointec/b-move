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
#ifndef RINDEX_H
#define RINDEX_H

#include "move.h"

#include "plcp.h"
#include "rindexhelpers.h"
#include "search.h"
#include "tkmer.h"

#include <fstream>  // used for reading in files
#include <iostream> // used for printing
#include <google/sparse_hash_map>
#include <cstdint>

class RIndex;
typedef bool (RIndex::*BRExtraCharPtr)(length_t, const BRSample&, BRSample&) const;

class RIndex {
    private:
        const std::string baseFile; // the basefile of the reference text
        length_t textLength;        // the length of the text

        std::vector<length_t> counts;   // the counts array of the reference genome
        Alphabet<ALPHABET> sigma;       // alphabet

        PLCP plcp;        // PLCP (permuted longest common prefix) array

        // First and last sa samples of each input interval, needed for keeping the toehold
        sdsl::int_vector<> samplesFirst;
        sdsl::int_vector<> samplesLast;
        sdsl::int_vector<> revSamplesFirst;
        sdsl::int_vector<> revSamplesLast;

        // Predecessor structures.
        SparseBitvec predFirst;
        SparseBitvec predLast;

        // Contain predecessor to run mappings.
        sdsl::int_vector<> firstToRun;
        sdsl::int_vector<> lastToRun;

        Move<ALPHABET> move;    // The Move data structure of the text
        Move<ALPHABET> moveR;   // The Move data structure of the reverse text

        // direction variables
        thread_local static Direction dir;          // the direction of the index
        
        thread_local static BRExtraCharPtr extraChar; // pointer to extra char method (for direction)

        // stacks for search schemes
        thread_local static std::vector<std::vector<BRPosExt>> stacks;  // stacks of nodes for the different partitions
        thread_local static std::vector<BitParallelED> matrices;        // alignment matrices for the different partitions
        thread_local static BitParallelED inTextMatrix;                 // alignment matrix for the whole read

        // sparse hash info
        const size_t wordSize = 10; // the size of the mers to be stored in a table
        google::sparse_hash_map<Kmer, BRSample, KmerHash> table; // hashtable that contains all wordSize-mers


        // ----------------------------------------------------------------------------
        // PREPROCESSING ROUTINES
        // ----------------------------------------------------------------------------

        /**
         * Private helper function that reads in all the necessary files
         * @param baseFile the baseFile of the files that will be read in
         * @param verbose if true the steps will we written to cout
         */
        void fromFiles(const std::string& baseFile, bool verbose);

        /**
         * Read a binary file and stores content in vector
         * @param filename File name
         * @param array output array (contents will be overwritten)
         * @returns True if successful, false otherwise
         */
        static bool readArray(const std::string& filename, std::vector<length_t>& array) {
            std::ifstream ifs(filename, std::ios::binary);
            if (!ifs) {
                return false;
            }
            ifs.seekg(0, std::ios::end);
            array.resize(ifs.tellg() / sizeof(length_t));
            ifs.seekg(0, std::ios::beg);
            ifs.read((char*)&array[0], array.size() * sizeof(length_t));

            return true;
        }

        /**
         * Populate the hash table
         * @param verbose if steps are written to cout
         */
        void populateTable(bool verbose);

        // ----------------------------------------------------------------------------
        // ROUTINES FOR ACCESSING DATA STRUCTURE
        // ----------------------------------------------------------------------------
        /**
         * Finds the text position given the position that succeeds it in the SA,
         * in other words, finds SA[i-1] given SA[i]
         * @param pos Position in the text
         * @returns Position in the text that precedes given position in the SA
         */
        length_t phi(length_t pos) const;

        /**
         * Finds the text position given the position that precedes it in the SA,
         * in other words, finds SA[i+1] given SA[i]
         * @param pos Position in the text
         * @returns Position in the text that succeeds given position in the SA
         */
        length_t phiInverse(length_t pos) const;

        /**
         * Get toehold for the given character to match.
         * @param range Range in the sa array in which to match char and get toehold.
         * @param c Character to match.
         * @return Index in the SA where the toehold starts.
        */
        length_t getToehold(const PositionRange& range, const length_t c) const;

        /**
         * Get toehold for the given character to match in the reverse text.
         * @param range Range in the sa array in which to match char and get toehold.
         * @param c Character to match.
         * @return Index in the SA where the toehold starts.
        */
        length_t getToeholdRev(const PositionRange& range, const length_t c) const;

        /**
         * Get an initial toehold for the full range in the BWT.
         * @returns A toehold for the full range in the BWT.
        */
        length_t getInitialToehold() const;

        /**
         * Function that returns the number of occurrences before an index of
         * the symbol at symbolindex in the alphabet
         * @param symbolIndex the index of the the symbol in the alphabet to
         * count the occrences of at index index
         * @param index the index whose entry for symbol in the occurrences table
         * is asked
         * @return the number of occurrences of the symbol before index in the
         * bwt
         */
        /*length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
            return fwdbwt.rank(symbolIndex, index);
        }*/

        /**
         * Same as getNumberOfOcc, but now in the bwt of the reversed text
         * @param symbolIndex the index in the alphabet to get the number of
         * occurrences of
         * @param index the index in the occurrences table, for which the number
         * of occurrences of the symbol at symbolindex is asked.
         * @return the number of occurrences at index index in the occurrences
         * table of the bwt of the reversed text for the symbol at symbolIndex
         * in the alphabet
         */
        /*length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
            return revbwt.rank(symbolIndex, index);
        }*/

        /**
         * Function that returns the number of occurrences before the index of
         * all symbols smaller than the symbol at symbolindex in the bwt
         * @param symbolIndex the index in the alphabet whose number of baseFile
         * occurrences is queried.
         * @param index the index whose entry for symbol in the prefixoccurrences
         * table is asked
         * @return the number of occurrences of symbols smaller than symbol at
         * symbolindex before index index in the bwt
         */
        /*length_t getNumberOfCumOcc(length_t symbolIndex, length_t index) const {
            return fwdbwt.cumOcc(symbolIndex, index);
        }*/

        /**
         * Function that returns the number of occurrences before the index of
         * all symbols smaller than the symbol at symbolindex in the bwt of the
         * reversed text
         * @param symbolIndex the index in the alphabet whose number of baseFile
         * occurrences is queried.
         * @param index the index whose entry for symbol in the
         * prefixoccurrences table of the reverse text is asked
         * @return the number of occurrences of symbols smaller than symbol at
         * symbolindex before index index in the bwt
         */
        /*length_t getNumberOfCumOccRev(length_t symbolIndex, length_t index) const {
            return revbwt.cumOcc(symbolIndex, index);
        }*/

        // ----------------------------------------------------------------------------
        //  APPROXIMATE PATTERN MATCHING
        // ----------------------------------------------------------------------------

        /**
         * Matches a search recursively with a depth first approach (each branch
         * of the tree is fully examined until the backtracking condition is
         * met) using the edit distance metric. This function uses all
         * optimizations for eliminating redundancy in the edit distance metric
         * @param search the search to follow
         * @param startMatch the approximate match found for all previous
         * partitions of the search
         * @param occ a datastructure with matches of the complete search, if  such
         * a match is found it will be added to this datastructure
         * @param parts the parts of the pattern
         * @param counters the performance counters
         * @param idx the index of the partition to match, defaults to 1 as an
         * exact search for the zeroth partition is assumed
         * @param descPrevDir the descendants of the previous direction,
         * defaults to empty vector
         * @param initPrevDir the initialization eds of the previous direction,
         * defaults to empty vector
         * @param descNotPrevDir the descendants of the other direction,
         * defaults to empty vector
         * @param initNotPrevDir the initialization eds of the other direction,
         * defaults to empty vector
         */
        void recApproxMatchEditOptimized(
            const Search& search, const BROcc& startMatch, BROccurrences& occ,
            const std::vector<Substring>& parts, Counters& counters,
            const int& idx = 1,
            const std::vector<BRPosExt>& descPrevDir = std::vector<BRPosExt>(),
            const std::vector<uint>& initPrevDir = std::vector<uint>(),
            const std::vector<BRPosExt>& descNotPrevDir = std::vector<BRPosExt>(),
            const std::vector<uint>& initNotPrevDir = std::vector<uint>());

        /**
         * Finds the ranges of cP using the principle explained in the paper of
         * Lahm
         * @param positionInAlphabet the position in alphabet of the character
         * that is added in the front
         * @param sampleOfP the sample of pattern P
         * @param childSample the sample of pattern cP
         */
        bool findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                             const BRSample& sampleOfP,
                                             BRSample& childSample) const;

        /**
         * Finds the ranges of Pc using the principle explained in the paper of
         * Lahm
         * @param positionInAlphabet the position in alphabet of the character c
         * that is added in the back
         * @param sampleOfP the sample of pattern P
         * @param childSample the sample of pattern Pc
         */
        bool findRangesWithExtraCharForward(length_t positionInAlphabet,
                                            const BRSample& sampleOfP,
                                            BRSample& childSample) const;

        /**
         * Helper function for the approximate matching. This function fills in
         * the matrix for the current node at the current row and goes deeper
         * for the next part is that is necessary The function true if the search
         * should backtrack.
         * @param bpED alignment matrix per part
         * @param clus the cluster corresponding to the final column of the
         * matrix
         * @param currentnode the node for which the matrix is filled in
         * @param s the current search
         * @param idx the idx of the current part
         * @param parts the parts of the pattern
         * @param occ Datastructure with the in-index and in-text occurrences, if an
         * occurrence is found it will be added to this datastructure
         * @param counters the performance counters
         * @param initOther eds of descendants in other direction
         * @param descOther descendants in other direction
         * @param remainingDesc the remaining descendants on the current
         * branch, that are already created but aren't checked yet and might
         * need to be checked for the next part, defaults to an empty vector
         * @return false if the search can continue along this branch for the
         * current part, true if the search can backtrack
         */
        bool branchAndBound(BitParallelED& bpED, BRCluster& clus,
                            const BRPosExt& currentNode, const Search& s,
                            const length_t& idx,
                            const std::vector<Substring>& parts, BROccurrences& occ,
                            Counters& counters, const std::vector<uint>& initOther,
                            const std::vector<BRPosExt>& descOther,
                            const std::vector<BRPosExt>& remainingDesc = {});

        /**
         * Goes deeper in a search if a valid approximate match is found in the
         * cluster
         * @param cluster the cluster to search for a valid approximate match
         * @param nIdx the idx of next part to research
         * @param s the search
         * @param parts the parts of the pattern
         * @param occ Datastructure with the in-index and in-text occurrences, if an
         * occurrence is found it will be added to this datastructure
         * @param cnts the performance counters
         * @param descsOtherD the descendants of the other direction,
         * defaults to empty vector
         * @param initOhterD the initialization eds of the other direction,
         * defaults to empty vector
         * @param remDesc the remaining descendants on the current
         * branch, that are already created but aren't checked yet and need to
         * be checked for the next part, defaults to an empty vector
         */
        void goDeeper(BRCluster& cluster, const length_t& nIdx, const Search& s,
                      const std::vector<Substring>& parts, BROccurrences& occ,
                      Counters& cnts,
                      const std::vector<BRPosExt>& descndansndansOtherD,
                      const std::vector<uint>& intitOtherD,
                      const std::vector<BRPosExt>& remDesc);


        /**
         * Converts a match in the suffix array to matches in the text.
         * @param matchInSA the match that will be converted
         * @returns a vector with the corresponding text occurrences
         */
        std::vector<TextOcc> convertToMatchesInText(const BROcc& matchInSA) const;

        /**
         * Find the begin positions of the range that needs to be verified in-text
         * @param rangeSA the range over the SA of the partial in-index match that
         * needs to be verified in text
         * @param startDiff the highest possible difference between the start of the
         * partial match and that of a valid full match in the text
         * @param shift the shift of the current node (default = 0)
         * @returns the lowest possible begin position for each of the values in the
         * range
         */
        std::vector<length_t> getBeginPositions(const BRSample& sample,
                                                length_t startDiff,
                                                length_t shift = 0) const {

            std::vector<length_t> positions = textPositionsFromSample(sample);

            for (length_t pos: positions) {
                pos += shift - startDiff;
            }

            return positions;
        }

        // ----------------------------------------------------------------------------
        // EXTEND ROUTINES
        // ----------------------------------------------------------------------------

        /**
         * Pushes all the children corresponding to the node with ranges equal
         * to ranges
         * @param sample the sample to get the children of
         * @param stack the stack to push the children on
         * @param counters the performance counters
         * @param row the row of the parentNode (defaults to 0)
         */
        void extendBRPos(const BRSample& sample, std::vector<BRPosExt>& stack,
                        Counters& counters, length_t row = 0) const;

        /**
         * Pushes all the children corresponding to the this position
         * @param pos the position to get the children of
         * @param stack the stack to push the children on
         * @param counters the performance counters
         */
        void extendBRPos(const BRPosExt& pos, std::vector<BRPosExt>& stack,
                        Counters& counters) const;

    public:

        size_t printMemSize() {
            // textLength, baseFile
            size_t totalSize = sizeof(textLength);

            std::ofstream devnull("/dev/null");

            size_t countsSize = counts.size() * sizeof(length_t);
            std::cout << "counts: " << countsSize << "\n";
            totalSize += countsSize;
            size_t sigmaSize = NUM_CHAR * sizeof(int) + sigma.size()*2 - 1;
            std::cout << "sigma: " << sigmaSize << "\n";
            totalSize += sigmaSize;

            totalSize += plcp.printMemSize();

            totalSize += move.printMemSize(devnull);
            totalSize += moveR.printMemSize(devnull);

            size_t samplesFirstSize = samplesFirst.serialize(devnull);
            cout << "samplesFirst: " << samplesFirstSize << "\n";
            totalSize += samplesFirstSize;
            size_t samplesLastSize = samplesLast.serialize(devnull);
            cout << "samplesLast: " << samplesLastSize << "\n";
            totalSize += samplesLastSize;
            size_t revSamplesFirstSize = revSamplesFirst.serialize(devnull);
            cout << "revSamplesFirst: " << revSamplesFirstSize << "\n";
            totalSize += revSamplesFirstSize;
            size_t revSamplesLastSize = revSamplesLast.serialize(devnull);
            cout << "revSamplesLast: " << revSamplesLastSize << "\n";
            totalSize += revSamplesLastSize;

            size_t firstToRunSize = firstToRun.serialize(devnull);
            cout << "firstToRun: " << firstToRunSize << "\n";
            totalSize += firstToRunSize;
            size_t lastToRunSize = lastToRun.serialize(devnull);
            cout << "lastToRun: " << lastToRunSize << "\n";
            totalSize += lastToRunSize;
            size_t predFirstSize = predFirst.serialize(devnull);
            cout << "predFirst: " << predFirstSize << "\n";
            totalSize += predFirstSize;
            size_t predLastSize = predLast.serialize(devnull);
            cout << "predLast: " << predLastSize << "\n";
            totalSize += predLastSize;

            std::cout << "total r index size: " << totalSize << "\n";
            return totalSize;
        }

        // ----------------------------------------------------------------------------
        // INITIALIZATION ROUTINES
        // ----------------------------------------------------------------------------

        /**
         * Constructor
         * @param baseFile base filename of the files that contain the info
         * @param verbose write info to cout
         */
        RIndex(const std::string& baseFile, bool verbose = true, length_t wordSize = 10) 
            : baseFile(baseFile), wordSize(wordSize) {
            
            // read in files
            fromFiles(baseFile, verbose);

            // populate table
            populateTable(verbose);
        }

        /**
         * Get the complete range of this index
         * @returns an SARangePair with both ranges the complete range of the
         * index
         */
        SAPositionRangePair getCompleteRange() const {
            return SAPositionRangePair(move.fullRange(), moveR.fullRange());
        }

        /**
         * Get the initial sample of the complete range over the index
         * @return a BRSample with complete ranges over the index
        */
        BRSample getInitialSample() const {
            SAPositionRangePair completeRanges = getCompleteRange();
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
            return BRSample(completeRanges, 0, 0, 0);
#else
            return BRSample(completeRanges, getInitialToehold(), 0, 0);
#endif
        }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATASTRUCTURE
    // ----------------------------------------------------------------------------

        /**
         * Get the length of the text
        */
        const length_t getTextSize() const {
            return textLength;
        }

        /**
         * @returns the wordsize of the mers stored in the table
         */
        length_t getWordSize() const {
            return wordSize;
        }

        /**
         * @returns the ranges of a single character in this index
         */
        SAPositionRangePair getPositionRangeOfSingleChar(char c) const {
            assert(sigma.inAlphabet(c));
            unsigned int i = sigma.c2i(c);

            PositionRange fullRange = move.fullRange();
            PositionRange fullRevRange = moveR.fullRange();

            if (i < sigma.size() - 1) {
                PositionRange range(Position::getPositionWithIndex(counts[i], fullRange.getBeginPos()), Position::getPositionWithIndex(counts[i+1]-1, fullRange.getEndPos()));
                range.setRunIndicesValid(false);
                PositionRange revRange(Position::getPositionWithIndex(counts[i], fullRevRange.getBeginPos()), Position::getPositionWithIndex(counts[i+1]-1, fullRevRange.getEndPos()));
                revRange.setRunIndicesValid(false);

                return SAPositionRangePair(range, revRange);
            }

            PositionRange range(Position::getPositionWithIndex(counts[i], fullRange.getBeginPos()), Position::getPositionWithIndex(textLength-1, fullRange.getEndPos()));
            range.setRunIndicesValid(false);
            PositionRange revRange(Position::getPositionWithIndex(counts[i], fullRevRange.getBeginPos()), Position::getPositionWithIndex(textLength-1, fullRevRange.getEndPos()));
            revRange.setRunIndicesValid(false);

            return SAPositionRangePair(range, revRange);
        }

        /**
         * @param c character
         * @return Sample of a single character in the index
        */
        BRSample getSampleOfSingleChar(char c) const {
            SAPositionRangePair rangeOfSingleChar = getPositionRangeOfSingleChar(c);
            return BRSample(rangeOfSingleChar, 0 /*Correct?*/, 0, 1);
        }

        /**
         * Looks up the sample corresponding to p in the hashtable. Assumes
         * p is of size wordSize
         * @param p, the substring to find the ranges of
         * @returns the ranges corresponding to substring p, if no pair can be
         * found returns empty ranges
         */
        BRSample lookUpInKmerTable(const Substring& p) const {
            Kmer k(p.tostring());

            auto it = table.find(k);
            if (it != table.end()) {
                return it->second;
            }

            return BRSample();
        }

        // ----------------------------------------------------------------------------
        // ROUTINES FOR EXACT MATCHING
        // ----------------------------------------------------------------------------

        /**
         * Calculates the exact matches to the string in the index and returns them
         * with their CIGAR string
         * @param s the string to match in the reference genome
         * @param counters the performance counters
         * @returns a sorted vector containing the start positions of all exact
         * substring matches of s in the reference sequence
         */
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
        size_t exactMatchesOutput(const std::string& s,
                                                Counters& counters) {
#else
        std::vector<TextOcc> exactMatchesOutput(const std::string& s,
                                                Counters& counters) {
#endif
            if (s.size() == 0) {
                return {};
            }
            BRSample sample = getInitialSample();
            assert(dir == FORWARD);

            std::vector<TextOcc> tOcc;

            // bool broken = false;
            for (const char& c : s) {
                assert(sigma.c2i(c) < ALPHABET && sigma.c2i(c) >= 0);
                findRangesWithExtraCharForward(sigma.c2i(c), sample, sample);
                counters.incNodeCounter();
                if (!sample.isValid()) {
                    break;
                }
            }
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
            return sample.width();
        }
#else

            // THE CIGAR for the matches
            std::vector<std::pair<char, uint>> CIGAR = {
                std::make_pair('M', s.size())};

            const auto& p = textPositionsFromSample(sample);
            // if (broken && sample.isValid()) {
            //     // for each of the positions check if the substring in the text
            //     // equals the string

            //     for (const auto& pos : p) {
            //         if (pos + s.size() <= textLength &&
            //             text.compare(pos, s.size(), s) == 0) {
            //             tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR);
            //         }
            //     }
            //     counters.inTextStarted += p.size();
            //     counters.usefulCigarsInText += tOcc.size();
            //     counters.abortedInTextVerificationCounter += p.size() - tOcc.size();
            //     counters.cigarsInTextVerification += tOcc.size();
            // } else {
                // No verification needed, just add the cigar
                for (const auto& pos : p) {
                    tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR);
                }
                counters.cigarsInIndex += p.size();
            // }

            counters.totalReportedPositions += tOcc.size();
            // sort the vector and return
            std::sort(tOcc.begin(), tOcc.end());
            return tOcc;
        }
#endif

        /**
         * This function matches a string exactly and returns the ranges in the
         * sa and saRev
         * @param pattern the string to match
         * @param counters the performance counters
         * @returns the pair of ranges of this pattern
         */
        BRSample matchStringBidirectionally(const Substring& pattern,
                                               Counters& counters) const {
            return matchStringBidirectionally(pattern, getInitialSample(),
                                              counters);
        }

        /**
         * This function matches a string exactly starting form startRange
         * @param pattern the string to match
         * @param startSample the sample to search in
         * @param counters the performance counters
         * @returns the pair of ranges
         */
        BRSample matchStringBidirectionally(const Substring& pattern,
                                            BRSample startSample,
                                            Counters& counters) const;

        /**
         * Adds one character and updates the range. If the character can't be
         * added the range will be set to an empty range
         * @param c the character to be added (in the current direction of the
         * index)
         * @param sample the sample to extend
         * @param counters the performance counters
         */
        bool addChar(const char& c, BRSample& sample, Counters& counters) const;

        // ----------------------------------------------------------------------------
        // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
        // ----------------------------------------------------------------------------

        /**
         * Resizes to the required number of stacks and reserves space on each
         * stack, so that each stack can match the entire pattern
         * @param number, the number of stacks required
         * @param size, the size of the pattern
         */
        void reserveStacks(const length_t number, const length_t size) {
            stacks.resize(number);
            length_t stackSize = size * sigma.size();
            for (auto& stack : stacks) {
                stack.reserve(stackSize);
            }
        }
        /**
         * Reset the in-text matrices to be empty matrices
         * @param number the number of partitions
         */
        void resetMatrices(const length_t number) {
            matrices.resize(2 * number);

            for (auto& matrix : matrices) {
                matrix.reset();
            }
        }

        /**
         * Set the sequence for the in-text verification matrix
         * @param s the sequence to be set
         */
        void setInTextMatrixSequence(const Substring& s) {
            inTextMatrix.setSequence(s);
        }

        // ----------------------------------------------------------------------------
        // ROUTINES FOR APPROXIMATE MATCHING
        // ----------------------------------------------------------------------------

        /**
         * Matches the pattern approximately. All matches are at most a certain
         * edit distance away from the pattern
         * @param pattern the pattern to match
         * @param maxED the maximum edit distance
         * @param counters the performance counters
         * @returns a vector with matches which contain a range (the range of
         * the text that matched) and the edit distance this substring is away
         * from the pattern
         */
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
        const size_t approxMatchesNaive(const std::string& pattern,
                                                    length_t maxED,
                                                    Counters& counters);

#else
        const std::vector<TextOcc> approxMatchesNaive(const std::string& pattern,
                                                    length_t maxED,
                                                    Counters& counters);

#endif
        /**
         * Sets the search direction of the r-index
         * @param d the direction to search in, either FORWARD or BACKWARD
         */
        void setDirection(Direction d) {
            dir = d;
            extraChar = (d == FORWARD) ? &RIndex::findRangesWithExtraCharForward
                                       : &RIndex::findRangesWithExtraCharBackward;
        }

        /**
         * Matches a search recursively with a depth first approach (each branch
         * of the tree is fully examined until the backtracking condition is
         * met) using hamming distance metric
         * @param search the search to follow
         * @param startMatch the approximate match found for all previous parts
         * of the search
         * @param occ a datastructure with matches of the complete search, if  such
         * a match is found it will be added to this datastructure
         * @param counters the performance counters
         * @param idx the index of the partition to match, defaults to 1 as an
         * exact search for the zeroth part is assumed
         */
        void recApproxMatchHamming(const Search& s, const BROcc& startMatch,
                                BROccurrences& occ,
                                const std::vector<Substring>& parts,
                                Counters& counters, const int& idx = 1);

        /**
         * Entry to the recusive approximate matching procedure for the edit
         * distance
         * @param search the search to follow
         * @param startMatch the match containing the SA ranges corresponding to
         * the match of the first part of the search
         * @param occ a datastructure with matches of the complete search, if  such
         * a match is found it will be added to this datastructure
         * @param parts the parts of the pattern
         * @param counters the performance counters
         * @param idx the index of the partition to match, defaults to 1 as an
         * exact search for the zeroth partition is assumed
         */
        void recApproxMatchEditOptimizedEntry(const Search& search,
                                            const BROcc& startMatch,
                                            BROccurrences& occ,
                                            const std::vector<Substring>& parts,
                                            Counters& counters,
                                            const int& idx = 1) {

            if (startMatch.getRanges().width() > 0) {
                counters.approximateSearchStarted++;
                recApproxMatchEditOptimized(search, startMatch, occ, parts,
                                            counters, idx);
                return;
            }

            const auto& partialStarts = getBeginPositions(startMatch.getSample(), 
                parts[search.getLowestPartProcessedBefore(idx)].begin() + search.getMaxED());

            for (auto pos: partialStarts) {
                occ.addTextOcc(Range(pos, pos + startMatch.getDepth()), startMatch.getDistance());
            }
        }
            

        /**
         * Matches a search recursively with a depth first approach (each branch
         * of the tree is fully examined until the backtracking condition is
         * met) using the edit distance metric. This function does not use any
         * optimizations for eliminating redundancy in the edit distance metric.
         * It simply matches the current part starting from startrange and each
         * node found that has an edit distance between the lower and upper
         * bound is used to start a search for the next part
         * @param s the search to follow
         * @param startMatch the approximate match found for all previous
         * partitions of the search
         * @param occ a datastructure with matches of the complete search, if
         * such a match is found it will be added to this datastructure
         * @param parts the parts of the pattern
         * @param counters the performance counters
         * @param idx the index of the partition to match, defaults to 1 as an
         * exact search for the zeroth partition is assumed
         */
        void recApproxMatchEditNaive(const Search& s, const BROcc& startMatch,
                                    BROccurrences& occ,
                                    const std::vector<Substring>& parts,
                                    Counters& counters, const int& idx);
        // ----------------------------------------------------------------------------
        // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
        // ----------------------------------------------------------------------------

        /**
         * Get all positions in the text from a sample
         * @param sample The BRSample
         * @return A vector with text positions
        */
        std::vector<length_t> textPositionsFromSample(const BRSample& sample) const;

        /**
         * Find unqiue text occurrences for the hamming distance
         */
        std::vector<TextOcc> getTextOccHamming(BROccurrences& occ,
                                            Counters& counters) const {
            counters.totalReportedPositions += occ.textOccSize();
            // erase in-index doubles
            occ.eraseDoublesIndex();

            // Get the in-index occurrences
            const auto& broccs = occ.getIndexOccurrences();

            // The size of the match
            length_t size = (broccs.empty()) ? 0 : broccs[1].getDepth();

            for (const auto& brOcc : broccs) {
                // Get the range
                Range saRange(
                    brOcc.getRanges().getRangeSA().getBeginPos().posIndex,
                    brOcc.getRanges().getRangeSA().getEndPos().posIndex + 1);
                counters.totalReportedPositions += saRange.width();

                std::vector<TextOcc> textoccs = convertToMatchesInText(brOcc);

                for (TextOcc textocc: textoccs) {
                    std::vector<std::pair<char, uint>> CIGAR = {std::make_pair('M', size)};
                    textocc.setCigar(CIGAR);
                    occ.addTextOcc(textocc);
                }
            }

            // remove doubles
            occ.eraseDoublesText();
            // Generate string output
            occ.generateOutput();

            return occ.getTextOccurrences();
        }

        /**
         * Find the unique text occurrences
         * @param occ occurrences, both in-text and in-index
         * @param maxED the maximal allowed edit distance
         * @param counters performance counters
         */
        std::vector<TextOcc> getUniqueTextOccurrences(BROccurrences& occ,
                                                    const length_t& maxED,
                                                    Counters& counters);
};

#endif