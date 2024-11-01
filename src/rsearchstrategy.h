/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Lore Depuydt <lore.depuydt@ugent.be> and        *
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

#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#include "definitions.h"    // for length_t, MAX_K, MappingMode, Distance...
#include "rindexhelpers.h"   // for TextOcc, Occurrences
#include "rindex.h" // for the IndexInterface
#include "reads.h"          // for ReadBundle
#include "search.h"         // for Search, SearchScheme
#include "substring.h"      // for Substring

#include <algorithm> // for max, move, min
#include <cassert>   // for assert
#include <cstdint>   // for uint16_t
#include <ios>       // for ifstream, basic_ios
#include <iterator>  // for move_iterator, back_insert_iterator
#include <memory>    // for allocator, allocator_traits<>::value_type
#include <numeric>   // for accumulate
#include <set>       // for set
#include <stdexcept> // for runtime_error, invalid_argument
#include <string>    // for string, operator+, char_traits, to_string
#include <utility>   // for pair, move
#include <vector>    // for vector, _Bit_iterator

#include <sys/stat.h> // for stat on POSIX systems
#include <unistd.h>   // For access on POSIX systems
#define PATH_SEPARATOR "/"

#define Distribution std::vector<int>

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is an abstract class. Every derived class should be able to create
// searches for a given value of k. This abstract base class handles the
// partitioning (either with values provided in the derived class or default
// uniform values) and approximate matching (either hamming or edit distance)
class SearchStrategy;

// Pointer to a partition function
typedef void (SearchStrategy::*PartitionPtr)(const std::string&,
                                             std::vector<Substring>&,
                                             const int&, const int&,
                                             std::vector<SARangePair>&,
                                             Counters&) const;

// Pointer to function that starts the index on a particular search
typedef void (SearchStrategy::*StartIdxPtr)(const Search&, const FMOcc&,
                                            Occurrences&,
                                            std::vector<Substring>&, Counters&,
                                            const int&) const;

// Pointer to match function (ALL or BEST scenario)
typedef void (SearchStrategy::*MatchPtrSE)(ReadBundle&, const length_t,
                                           Counters&, std::vector<TextOcc>&);

// pointer to filter function (hamming or edit distance)
typedef std::vector<TextOcc> (SearchStrategy::*FilterPtr)(
    Occurrences&, length_t, Counters&, const ReadBundle& bundle) const;

// Pointer to function that generates SAM lines for single-end
typedef void (SearchStrategy::*GenerateOutputSEPtr)(ReadBundle&, size_t, int,
                                                    std::vector<TextOcc>&,
                                                    Counters& counters) const;
typedef TextOcc (*CreateUnmappedSEPtr)(const ReadBundle&);

/**
 * Abstract class for searching a pattern in an index. Supports both all- and
 * best-alignment, as well as single-end mapping. Is based on search
 * schemes for different values of k
 */
class SearchStrategy {

    // ATTRIBUTES
    //-------------------------------------------------------------------------------
  protected:
    IndexInterface&
        index;        // reference to the index of the text that is searched
    std::string name; // the name of the particular search strategy
    DistanceMetric distanceMetric; // which distance metric to use
  private:
    // helper type for best-matching in paired end alignment. The bool
    // indicates whether the distance (= index in OccVector) has been
    // processed and the second element is the approximate matches for this
    // distance
    typedef std::pair<bool, std::vector<TextOcc>> BoolAndVector;
    typedef std::vector<BoolAndVector> OccVector;

    // variables for getting info about strategy used
    PartitionStrategy partitionStrategy; // the partitioning strategy
    MappingMode mode;                    // whether to find all or best matches

    // pointers for correct partitioning and correct distance metric
    PartitionPtr partitionPtr; // pointer to the partition method
    StartIdxPtr startIdxPtr;   // pointer to start method (hamming or
                               // (naive/optimized) edit distance)

    // pointers to functions for mapping mode (ALL or BEST)
    MatchPtrSE matchPtrSE; // pointer to the correct match function SE
                           // alignments

    // pointers for correct filtering and in-text verification
    FilterPtr
        filterPtr; // pointer to the correct filter (edit or hamming based)
    // pointer for generating Single End output lines (SAM or RHS)
    GenerateOutputSEPtr generateOutputSEPtr;
    CreateUnmappedSEPtr createUnmappedSEPtr;

    bool samOutput = false;
#ifndef NO_REPORT
    bool generateSAMLines = true; // whether to generate SAM lines, currently
                                  // only false when mapping single-end to infer
                                  // the fragment size and orientation
#else
    bool generateSAMLines = false; // whether to generate SAM lines, currently
                                  // only false when mapping single-end to infer
                                  // the fragment size and orientation

#endif
    bool unmappedSAM = true; // whether to generate SAM lines for unmapped reads
    bool xaTag = false; // whether to use the XA tag for secondary alignments

    // ----------------------------------------------------------------------------

    // PROTECTED FUNCTIONS
    // ----------------------------------------------------------------------------

  protected:
    /**
     * Constructor
     * @param index the index to be used
     * @param p the partitioning strategy to be used
     * @param distanceMetric the distance metric to be used
     * @param mode whether all or best matches should be found
     * @param sm whether single end mode should be used
     */
    SearchStrategy(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric distanceMetric, MappingMode mode,
                   SequencingMode sm);

    /**
     * Copy constructor
     */
    SearchStrategy(const SearchStrategy& other) = default;

    // PROTECTED FUNCTIONS: PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual uint32_t calculateNumParts(unsigned int maxED) const = 0;

    // PROTECTED FUNCTIONS: PARTITIONING (STATIC)
    // ----------------------------------------------------------------------------

    /**
     * Default function that retrieves the begin positions for  static
     * partitioning. If derived class does not implement this function then
     * uniform positions are given.
     * @param numParts  how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this
     * function)
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double> getBegins(const int& numParts,
                                                const int& maxScore) const {
        std::vector<double> b;
        double u = 1.0 / numParts;
        for (int i = 1; i < numParts; i++) {
            b.push_back(i * u);
        }
        return b;
    }
    // ----------------------------------------------------------------------------

    // PROTECTED FUNCTIONS: PARTITIONING (DYNAMIC)
    // ----------------------------------------------------------------------------

    /**
     * Default function that retrieves the seeding positions for dynamic
     * partitioning. If derived class does not implement this function then
     * uniform seeds are given.
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this default
     * function)
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double>
    getSeedingPositions(const int& numParts, const int& maxScore) const {

        double u = 1.0 / (numParts - 1);
        std::vector<double> s;
        for (int i = 1; i < numParts - 1; i++) {
            s.push_back(i * u);
        }
        return s;
    }

    /**
     * Helper function for dynamic partitioning. Seeds the parts.
     * @param pattern the pattern to partition
     * @param parts empty vector tro which the seeds are added
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed edit or hamming distance
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @returns the number of characters used by the seeding operation
     */
    int seed(const std::string& pattern, std::vector<Substring>& parts,
             const int& numParts, const int& maxScore,
             std::vector<SARangePair>& exactMatchRanges) const;
    /**
     * Function that retrieves the weights for dynamic partitioning.
     * If derived class does not implement this function then uniform weights
     * are given with a preference for the two edges.
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this default
     * function)
     * @returns vector with weights
     */
    virtual const std::vector<int> getWeights(const int& numParts,
                                              const int& maxScore) const {
        std::vector<int> w(numParts, 1);
        w.front() = 2;
        w.back() = 2;
        return w;
    }
    // ----------------------------------------------------------------------------

    // PROTECTED FUNCTIONS: ACCESS
    // ----------------------------------------------------------------------------

    /**
     * Find the maximum supported distance score of this strategy.
     */
    virtual length_t getMaxSupportedDistance() const = 0;

    /**
     * Static function that calculates the number of elements in an OccVector.
     */
    static length_t numElements(const OccVector& ov) {
        length_t num = 0;
        for (const auto& p : ov) {
            num += p.second.size();
        }
        return num;
    }

    /**
     * Static helper function for matchApproxBestPlusX. Combines the forward and
     * reverse complemented matches, removes duplicates and only keeps matches
     * that have a distance score of at most max.
     * @param ovFW the forward occurrences
     * @param ovRC the reverse complement occurrences
     * @param best the best found distance score
     * @param max the maximal allowed distance score
     * @returns a vector with all unique occurrences, sorted on distance score
     */
    static std::vector<TextOcc> combineOccVectors(OccVector& ovFW,
                                                  OccVector& ovRC,
                                                  length_t best, length_t max);
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS
    // ----------------------------------------------------------------------------

  private:
    // PRIVATE FUNCTIONS: PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, either by uniform range or
     * uniform size
     * @param pattern the pattern to be split
     * @param parts the vector containing the substrings of this pattern,
     * will be cleared and filled during the execution of this method. If
     * the splitting fails for some reason, the vector will be empty
     * @param numParts  how many parts are needed
     * @param maxScore the maximum allowed edit distance
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partition(const std::string& pattern, std::vector<Substring>& parts,
                   const int& numParts, const int& maxScore,
                   std::vector<SARangePair>& exactMatchRanges,
                   Counters& counters) const;

    /**
     * Helper function for dynamic partitioning. This function extends the parts
     * so that nothing of the pattern is not allocated to any part. This does
     * not keep track of the ranges over the suffix array, so should only be
     * called if this does not matter (e.g. when the parts that can be extended
     * all correspond to empty ranges)
     * @param pattern the pattern that is split
     * @param parts  the current parts, they are updated so that all characters
     * of pattern are part of exactly one part
     */
    void extendParts(const std::string& pattern,
                     std::vector<Substring>& parts) const;

    /**
     * Helper function for uniform and static partitioning. This function
     * calculates the exact match ranges of the parts
     * @param parts the parts of the pattern (they must be set before calling)
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of each part
     * @param counters the performance counters
     */
    void calculateExactMatchRanges(std::vector<Substring>& parts,
                                   std::vector<SARangePair>& exactMatchRanges,
                                   Counters& counters) const;

    // PRIVATE FUNCTIONS: PARTITIONING (UNIFORM)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each part has the
     * same size
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partitionUniform(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PARTITIONING (STATIC)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each search carries
     * the same weight (on average)
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partitionOptimalStatic(const std::string& pattern,
                                std::vector<Substring>& parts,
                                const int& numParts, const int& maxScore,
                                std::vector<SARangePair>& exactMatchRanges,
                                Counters& counters) const;

    /**
     * Helper function for optimal static partitioning. This function
     * creates the optimal static parts
     * @param pattern the pattern to partition
     * @param part  empty vector to which the parts are added
     * @param numParts how many parts there need to be in the partition
     */
    void setParts(const std::string& pattern, std::vector<Substring>& parts,
                  const int& numParts, const int& maxScore) const;

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PARTITIONING (DYNAMIC)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts how many parts are needed
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionDynamic(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: MATCHING
    // ----------------------------------------------------------------------------
    /**
     * Executes the search recursively. If U[0] != 1, then the search will
     * start at pi[0], else the search will start with idx i and U[i]!=0 and
     * U[j]=0 with j < i. The search will not be executed if the U[0]== 0 and
     * exactMatchRanges[pi[0]].width() <= switch point of the index
     * @param s the search to follow
     * @param parts the parts of the pattern
     * @param occ Data structure containing in-index and in-text occurrences.
     * If during the search occurrences are found they will be added to this
     * data structure.
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts
     * @param counters the performance counters
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     Occurrences& occ,
                     const std::vector<SARangePair>& exactMatchRanges,
                     Counters& counters) const;

    /**
     * Matches a sequence using searches with a maximal allowed distance of k
     * Found occurrences either in-text or in-index are added to the occs
     * data structure.
     * Warning: this function assumes that the index is in the correct mode
     * (revCompl or forward strand) and that the in-text verification matrix for
     * this mode has been set.
     * @param seq the sequence to match
     * @param k the maximal allowed edit distance
     * @param counters the performance counters
     * @param occs the data structure containing the occurrences.
     * @param minED the minimal  distance of all considered matches
     */
    void matchWithSearches(const std::string& seq, const length_t k,
                           Counters& counters, Occurrences& occs,
                           const length_t minED = 0);

    /**
     * Maps a read to the index using this strategy, along the correct strand.
     * @param read the read to map (in the correct orientation)
     * @param maxED the maximal allowed edit distance
     * @param counters the performance counters
     * @param pairStatus whether this is the first or second read of a pair
     * @param strand the strand of the read
     * @param minD the minimum allowed distance
     * @returns a vector with approximate occurrences in the text of this read,
     * flagged with the correct strand and pairStatus.
     */
    std::vector<TextOcc> mapRead(const std::string& read, length_t maxED,
                                 Counters& counters, PairStatus pairStatus,
                                 Strand strand, length_t minD = 0) {
        Occurrences occurrences;
        // set the index in the correct mode
        index.setIndexInMode(strand, pairStatus);
        if (maxED == 0) {
            std::vector<TextOcc> occ;
            index.exactMatchesOutput(read, counters, occ);
            return occ;
        }

        // match the sequence and return the filtered occurrences
        matchWithSearches(read, maxED, counters, occurrences, 0);
        // filter out redundant matches
        const auto& bundle = ReadBundle::createBundle(read, strand);
        auto textOcc = (this->*filterPtr)(occurrences, maxED, counters, bundle);

        // remove elements under minD
        textOcc.erase(std::remove_if(textOcc.begin(), textOcc.end(),
                                     [minD](const TextOcc& elem) {
                                         return elem.getDistance() < minD;
                                     }),
                      textOcc.end());
        return textOcc;
    }

    // PRIVATE FUNCTIONS: MATCHING (HAMMING)
    // ----------------------------------------------------------------------------

    /**
     * Starts the index with hamming distance
     * @param s the search to follow
     * @param startMatch the partial match to start from.
     * @param occ Data structure containing in-index and in-text occurrences. If
     * during the search occurrences are found they will be added to this
     * @param parts the parts of the pattern
     * @param idx the index in the search to match next
     */
    void startIndexHamming(const Search& s, const FMOcc& startMatch,
                           Occurrences& occ, std::vector<Substring>& parts,
                           Counters& counters, const int& idx) const {
        index.recApproxMatchHamming(s, startMatch, occ, parts, counters, idx);
    }

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: MATCHING (EDIT)
    // ----------------------------------------------------------------------------

    /**
     * Starts the index with edit distance and optimized alignment for the
     * edit distance metric
     * @param s the search to follow
     * @param startMatch the partial match to start from.
     * @param occ Data structure containing in-index and in-text occurrences. If
     * during the search occurrences are found they will be added to this
     * @param parts the parts of the pattern
     * @param idx the index in the search to match next
     */
    void startIndexEdit(const Search& s, const FMOcc& startMatch,
                        Occurrences& occ, std::vector<Substring>& parts,
                        Counters& counters, const int& idx) const {
        index.recApproxMatchEditEntry(s, startMatch, occ, parts, counters, idx);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END - ALL
    // ----------------------------------------------------------------------------
    /**
     * Matches a pattern approximately using this strategy in, finding all
     * approximate occurrences.
     * @param bundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     *
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     *
     */
    virtual void matchApproxAllMap(ReadBundle& bundle, length_t maxED,
                                   Counters& counters,
                                   std::vector<TextOcc>& result);
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END - BEST (+X)
    // ----------------------------------------------------------------------------

    /**
     * Helper function for matchApproxBestPlusX. Checks if the occurrences are
     * valid and assigns sequences to them if possible
     * @param occVector the vector of occurrences per edit distance
     * @param best the best edit distance found so far, can be updated
     * @param l the current edit distance
     * @param counters the performance counters
     * @param cutOff the maximal allowed edit distance
     * @param seq the sequence to which the occurrences should match, needed for
     * CIGAR string generation and trimming
     */
    void checkAlignments(OccVector& occVector, uint32_t& best, uint32_t l,
                         Counters& counters, uint32_t cutOff,
                         const std::string& seq) const;
    /**
     * Matches a pattern approximately using this strategy in, finding all
     * approximate occurrences.
     * @param readBundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     * @param minIdentity the minimal identity of all considered matches
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     */
    virtual void matchApproxBestPlusX(ReadBundle& readBundle, length_t x,
                                      Counters& counters,
                                      const length_t minIdentity,
                                      std::vector<TextOcc>& result);

    /**
     * Matches a pattern approximately using this strategy, finding all
     * best approximate occurrences.
     * @param readBundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     * @param minIdentity the minimal identity of all considered matches
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     */
    virtual void matchApproxBestMap(ReadBundle& readBundle,
                                    const length_t minIdentity,
                                    Counters& counters,
                                    std::vector<TextOcc>& result) {
        return matchApproxBestPlusX(readBundle, 0, counters, minIdentity,
                                    result);
    }

    // ----------------------------------------------------------------------------
    // PRIVATE FUNCTIONS: BEST
    // ----------------------------------------------------------------------------
    /**
     * Process the sequence.
     * @param seq the sequence to process
     * @param status the status of the read (first or second)
     * @param strand the strand of the read
     * @param maxDist the maximum distance allowed
     * @param vector the vector with the occurrences
     * @param counters the performance counters
     * @returns true if any of the distances between 0 and maxDist have an
     * occurrence
     */
    bool processSeq(const std::string& seq, const PairStatus status,
                    const Strand strand, const length_t maxDist,
                    OccVector& vector, Counters& counters);

    /**
     * Handles trimmed occurrences (whose distance is now larger than before)
     * @param trimmedIdsSet the list with the unique ids of the trimmed
     * occurrences in vector[oDist]
     * @param oDist the original distance of the occurrences
     * @param vector the vector with the occurrences for each distance
     */
    void handleTrimmedOccs(std::vector<length_t>& trimmedIdsSet,
                           const length_t oDist, OccVector& vector);

    /**
     * Helper function for findBestMapping and pairDiscordantlyBest. Maps the
     * stratum if it has not been mapped yet.
     * @param v the vector with occurrences and a boolean flag indicating if the
     * stratum has been mapped
     * @param maxD the maximal distance of this stratum
     * @param counters the performance counters
     * @param seq the sequence to map
     * @param status the status within the pair (first or second)
     * @param strand the strand of the the sequence
     */
    void mapStratum(OccVector& v, length_t maxD, Counters& counters,
                    const std::string& seq, PairStatus status, Strand strand) {
        if (!v[maxD].first) {
            v[maxD].second = mapRead(seq, maxD, counters, status, strand, maxD);
            v[maxD].first = true;
        }
    }

    /**
     * Helper function for pairDiscordantlyBest if no discordant pair could be
     * found. Finds all occurrences in the best stratum for the given read
     * bundle.
     * @param fw the forward occurrences (can be updated). If some strata have
     * already been processed, this should be reflected in this  OccVector.
     * @param rc the reverse complemented occurrences (can be updated). If some
     * strata have already been processed, this should be reflected in this
     * OccVector.
     * @param bundle the read bundle
     * @param counters the performance counters
     * @param status the status within the pair (first or second)
     */
    std::vector<TextOcc> findBestMapping(OccVector& fw, OccVector& rc,
                                         const ReadBundle& bundle,
                                         Counters& counters, PairStatus status);

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PROCESSING SINGLE ENDED OCCS FROM
    // INFERENCE
    // ----------------------------------------------------------------------------

    /**
     * Static function to process the single ended matches in a best-map
     * scenario. This will add the occurrences to the correct OccVector for
     * further pairing. After the function the elements of matches1 and matches2
     * will be in an undefined state.
     * @param matchesSE1 the matches of the first read
     * @param matchesSE2 the matches of the second read
     * @param fw1 the forward matches of the first read (can be updated)
     * @param rc1 the reverse complemented matches of the first read (can be
     * updated)
     * @param fw2 the forward matches of the second read (can be updated)
     * @param rc2 the reverse complemented matches of the second read (can be
     * updated)
     * @param read2done whether the second read has been processed
     */
    static void addSingleEndedForBest(std::vector<TextOcc>& matchesSE1,
                                      std::vector<TextOcc>& matchesSE2,
                                      OccVector& fw1, OccVector& rc1,
                                      OccVector& fw2, OccVector& rc2,
                                      bool read2done);

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: POST PROCESSING
    // ----------------------------------------------------------------------------
    /**
     * Filter the occurrences (both in-index and in-text) based on hamming
     * distance. The in-index occurrences are first converted to in-text
     * occurrences. The function also ensures that each in-text occurrence has a
     * CIGAR string.
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     * @param bundle the read bundle with info about the current read (dummy
     * because filter for edit distance needs it).
     * @returns the filtered in-text occurrences with CIGAR strings
     */
    std::vector<TextOcc> filterHamming(Occurrences& occ, length_t maxED,
                                       Counters& counters,
                                       const ReadBundle& bundle) const {
        return index.getTextOccHamming(occ, counters);
    }

    /**
     * Filter the occurrences (both in-index and in-text) based on edit
     * distance. The in-index occurrences are first converted to in-text
     * occurrences. The function also ensures that each in-text occurrence has a
     * CIGAR string.
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     * @param bundle the read bundle with info about the current read
     * @returns the filtered in-text occurrences ()
     */
    std::vector<TextOcc>
    filterEditWithoutCIGARCalculation(Occurrences& occ, length_t maxED,
                                      Counters& counters,
                                      const ReadBundle& bundle) const {
        return index.getUniqueTextOccurrences(occ, maxED, counters, bundle);
    }

    /**
     * Assigns a sequence to an occurrence.
     * Returns how and if the sequence was found.
     * @param occ the occurrence to assign a sequence to
     * @param counters the performance counters
     * @param maxED the maximal allowed edit distance
     * @param pattern the pattern that was matched
     * @returns how and if the sequence was found
     */
    SeqNameFound assignSequenceAndCIGAR(TextOcc& occ, Counters& counters,
                                        length_t maxED,
                                        const std::string& pattern) const {
        return assignSequence(occ, counters, maxED, pattern);
    }

    /**
     * Assigns a reference sequence to an occurrence if it has not been
     * assigned. The range of the occurrence is also changed to show its
     * position within the assigned sequence
     * @param occ the occurrence to assign a sequence to
     * @param counters the performance counters
     * @param maxED the maximal allowed edit distance
     * @param pattern the pattern that was matched
     * @returns true if a sequence was assigned, false if no sequence could
     * be assigned
     */
    SeqNameFound assignSequence(TextOcc& occ, Counters& counters,
                                length_t maxED,
                                const std::string& pattern) const {

        if (!occ.isSeqNameChecked()) {
            index.setIndexInMode(occ.getStrand(), occ.getPairStatus());

            length_t seqID;
            // if a sequence name is found the coordinates of the occurrence
            // will change
            auto found = index.findSeqName(occ, seqID, counters, maxED,
                                           distanceMetric, pattern);
            if (found != SeqNameFound::NOT_FOUND) {
                occ.setAssignedSequence(found, seqID);
            } else {
                occ.setAssignedSequenceNotFound();
            }
        }

        return occ.getSeqNameFound();
    }

    TextOcc createUnmappedRecordSE(ReadBundle& bundle) const {
        return createUnmappedSEPtr(bundle);
    }

    /**
     * Generates the SAM line for the first occurrence and put the other
     * occurrences in the XA tag.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score
     * @param minScore the best score of the occurrences
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_SAM_XATag(ReadBundle& bundle, size_t nHits, int minScore,
                              std::vector<TextOcc>& occs,
                              Counters& counters) const {
#ifndef NO_REPORT
        occs.front().generateSAMSingleEndXA(bundle, nHits, minScore,
                                            {occs.begin() + 1, occs.end()},
                                            index.getSeqNames());
        counters.inc(Counters::DROPPED_UNIQUE_MATCHES, occs.size() - 1);
        occs.resize(1);
#endif
    }

    /**
     * Generates the SAM lines for all occurrences. The first occurrence is the
     * primary occurrence and the rest are secondary.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score
     * @param minScore the best score of the occurrences
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_SAM(ReadBundle& bundle, size_t nHits, int minScore,
                        std::vector<TextOcc>& occs, Counters& counters) const {
#ifndef NO_REPORT
        occs.front().generateSAMSingleEndFirst(bundle, nHits, minScore,
                                               index.getSeqNames());
        for (length_t i = 1; i < occs.size(); i++) {
            occs[i].generateSAMSingleEndNotFirst(bundle.getSeqID(), nHits,
                                                 minScore, index.getSeqNames());
        }
#endif
    }

    /**
     * Generates the RHS line for the occurrence.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score (dummy parameter)
     * @param minScore the best score of the occurrences (dummy parameter)
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_RHS(ReadBundle& bundle, size_t nHits, int minScore,
                        std::vector<TextOcc>& occs, Counters& counters) const {
        // for each distinct assigned sequence + distance combination keep only
        // the first ones

        // Sort the vector using a lambda comparator
        std::sort(occs.begin(), occs.end(),
                  [](const TextOcc& lhs, const TextOcc& rhs) {
                      if (lhs.getDistance() != rhs.getDistance())
                          return lhs.getDistance() < rhs.getDistance();
                      return lhs.getAssignedSequenceID() <
                             rhs.getAssignedSequenceID();
                  });

        // Use std::unique with a lambda to remove consecutive duplicates
        auto uniqueEnd =
            std::unique(occs.begin(), occs.end(),
                        [](const TextOcc& lhs, const TextOcc& rhs) {
                            return lhs.getDistance() == rhs.getDistance() &&
                                   lhs.getAssignedSequenceID() ==
                                       rhs.getAssignedSequenceID();
                        });

        // Resize the vector to remove the redundant elements
        occs.erase(uniqueEnd, occs.end());

        occs.front().generateRHSSingleEnd(
            bundle, {occs.begin() + 1, occs.end()}, index.getSeqNames());
        counters.inc(Counters::DROPPED_UNIQUE_MATCHES, occs.size() - 1);
        occs.resize(1);
    }

    /**
     * Generates the SAM lines for the occurrences. This function will assign
     * sequences to the occurrences sometimes with the use of trimming (if
     * possible).
     * @param occs the occurrences in the text
     * @param bundle the read bundle with info about the curren read
     * @param counters the performance counters
     * @param cutOff the maximum allowed distance for an occurrence
     */
    void generateOutputSingleEnd(std::vector<TextOcc>& occs, ReadBundle& bundle,
                                 Counters& counters, length_t cutOff) const;

    /**
     * Generates the RHS lines for the occurrences. This function will assign
     * sequences to the occurrences sometimes with the use of trimming (if
     * possible).
     * @param occs the occurrences in the text
     * @param bundle the read bundle with info about the curren read
     * @param counters the performance counters
     * @param cutOff the maximum allowed distance for an occurrence
     */
    void generateRHSSingleEnd(std::vector<TextOcc>& occs, ReadBundle& bundle,
                              Counters& counters, length_t cutOff) const;

    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------

  public:
    // PUBLIC FUNCTIONS
    // ----------------------------------------------------------------------------

    virtual ~SearchStrategy() {
    }

    // PUBLIC FUNCTIONS: ACCESS
    // ----------------------------------------------------------------------------

    /**
     * Retrieves the name of this strategy, derived classes should set a
     * meaningful name
     */
    std::string getName() const {
        return name;
    }
    /**
     * Retrieve the partitioning strategy in string format
     */
    std::string getPartitioningStrategy() const;

    /**
     * Retrieve the distance metric in string format
     */
    std::string getDistanceMetric() const;

    /**
     * Retrieve the mapping in string format
     */
    std::string getMappingModeString() const;

    /**
     * Retrieve the mapping mode
     */
    MappingMode getMappingMode() const {
        return mode;
    }

    /**
     * Calculate the maximal allowed distance for a given minimal identity
     * and sequence size
     * @param minIdentity the minimal identity
     * @param seqSize the size of the sequence
     * @returns the maximal allowed edit distance
     */
    length_t getMaxED(length_t minIdentity, length_t seqSize) const {
        assert(minIdentity >= 50 && minIdentity <= 100);

        length_t cutOff = (seqSize * (100 - minIdentity)) / 100;
        const length_t distanceCutoff = BEST_CUTOFF_COLUMBA;

        return std::min(std::min(distanceCutoff, getMaxSupportedDistance()),
                        cutOff);
    }

    /**
     * Turn on or off the generation of SAM lines
     */
    void setGenerateSAMLines(bool generateSAMLines) {
        this->generateSAMLines = generateSAMLines;
    }

    /**
     * Set whether the output should be SAM format or the custom Read Hit
     * Summary Format
     */
    void setSamOutput(bool samOutput) {
        this->samOutput = samOutput;
        if (samOutput) {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM;
            createUnmappedSEPtr = &TextOcc::createUnmappedSAMOccurrenceSE;
        } else {
            generateOutputSEPtr = &SearchStrategy::generateSE_RHS;
            createUnmappedSEPtr = &TextOcc::createUnmappedRHSOccurrenceSE;
        }
    }

    /**
     * Set the mapping mode to use
     * @param mode the mapping mode to use
     */
    void setMappingMode(MappingMode mode) {
        this->mode = mode;
        if (mode == ALL) {
            matchPtrSE = &SearchStrategy::matchApproxAllMap;
        } else {
            // BEST
            matchPtrSE = &SearchStrategy::matchApproxBestMap;
        }
    }

    /**
     * Set whether to use the XA tag for secondary alignments or use separate
     * SAM lines for them.
     * @param xaTag whether to use the XA tag for secondary alignments
     */
    void setXATag(bool xaTag) {
        this->xaTag = xaTag;
        if (xaTag) {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM_XATag;
        } else {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM;
        }
    }

    /**
     * Set whether unmapped SAM records should not be generated
     * @param unmappedSam whether unmapped SAM records should not be generated
     */
    void setUnmappedSam(bool unmappedSam) {
        this->unmappedSAM = unmappedSam;
    }

    /**
     * Whether the search scheme supports this distance score
     * @param maxScore the maximal score
     * @returns true if the search scheme supports this distance score
     */
    virtual bool supportsDistanceScore(const int& maxScore) const = 0;

    /**
     * Whether the search scheme supports the best mapping mode
     * @param max the maximal score for best mapping mode (output)
     * @returns true if the search scheme supports the best mapping mode
     */
    bool supportsBestMapping(int& max) const {
        max = 0;
        for (length_t i = 0; i < MAX_K; i++) {
            if (supportsDistanceScore(i)) {
                return (max == 0) ? false : true;
            }
            max = i;
        }
        return (max == 0) ? false : true;
    }
    // ----------------------------------------------------------------------------

    // PUBLIC FUNCTIONS: MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Creates all searches for this specific strategy. This is strategy
     * dependent
     * @param maxED the maximal allowed edit distance for the aligning
     * @param exactMatchRanges the ranges for the exact matches of the parts
     */
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& exactMatchRanges) const = 0;

    /**
     * Matches a pattern approximately using this strategy in the current
     * mode. WARNING: make sure that the orientation is set correctly before
     * calling. (@see setOrientation)
     * @param readBundle the read bundle with info about the curren read
     * @param maxEDOrIdentity the maximal allowed edit distance (or  hamming
     * distance) in case of all-map mode and the minimal identity in case of
     * best-map mode
     * @param counters the performance counters
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     *
     */
    virtual void matchApprox(ReadBundle& readBundle, length_t maxEDOrIdentity,
                             Counters& counters, std::vector<TextOcc>& result) {
        (this->*matchPtrSE)(readBundle, maxEDOrIdentity, counters, result);
    }
    // ----------------------------------------------------------------------------

    // PUBLIC FUNCTIONS: PAIRING
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
};

// ============================================================================
// CLASS CUSTOM SEARCHSTRATEGY
// ============================================================================

class CustomSearchStrategy;

// Pointer to the correct getBegins() function for static partitioning
typedef const std::vector<double> (CustomSearchStrategy::*GetBeginsPtr)(
    const int& numParts, const int& maxScore) const;

// Pointer to the correct getSeedingPositions() function for dynamic
// partitioning
typedef const std::vector<double> (
    CustomSearchStrategy::*GetSeedingPositionsPtr)(const int& numParts,
                                                   const int& maxScore) const;
// Pointer to the correct getWeights() function for dynamic partitioning
typedef const std::vector<int> (CustomSearchStrategy::*GetWeightsPtr)(
    const int& numParts, const int& maxScore) const;

/**
 * This is a derived class of SearchStrategy. It creates a custom scheme
 * using files provided by the user. It takes a specified folder in which a
 * file "name.txt", containing the name of the scheme on the first line, and
 * for each supported distance score a subfolder exists. Such a subfolder
 * has as name the distance score. Each subfolder must contain at least a
 * file "searches.txt". The different searches of the scheme for this
 * distance score should be written on separate lines of this file. Each
 * search consists out of three arrays, pi, L and U, the arrays are
 * separated by a single space. Each array is written between curly braces
 * {} and the different values are separated. The pi array must be
 * zero-based.
 *
 * The different sub-folders can also contain files for static and dynamic
 * partitioning. The file "static_partitioning.txt" should consist out of
 * one line of space- separated percentages (between 0 and 1 - exclusive).
 * These percentages point to start positions of the second to the last part
 * of the partitioned pattern (relative to the size of the pattern). Hence,
 * if a search scheme partitions a pattern in k parts, then k+1 percentages
 * should be provided. The file "dynamic_partitioning.txt" should consist
 * out of two lines. The first line contains k-1 space- separated
 * percentages. These percentages are the seeding positions of the middle
 * parts (the first and last part are seeded at the begin and end and thus
 * do not need a percentage). Note that this line can be empty if the
 * pattern is partitioned into 2 parts. The second line contains k integers,
 * where k is the number of parts. Each integer corresponds to the weight
 * given to that part in the dynamic partitioning process.
 */
class CustomSearchStrategy : public SearchStrategy {
  private:
    std::vector<std::vector<Search>>
        schemePerED; // the search schemes for each distance score,
    std::vector<bool>
        supportsMaxScore; // if a particular distance score is supported

    // helpers for static partitioning
    std::vector<std::vector<double>>
        staticPositions;                     // the static positions per k
    std::vector<GetBeginsPtr> beginsPointer; // pointer to the correct
                                             // getBegins() function,
                                             // either default or custom

    // helpers for dynamic partitioning
    std::vector<std::vector<double>>
        seedingPositions; // the seeds for dynamic partitioning per

    std::vector<std::vector<int>>
        weights; // the weights for dynamic partitioning per score

    std::vector<GetSeedingPositionsPtr>
        seedingPointer; // pointer to the correct
                        // getSeedingPositions() function,
                        // either default or custom

    std::vector<GetWeightsPtr>
        weightsPointers; // pointer to the correct getWeights() function,
                         // either default or custom

    /**
     * Retrieves the search scheme from a folder, also checks if the scheme
     * is valid
     * @param pathToFolder the path to the folder containing the search
     * scheme
     * @param verbose if the sanity check should be verbose
     */
    void getSearchSchemeFromFolder(std::string pathToFolder, bool verbose);

    /**
     * If the values provided for dynamic partitioning for the given max
     * score are valid (i.e. strictly increasing and between 0 and 1). Will
     * throw a runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckDynamicPartitioning(const int& maxScore) const;

    /**
     * If the values provided for static partitioning for the given max
     * score are valid (i.e. strictly increasing and between 0 and 1). Will
     * throw a runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckStaticPartitioning(const int& maxScore) const;

    /**
     * Parse the search from a line.
     * @param line the line to parse
     * @param idx the index of the search
     * @returns the parsed line as a search, if the line is not valid a
     * runtime error will be thrown.
     */
    Search makeSearch(const std::string& line, length_t idx) const;

    /**
     * Parses an array from a string.
     * @param vectorString the string to parse
     * @param vector the vector with the parsed array as values
     */
    void getVector(const std::string& vectorString,
                   std::vector<length_t>& vector) const;

    /**
     * Checks whether the connectivity property is satisfied for all
     * searches. Will throw a runtime error if one of these is not satisfied
     * @param verbose if the information about which search covers which
     * error distribution should be written to standard out
     */
    void sanityCheck(bool verbose) const;

    // static partitioning
    /**
     * Gets the static positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double> getBeginsDefault(const int& numParts,
                                               const int& maxScore) const {
        return SearchStrategy::getBegins(numParts, maxScore);
    }

    /**
     * Gets the static positions in the custom manner (i.e. with values
     * provided by user in "static_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double> getBeginsCustom(const int& numParts,
                                              const int& maxScore) const {
        return staticPositions[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the begin positions
     * in the custom way if values were provided in a
     * "static_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*beginsPointer[maxScore - 1])(numParts, maxScore);
    }

    // dynamic partitioning
    /**
     * Gets the seeding positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsDefault(const int& numParts, const int& maxScore) const {
        return SearchStrategy::getSeedingPositions(numParts, maxScore);
    }

    /**
     * Gets the seeding positions in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsCustom(const int& numParts, const int& maxScore) const {
        return seedingPositions[maxScore - 1];
    }

    /**
     * Overridden function of the base class. Retrieves the seeding
     * positions in the custom way if values were provided in a
     * "dynamic_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*seedingPointer[maxScore - 1])(numParts, maxScore);
    }

    /**
     * Gets the seeding positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<int> getWeightsDefault(const int& numParts,
                                             const int& maxScore) const {
        return SearchStrategy::getWeights(numParts, maxScore);
    }

    /**
     * Gets the weights in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<int> getWeightsCustom(const int& numParts,
                                            const int& maxScore) const {
        return weights[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the weights
     * positions in the custom way if values were provided in a
     * "dynamic_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*weightsPointers[maxScore - 1])(numParts, maxScore);
    }

  public:
    CustomSearchStrategy(IndexInterface& index, const std::string& pathToFolder,
                         PartitionStrategy p, DistanceMetric metric,
                         MappingMode mode, SequencingMode sMode,
                         bool verbose = false)
        : SearchStrategy(index, p, metric, mode, sMode) {

        // resize and fill the vectors
        schemePerED.resize(MAX_K);
        supportsMaxScore.resize(MAX_K, false);
        staticPositions.resize(MAX_K);
        beginsPointer.resize(MAX_K, &CustomSearchStrategy::getBeginsDefault);
        seedingPositions.resize(MAX_K);
        weights.resize(MAX_K);
        seedingPointer.resize(
            MAX_K, &CustomSearchStrategy::getSeedingPositionsDefault);
        weightsPointers.resize(MAX_K, &CustomSearchStrategy::getWeightsDefault);

        // read the schemes
        getSearchSchemeFromFolder(pathToFolder, verbose);
    }

    uint32_t calculateNumParts(unsigned int maxED) const override {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1][0].getNumParts();
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges = {}) const override {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1];
    }
    length_t getMaxSupportedDistance() const override {
        // find biggest value for which supportsMaxScore[0: value -1] is
        // true
        length_t max = 0;
        for (length_t i = 0; i < MAX_K; i++) {
            if (supportsMaxScore[i]) {
                max++;
            } else {
                break;
            }
        }
        return max;
    }

    bool supportsDistanceScore(const int& maxScore) const override {
        return supportsMaxScore[maxScore - 1];
    }
};

// ============================================================================
// CLASS MULTIPLE SCHEMES
// ============================================================================

/**
 * @class MultipleSchemes
 *
 * @brief A class representing a collection of search schemes for
 * approximate string matching.
 *
 * This class manages multiple search schemes. It calculates the number of parts
 * and provides a method to choose the search scheme to use based on the number
 * of exact matches of the parts.
 */
class MultipleSchemes {
  private:
    unsigned int k; // the maximum number of errors for these schemes
    std::vector<SearchScheme> schemes; // the list of schemes
    unsigned int numParts = 0; // the number of parts (p) for this scheme

    /**
     * @brief Private helper function to read search schemes from a folder.
     *
     * This function reads search schemes from text files in a specified
     * folder. Each file should be named "scheme1.txt," "scheme2.txt," and
     * so on.
     *
     * @param pathToFolder The path to the folder containing the scheme
     * files.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    void getSchemesFromFolder(const std::string& pathToFolder, bool verbose) {
        // in this folder there should be several files:
        // scheme1.txt, scheme2.txt, scheme3.txt, ...

        std::string base = "scheme";
        int x = 1; // Starting value of x

        while (true) {
            std::string filename = base + std::to_string(x) + ".txt";
            std::string filePath = pathToFolder + PATH_SEPARATOR + filename;

            std::ifstream file(filePath);
            if (!file.is_open()) {
                // File does not exist, break the loop
                break;
            }

            // Read in the scheme
            schemes.emplace_back(SearchScheme::readScheme(file, filename, k));

            // Increment x for the next iteration
            ++x;
        }

        if (!schemes.empty()) {

            numParts = schemes.front().getNumParts();

            for (const auto& scheme : schemes) {
                if (scheme.getNumParts() != numParts) {
                    throw std::runtime_error(
                        "Not all schemes have same amount of parts in: " +
                        pathToFolder);
                }
            }
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemes class.
     *
     * Initializes the MultipleSchemes object by reading search schemes from
     * a folder.
     *
     * @param pathToFolder The path to the folder containing the scheme
     * files.
     * @param k The maximum number of errors for these schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    MultipleSchemes(const std::string& pathToFolder, unsigned int k,
                    bool verbose = false)
        : k(k) {
        getSchemesFromFolder(pathToFolder, verbose);
    }

    /**
     * @brief Constructor for an MultipleSchemes instance.
     *
     * Initializes the MultipleSchemes object by creating an empty list of
     * schemes.
     *
     * @param k The maximum number of errors for these schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    MultipleSchemes(unsigned int k, bool verbose = false) : k(k) {
    }

    /**
     * Create the multiple schemes from a single strategy by mirroring the
     * pi-strings
     * @param strategy the strategy to take the scheme from
     * @param k the maximum number of errors
     */
    MultipleSchemes(const SearchStrategy& strategy, unsigned int k) : k(k) {
        std::vector<SARangePair> dummyRanges;
        const auto& searches = strategy.createSearches(k, dummyRanges);
        SearchScheme scheme(searches, k);
        numParts = scheme.getNumParts();
        SearchScheme reversedScheme = scheme.mirrorPiStrings();
        schemes.push_back(scheme);
        schemes.push_back(reversedScheme);
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit
     * distance. This function should only be called with maxED == k.
     *
     * @param maxED The maximum edit distance for which to calculate the
     * number of parts.
     * @return The number of parts (p) corresponding to the maximum edit
     * distance.
     */
    uint32_t calculateNumParts(unsigned int maxED) const {
        assert(maxED == k);
        return numParts;
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges (=number of exact matches of the searches).
     *
     * This method creates searches based on a given maximum edit distance
     * and a set of SARangePair ranges. It iterates over the schemes and
     * chooses the scheme for which the starting part of the critical search
     * has the least amount of exact matches.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search
     * ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED == k);
        assert(!isEmpty());
        // find the scheme for which the critical search has the smallest
        // exact range to start from

        // if the total number of exact matches of the parts is smaller than
        // the number of parts dynamic selection has too much overhead
        unsigned int total =
            std::accumulate(ranges.begin(), ranges.end(), 0,
                            [](unsigned int sum, const SARangePair& range) {
                                return sum + range.width();
                            });

        if (total <= numParts) {
            return schemes[0].getSearches();
        }

        int minIndex = 0;
        unsigned int minValue =
            ranges[schemes[0].getCriticalPartIndex()].width();

        for (unsigned int i = 1; i < schemes.size(); ++i) {
            unsigned int criticalPartIndex = schemes[i].getCriticalPartIndex();
            if (ranges[criticalPartIndex].width() < minValue) {
                minValue = ranges[criticalPartIndex].width();
                minIndex = i;
            }
        }
        return schemes[minIndex].getSearches();
    }

    /**
     * @brief find out whether the list of schemes is empty
     *
     * @returns true if the list of schemes is empty.
     */
    bool isEmpty() const {
        return schemes.empty();
    }

    /**
     * Add a single scheme to the list of schemes. The scheme to add must be
     * defined for the same number of maximal errors and the same number of
     * parts.
     * @param scheme the scheme to add
     */
    void addScheme(const SearchScheme& scheme) {
        if (scheme.getNumParts() != numParts) {
            throw std::runtime_error(
                "The number of parts of the scheme does not match the number "
                "of parts of the other schemes");
        }
        if (scheme.getK() != k) {
            throw std::runtime_error(
                "The number of errors of the scheme does not match the number "
                "of errors of the other schemes");
        }
        schemes.push_back(scheme);
    }
};
// ============================================================================
// CLASS MULTIPLE SCHEME STRATEGY
// ============================================================================

/**
 * @class MultipleSchemesStrategy
 *
 * @brief A class representing a strategy for handling multiple search
 * schemes with different edit distances.
 *
 * This class is designed to manage multiple search schemes with varying
 * edit distances. It inherits from the SearchStrategy class and provides
 * methods for calculating the number of parts and creating searches based
 * on a specified edit distance and search ranges. The Multiple Schemes will
 * dynamically choose what scheme to use.
 */
class MultipleSchemesStrategy : public SearchStrategy {
  private:
    std::vector<MultipleSchemes>
        schemesPerED; // A vector of MultipleSchemes objects, one for each
                      // edit distance.

    static bool directoryExists(const std::string& path) {
        struct stat info;
        return (stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR));
    }

    static bool fileExists(const std::string& path) {
        std::ifstream file(path);
        return file.good();
    }

    /**
     * @brief Private helper function to read and initialize multiple search
     * schemes.
     *
     * This function reads and initializes multiple search schemes from a
     * specified folder for various edit distances. If an edit distance is
     * not provided the function stops reading.
     *
     * @param pathToFolder The path to the folder containing the search
     * schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    void readSchemes(const std::string& pathToFolder, bool verbose) {

        // get the name of the scheme
        std::string line;
        {
            std::ifstream ifs(pathToFolder + "name.txt");
            if (!ifs) {
                throw std::runtime_error(
                    "Problem reading: " + pathToFolder +
                    "name.txt\nDid you provide a directory to "
                    "a search scheme without a name file?");
            }
            std::getline(ifs, line);
            name = line;
            ifs.close();
        }

        // read in the schemes
        for (uint16_t k = 1; k <= MAX_K; k++) {
            std::string dir = pathToFolder + PATH_SEPARATOR + std::to_string(k);

            if (directoryExists(pathToFolder)) {
                std::string file = dir + PATH_SEPARATOR + "scheme1.txt";
                if (fileExists(file)) {
                    schemesPerED.emplace_back(MultipleSchemes(dir, k, verbose));
                    continue;
                }
            }
            // Add an empty search scheme
            schemesPerED.emplace_back(MultipleSchemes(k, verbose));
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemesStrategy class.
     *
     * Initializes the MultipleSchemesStrategy object by reading and
     * initializing search schemes for various edit distances.
     *
     * @param index The IndexInterface object used for searching.
     * @param pathToFolder The path to the folder containing the search
     * schemes.
     * @param p The partition strategy to use
     * @param metric The distance metric to use
     * @param mode The mapping mode to use
     * @param sMode The sequencing mode to use
     * @param verbose A flag indicating whether to display verbose output
     * during reading (default: false).
     */
    MultipleSchemesStrategy(IndexInterface& index,
                            const std::string& pathToFolder,
                            PartitionStrategy p, DistanceMetric metric,
                            MappingMode mode, SequencingMode sMode,
                            bool verbose = false)
        : SearchStrategy(index, p, metric, mode, sMode) {

        readSchemes(pathToFolder, verbose);
    }

    /**
     * Construct a multiple schemes strategy from a single strategy. This will
     * take the schemes in the strategy and mirror them and add both versions to
     * the multiple schemes.
     * @param strategy the strategy to take the scheme from
     * @param verbose if the constructor should be verbose
     */
    MultipleSchemesStrategy(SearchStrategy* strategy, bool verbose = false)
        : SearchStrategy(*strategy) {
        schemesPerED.reserve(MAX_K);
        for (uint16_t k = 1; k <= MAX_K; k++) {
            if (strategy->supportsDistanceScore(k)) {
                // create the multiple schemes
                schemesPerED.emplace_back(MultipleSchemes(*strategy, k));
            } else {
                // add an empty scheme
                schemesPerED.emplace_back(MultipleSchemes(k, verbose));
            }
        }
    }

    void addScheme(const SearchScheme& scheme) {
        unsigned int k = scheme.getK();
        if (k > MAX_K) {
            throw std::runtime_error("The maximum number of errors is " +
                                     std::to_string(MAX_K));
        }
        if (schemesPerED.size() < k) {
            for (uint16_t i = schemesPerED.size(); i < k; i++) {
                schemesPerED.emplace_back(MultipleSchemes(i + 1));
            }
        }
        schemesPerED[k - 1].addScheme(scheme);
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit
     * distance.
     *
     * @param maxED The maximum edit distance for which to calculate the
     * number of parts.
     * @return The number of parts (p) corresponding to the specified
     * maximum edit distance.
     */
    uint32_t calculateNumParts(unsigned int maxED) const override {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].calculateNumParts(maxED);
    }

    length_t getMaxSupportedDistance() const override {
        return schemesPerED.size();
    }
    bool supportsDistanceScore(const int& maxScore) const override {
        return (length_t)maxScore - 1 < schemesPerED.size() &&
               !schemesPerED[maxScore - 1].isEmpty();
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges.
     *
     * This method creates searches based on a given maximum edit distance
     * and a set of SARangePair ranges.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search
     * ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].createSearches(maxED, ranges);
    }
};

// ============================================================================
// CLASS NaiveBackTrackingStrategy
// ============================================================================

/**
 * Search strategy that uses naive backtracking for approximate string matching.
 */
class NaiveBackTrackingStrategy : public SearchStrategy {
  private:
    const std::vector<std::vector<Search>> searches = createSearchesVector();
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {

        return searches[maxED - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return MAX_K;
    }
    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore <= (int)getMaxSupportedDistance();
    }

    static std::vector<std::vector<Search>> createSearchesVector() {
        std::vector<std::vector<Search>> searches;
        searches.reserve(MAX_K - 1);
        for (length_t k = 1; k < MAX_K; k++)
            searches.push_back({Search::makeSearch({0}, {0}, {k}, 0)});
        return searches;
    }

  public:
    NaiveBackTrackingStrategy(IndexInterface& index, PartitionStrategy p,
                              DistanceMetric metric, MappingMode mode,
                              SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "Naive backtracking";
    }
};

// ============================================================================
// HARDCODED CUSTOM CLASSES
// ============================================================================

/**
 * Hardcoded class for the Kucherov  k+1 strategy.
 */
class KucherovKPlus1 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 1}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}, 0),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 1, 2}, 1),
        Search::makeSearch({1, 0, 2}, {0, 0, 1}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4},
                           0),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4},
                           1),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 3, 4},
                           2),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 3, 4},
                           3),
        Search::makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4},
                           4),
        Search::makeSearch({2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4},
                           5),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4},
                           6),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4},
                           7)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.57}, {0.38, 0.65}, {0.38, 0.55, 0.73}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {39, 10, 40}, {400, 4, 5, 400}, {100, 5, 1, 6, 105}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.41, 0.7}, {0.25, 0.50, 0.75}, {0.27, 0.47, 0.62, 0.81}};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED >= 1);
        assert(maxED <= 4);

        return schemePerED[maxED - 1];
    }
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    KucherovKPlus1(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode mode,
                   SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "KUCHEROV K + 1";
    };
    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

/**
 * Hardcoded class for the Kucherov  k+2 strategy.
 */
class KucherovKPlus2 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}, 0),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}, 1),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}, 2),
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2}, 3)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3},
                           3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}, 1),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 2, 4, 4}, 2),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2},
                           {0, 1, 1, 3, 4, 4}, 3),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3},
                           {0, 1, 1, 2, 4, 4}, 4),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3},
                           {0, 0, 3, 3, 4, 4}, 5),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3},
                           {0, 0, 3, 3, 4, 4}, 6),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 2, 4, 4, 4}, 7),
        Search::makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4},
                           {0, 0, 2, 2, 4, 4}, 8),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 1, 4, 4, 4}, 9)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.48, 0.55}, {0.4, 0.63, 0.9}, {0.34, 0.5, 0.65, 0.7}};

    const std::vector<std::vector<int>> weights = {{11, 10, 1},
                                                   {400, 4, 1, 800},
                                                   {6, 3, 2, 1, 1},
                                                   {52, 42, 16, 14, 1, 800}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.47, 0.94},
        {0.35, 0.50, 0.65},
        {0.22, 0.44, 0.66, 0.88},
        {0.18, 0.37, 0.53, 0.69, 0.83}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    KucherovKPlus2(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode mode,
                   SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "KUCHEROV K + 2";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

/**
 * Hardcoded class for the optimal strategy by Kianfar et al.
 */
class OptimalKianfar : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 1}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}, 0),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}, 1),
        Search::makeSearch({1, 2, 0}, {0, 1, 1}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}, 0),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3}, 2)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 4}, {0, 3, 3, 4, 4},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {2, 2, 3, 3, 4},
                           1),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 3, 3}, {0, 0, 4, 4, 4},
                           2)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.50}, {0.34, 0.66}, {0.42, 0.56, 0.67}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {10, 1, 5}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.30, 0.60}, {0.17, 0.69, 0.96}, {0.2, 0.5, 0.6, 0.8}};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        if (maxED < 1 || maxED > 5) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        return schemePerED[maxED - 1];
    }

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    OptimalKianfar(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode m, SequencingMode sMode)
        : SearchStrategy(index, p, metric, m, sMode) {
        name = "OPTIMAL KIANFAR";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

/**
 * A concrete derived class of SearchStrategy.The strategy here is founded
 * on this observation: if x errors are allowed and the pattern is divided
 * up in (x + 2) parts then every match with max x errors contains a seed
 * consisting of n parts, where the first and last part of the seed contain no
 * errors and all parts between these contain exactly one error.(2 <= n <= x +
 * 2) */
class O1StarSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}, 0),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3},
                           3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 1),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 2),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 3),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 0, 4, 4, 4, 4}, 4),
    };

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.51, 0.93}, {0.34, 0.64, 0.88}, {0.28, 0.48, 0.63, 0.94}};

    const std::vector<std::vector<int>> weights = {
        {11, 10, 1}, {20, 11, 11, 10}, {3, 2, 2, 1, 1}, {1, 2, 2, 1, 2, 1}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.50, 0.96},
        {0.26, 0.64, 0.83},
        {0.22, 0.46, 0.67, 0.95},
        {0.19, 0.37, 0.57, 0.74, 0.96}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    O1StarSearchStrategy(IndexInterface& index, PartitionStrategy p,
                         DistanceMetric metric, MappingMode mode,
                         SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "01*0";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS PIGEONHOLE SEARCHSTRATEGY
// ============================================================================

/**
 * A concrete derived class of SearchStrategy. The strategy here is founded
 * on this observation: if x errors are allowed and the pattern is divided
 * up in (x
 * + 1) sections then every approximate match has an exact match with at
 * least one of the sections. The strategy iterates over the sections, it
 * tries to exactly match the current section, then approximately match the
 * pattern before this section and after the the pattern after this section
 * with the remaining edit distance.
 */
class PigeonHoleSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 2, 2}, 1),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           3),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           4)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    PigeonHoleSearchStrategy(IndexInterface& index, PartitionStrategy p,
                             DistanceMetric metric, MappingMode mode,
                             SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "PIGEON HOLE";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS MINU SEARCHSTRATEGY
// ============================================================================

/**
 * A concrete derived class of SearchStrategy. This search strategy contains
 * the search schemes created by Hato.
 */
class MinUSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 1, 1}, {0, 2, 2}, 0),
        Search::makeSearch({1, 0, 2}, {0, 0, 0}, {0, 1, 2}, 1),
        Search::makeSearch({2, 1, 0}, {0, 0, 2}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 1, 1, 1}, {0, 1, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 2}, {0, 1, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 1, 1, 3}, {0, 1, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 2, 2, 2}, {0, 2, 2, 4, 4},
                           0),
        Search::makeSearch({1, 2, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 4, 4},
                           1),
        Search::makeSearch({2, 1, 0, 3, 4}, {0, 1, 1, 1, 1}, {0, 1, 2, 4, 4},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 3}, {0, 1, 4, 4, 4},
                           3),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 1, 1, 1, 4}, {0, 1, 4, 4, 4},
                           4)};

    const std::vector<Search> ED5 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 2, 2, 2},
                           {0, 1, 3, 5, 5, 5}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5}, {0, 1, 1, 3, 3, 3},
                           {0, 1, 3, 5, 5, 5}, 1),
        Search::makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 5, 5}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5}, {0, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 5, 5}, 3),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 4},
                           {0, 1, 3, 5, 5, 5}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 5},
                           {0, 1, 3, 5, 5, 5}, 5)};

    const std::vector<Search> ED6 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6}, {0, 0, 2, 2, 2, 2, 6},
                           {0, 2, 2, 6, 6, 6, 6}, 0),
        Search::makeSearch({1, 2, 0, 3, 4, 5, 6}, {0, 1, 1, 1, 1, 1, 5},
                           {0, 1, 2, 6, 6, 6, 6}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0, 4},
                           {0, 1, 2, 6, 6, 6, 6}, 2),
        Search::makeSearch({3, 4, 5, 6, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 6, 6, 6}, 3),
        Search::makeSearch({4, 3, 5, 6, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 6, 6, 6}, 4),
        Search::makeSearch({5, 6, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2},
                           {0, 1, 3, 3, 6, 6, 6}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3},
                           {0, 1, 3, 3, 6, 6, 6}, 6)};
    const std::vector<Search> ED7 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7}, {0, 1, 1, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 1),
        Search::makeSearch({2, 3, 1, 0, 4, 5, 6, 7}, {0, 0, 0, 2, 2, 2, 2, 2},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7}, {0, 1, 1, 3, 3, 3, 3, 3},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 3),
        Search::makeSearch({4, 5, 6, 7, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 4},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 4),
        Search::makeSearch({5, 4, 6, 7, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1, 5},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 5),
        Search::makeSearch({6, 7, 5, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2, 6},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3, 7},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 7)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4,
                                                          ED5, ED6, ED7};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }

    length_t getMaxSupportedDistance() const override {
        return 7;
    }

  protected:
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED <= 7);
        assert(maxED >= 1);
        return schemePerED[maxED - 1];
    }

  public:
    MinUSearchStrategy(IndexInterface& index, PartitionStrategy p,
                       DistanceMetric metric, MappingMode mode,
                       SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "minU";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 7;
    }
};

/**
 * The default search strategy for Columba without dynamic selection. It uses
 * the MinU search schemes for k <=7 and the greedy search schemes for k > 7 and
 * k <= 13.
 */
class ColumbaSearchStrategy : public MinUSearchStrategy {
  public:
    ColumbaSearchStrategy(IndexInterface& index, PartitionStrategy p,
                          DistanceMetric metric, MappingMode mode,
                          SequencingMode sMode)
        : MinUSearchStrategy(index, p, metric, mode, sMode) {
        name = "Columba";
    };

    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        if (maxED <= 7) {
            return MinUSearchStrategy::createSearches(maxED, ranges);
        } else {
            return schemePerEDHigh[maxED - 8];
        }
    }

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 13;
    }

  private:
    const std::vector<Search> ED8 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 2, 3, 5, 7, 8, 8, 8, 8}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 4, 5, 7, 8, 8, 8, 8}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 6, 7, 8, 8, 8, 8}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5},
                           {0, 1, 3, 3, 7, 8, 8, 8, 8}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6},
                           {0, 1, 3, 4, 4, 8, 8, 8, 8}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 5, 5, 5, 8, 8, 8}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 5, 6, 6, 6, 8, 8}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 7, 7, 7, 7, 8}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8}, 8),
    };

    const std::vector<Search> ED9 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 9}, 1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 7, 7, 7, 7, 9, 9}, 2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 6, 6, 6, 6, 9, 9, 9}, 3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 5, 5, 5, 9, 9, 9, 9}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 3, 4, 5, 6, 7, 8},
                           {0, 1, 4, 4, 4, 9, 9, 9, 9, 9}, 5),
        Search::makeSearch({6, 7, 8, 9, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 3, 3, 7, 9, 9, 9, 9, 9}, 6),
        Search::makeSearch({7, 8, 9, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 2, 2, 5, 7, 9, 9, 9, 9, 9}, 7),
        Search::makeSearch({8, 9, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6},
                           {0, 1, 3, 5, 7, 9, 9, 9, 9, 9}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 1, 2, 2, 3, 4, 5, 6, 7},
                           {0, 1, 3, 5, 7, 9, 9, 9, 9, 9}, 9),
    };

    const std::vector<Search> ED10 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 2, 3, 5, 7, 9, 10, 10, 10, 10, 10}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 4, 5, 7, 9, 10, 10, 10, 10, 10}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 6, 7, 9, 10, 10, 10, 10, 10}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 3, 3, 7, 9, 10, 10, 10, 10, 10}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 3, 4, 4, 9, 10, 10, 10, 10, 10}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 5, 5, 5, 10, 10, 10, 10, 10}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 5, 6, 6, 6, 10, 10, 10, 10}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 7, 7, 7, 7, 10, 10, 10}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 10, 10}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 10}, 9),
        Search::makeSearch({10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10}, 10),
    };

    const std::vector<Search> ED11 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 11}, 1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 11, 11}, 2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 11, 11, 11}, 3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 10, 11, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 7, 7, 7, 7, 11, 11, 11, 11}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 10, 11, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 5, 6, 6, 6, 11, 11, 11, 11, 11}, 5),
        Search::makeSearch({6, 7, 8, 9, 10, 11, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 5, 5, 5, 11, 11, 11, 11, 11, 11}, 6),
        Search::makeSearch({7, 8, 9, 10, 11, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 3, 4, 4, 9, 11, 11, 11, 11, 11, 11}, 7),
        Search::makeSearch({8, 9, 10, 11, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 3, 3, 7, 9, 11, 11, 11, 11, 11, 11}, 8),
        Search::makeSearch({9, 10, 11, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 6, 7, 9, 11, 11, 11, 11, 11, 11}, 9),
        Search::makeSearch({10, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 4, 5, 7, 9, 11, 11, 11, 11, 11, 11}, 10),
        Search::makeSearch({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 2, 3, 5, 7, 9, 11, 11, 11, 11, 11, 11}, 11)};
    const std::vector<Search> ED12 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 2, 3, 5, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 1, 4, 5, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 6, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 3, 3, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 3, 4, 4, 9, 11, 12, 12, 12, 12, 12, 12}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 5, 5, 5, 11, 12, 12, 12, 12, 12, 12}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 5, 6, 6, 6, 12, 12, 12, 12, 12, 12}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 7, 7, 7, 7, 12, 12, 12, 12, 12}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 12, 12, 12, 12}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 12, 12, 12}, 9),
        Search::makeSearch({10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 12, 12}, 10),
        Search::makeSearch({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11, 12}, 11),
        Search::makeSearch({12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 5, 12, 12, 12, 12, 12, 12, 12}, 12)};

    const std::vector<Search> ED13 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 6, 13, 13, 13, 13, 13, 13, 13},
                           0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 5, 12, 12, 12, 12, 12, 12, 12, 13},
                           1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11, 13, 13},
                           2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 13, 13, 13},
                           3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 13, 13, 13, 13}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 10, 11, 12, 13, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 13, 13, 13, 13, 13}, 5),
        Search::makeSearch({6, 7, 8, 9, 10, 11, 12, 13, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 3, 7, 7, 7, 7, 13, 13, 13, 13, 13, 13}, 6),
        Search::makeSearch({7, 8, 9, 10, 11, 12, 13, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 5, 6, 6, 6, 13, 13, 13, 13, 13, 13, 13},
                           7),
        Search::makeSearch({8, 9, 10, 11, 12, 13, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 5, 5, 5, 11, 13, 13, 13, 13, 13, 13, 13},
                           8),
        Search::makeSearch({9, 10, 11, 12, 13, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 1, 3, 4, 4, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           9),
        Search::makeSearch({10, 11, 12, 13, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 3, 3, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           10),
        Search::makeSearch({11, 12, 13, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 2, 6, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           11),
        Search::makeSearch({12, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 1, 4, 5, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           12),
        Search::makeSearch({13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                           {0, 2, 3, 5, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           13)};

    const std::vector<std::vector<Search>> schemePerEDHigh = {ED8,  ED9,  ED10,
                                                              ED11, ED12, ED13};

    length_t getMaxSupportedDistance() const override {
        return 13;
    }
};

/**
 * A search strategy based on the minU seach schemes, and the greedy solution
 * seach schemes that uses dynamic selection. This is used by default in
 * Columba. It is based on ColumbaSearchStrategy but adds the search schemes
 * that have the critical search start with a different part.
 */
class DynamicColumbaStrategy : public MultipleSchemesStrategy {

  private:
    const static std::vector<Search> getMidSearch2() {
        return {Search::makeSearch({2, 1, 0}, {0, 1, 1}, {0, 2, 2}, 0),
                Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 1, 2}, 1),
                Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}, 2)};
    }
    // keep the middle search schemes for even k

    const static std::vector<Search> getMidSearch4() {
        return {Search::makeSearch({0, 1, 2, 3, 4}, {0, 1, 1, 1, 4},
                                   {0, 1, 4, 4, 4}, 0),
                Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 0, 0, 3},
                                   {0, 1, 4, 4, 4}, 1),
                Search::makeSearch({2, 3, 4, 1, 0}, {0, 1, 1, 1, 1},
                                   {0, 2, 2, 4, 4}, 2),
                Search::makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 0, 0},
                                   {0, 1, 2, 4, 4}, 3),
                Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 2, 2, 2},
                                   {0, 1, 2, 4, 4}, 4)};
    }

    const static std::vector<Search> getMidSearch6() {
        return {Search::makeSearch({0, 1, 2, 3, 4, 5, 6}, {0, 1, 1, 1, 1, 1, 5},
                                   {0, 1, 2, 6, 6, 6, 6}, 0),
                Search::makeSearch({1, 0, 2, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0, 4},
                                   {0, 1, 2, 6, 6, 6, 6}, 1),
                Search::makeSearch({2, 1, 0, 3, 4, 5, 6}, {0, 0, 2, 2, 2, 2, 6},
                                   {0, 2, 2, 6, 6, 6, 6}, 2),
                Search::makeSearch({3, 4, 5, 6, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2},
                                   {0, 1, 3, 3, 6, 6, 6}, 3),
                Search::makeSearch({4, 3, 5, 6, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3},
                                   {0, 1, 3, 3, 6, 6, 6}, 4),
                Search::makeSearch({5, 6, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0},
                                   {0, 1, 3, 3, 6, 6, 6}, 5),
                Search::makeSearch({6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1},
                                   {0, 1, 3, 3, 6, 6, 6}, 6)};
    }

    /**
     * @brief Create a new dynamic Columba strategy.
     */
    DynamicColumbaStrategy(SearchStrategy* strategy)
        : MultipleSchemesStrategy(strategy) {
    }

  public:
    /**
     * Statically create a pointer to dynamic columba strategy.
     * @param index the index to use
     * @param p the partition strategy to use
     * @param metric the distance metric to use
     * @param mode the mapping mode to use
     * @param sMode the sequencing mode to use
     */
    static DynamicColumbaStrategy*
    createDynamicColumbaStrategy(IndexInterface& index, PartitionStrategy p,
                                 DistanceMetric metric, MappingMode mode,
                                 SequencingMode sMode) {
        ColumbaSearchStrategy columbaStrategy(index, p, metric, mode, sMode);
        auto instance = new DynamicColumbaStrategy(&columbaStrategy);

        std::vector<Search> midSearch2 = {
            Search::makeSearch({2, 1, 0}, {0, 1, 1}, {0, 2, 2}, 0),
            Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 1, 2}, 1),
            Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}, 2)};
        instance->addScheme(SearchScheme(getMidSearch2(), 2));
        instance->addScheme(SearchScheme(getMidSearch4(), 4));
        SearchScheme scheme6 = SearchScheme(getMidSearch6(), 6);
        instance->addScheme(scheme6);
        instance->addScheme(scheme6.mirrorPiStrings());
        return instance;
    }
};

/**
 * A search strategy based on a custom search strategy provided by the user but
 * with dynamic selection between the schemes and their mirrored variants.
 * @see CustomSearchStrategy for the requirements of the definition of the
 * custom search strategy.
 */
class DynamicCustomStrategy : public MultipleSchemesStrategy {
  private:
    /**
     * @brief Create a new dynamic custom strategy.
     */
    DynamicCustomStrategy(SearchStrategy* strategy)
        : MultipleSchemesStrategy(strategy) {
    }

  public:
    /**
     * Statically create a pointer to the dynamic variant of a custom search
     * strategy. @see CustomSearchStrategy for the requirements in the folder.
     * @param index the index to use
     * @param pathToFolder the path to the folder where the custom strategy is
     * defined
     * @param p the partition strategy to use
     * @param metric the distance metric to use
     * @param mode the mapping mode to use
     * @param sMode the sequencing mode to use
     * @param verbose whether to print debug information
     */
    static DynamicCustomStrategy* createDynamicCustomSearchStrategy(
        IndexInterface& index, const std::string& pathToFolder,
        PartitionStrategy p, DistanceMetric metric, MappingMode mode,
        SequencingMode sMode, bool verbose = false) {
        CustomSearchStrategy customStrategy(index, pathToFolder, p, metric,
                                            mode, sMode, verbose);
        auto instance = new DynamicCustomStrategy(&customStrategy);
        return instance;
    }
};

#endif