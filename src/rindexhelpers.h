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

#ifndef INDEXHELPERS_H
#define INDEXHELPERS_H

#include "definitions.h" // for length_t, PairStatus, Strand, FIRST_IN...
#include "reads.h"       // for ReadBundle

#include <algorithm>          // for max, copy, sort, unique
#include <array>              // for array
#include <assert.h>           // for assert
#include <chrono>             // for duration, operator-, high_resolution_clock
#include <cmath>              // for log10, round
#include <cstdint>            // for uint16_t, uint32_t, uint64_t, int16_t
#include <ext/alloc_traits.h> // for __alloc_traits<>::value_type
#include <fmt/core.h>         // for format
#include <fmt/format.h>       // for to_string
#include <memory>             // for allocator, allocator_traits<>::value_type
#include <sstream>            // for ostream
#include <stddef.h>           // for size_t
#include <string>             // for string, operator+, char_traits, basic_...
#include <utility>            // for move, pair
#include <vector>             // for vector

#ifdef LF_BENCHMARK_FUNCTIONALITY
extern length_t LF_call_count;
extern length_t elapsed_LF;
#endif

#ifdef PHI_BENCHMARK_FUNCTIONALITY
extern length_t phi_call_count;
extern length_t locate_call_count;
extern length_t elapsed_phi;
extern std::chrono::duration<double> elapsed_locate;
#endif


// ============================================================================
// HELPERS
// ============================================================================

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

// ============================================================================
// CLASS RANGE
// ============================================================================

/**
 * A range. A range is a pair of integers (begin, end) where begin is the
 * beginning of the range and end is the end of the range (non-inclusive).
 */
class Range {
  protected:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * Constructor
     * @param b the beginning of the range
     * @param e  the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    Range() : begin(0), end(0) {
    }

    /**
     * @returns the beginning of the range
     */
    length_t getBegin() const {
        return begin;
    }

    /**
     * @returns the end of the range (non-inclusive)
     */
    length_t getEnd() const {
        return end;
    }
    /**
     * Check if this range is empty.
     * @returns true if the range is empty, false otherwise
     */
    bool empty() const {
        return end <= begin;
    }

    /**
     * Gets the width of the range (end - begin)
     * @returns the width of this range
     */
    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    /**
     * Operator overloading, two ranges are equal if their begin and end field
     * are equal.
     * @param o the range to compare to this
     * @returns true if this is equal to o
     */
    bool operator==(const Range& o) const {
        return o.getBegin() == begin && o.getEnd() == end;
    }
    /**
     * Operator overloading, to stream the range to an output stream.
     * @param os the output stream
     * @param r the range to print
     */
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
};

/**
 * Operator overloading. Outputs the range as [begin, end) to the output stream
 * @param output The output stream.
 * @param r The range to print.
 */
std::ostream& operator<<(std::ostream& output, const Range& r);

// ============================================================================
// CLASS MOVERANGE
// ============================================================================

class MoveRange : public Range {
  private:
    length_t beginRun; // run index corresponding to the beginning of the range
    length_t
        endRun; // run index corresponding to the end of the range (inclusive)
    bool runIndicesValid;

  public:
    /**
     * Constructor
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     * @param bRun, the run index corresponding to the beginning of the range
     * @param eRun, the run index corresponding to the end of the range
     * (inclusive)
     * @param runIndicesValid, true if the run indices are valid, false
     * otherwise
     */
    MoveRange(length_t b, length_t e, length_t bRun, length_t eRun,
              bool runIndicesValid = true)
        : Range(b, e), beginRun(bRun), endRun(eRun),
          runIndicesValid(runIndicesValid) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    MoveRange() : Range(), beginRun(0), endRun(0), runIndicesValid(false) {
    }

    /**
     * @returns the run index corresponding to the beginning of the range
     */
    length_t getBeginRun() const {
        return beginRun;
    }

    /**
     * @returns the run index corresponding to the end of the range (inclusive)
     */
    length_t getEndRun() const {
        return endRun;
    }

    /**
     * Sets the run index corresponding to the beginning of the range
     * @param bRun the run index to set
     */
    void setBeginRun(length_t bRun) {
        beginRun = bRun;
    }

    /**
     * Sets the run index corresponding to the end of the range
     * @param eRun the run index to set
     */
    void setEndRun(length_t eRun) {
        endRun = eRun;
    }

    /**
     * @brief Set the Begin object
     *
     * @param b
     */
    void setBegin(length_t b) {
        begin = b;
    }

    /**
     * @brief Set the End object
     *
     * @param e
     */
    void setEnd(length_t e) {
        end = e;
    }

    /**
     * @brief Set the range to empty values
     *
     */
    void setEmpty() {
        begin = 0;
        end = 0;
        beginRun = 0;
        endRun = 0;
        runIndicesValid = false;
    }

    /**
     * @returns true if the run indices are valid, false otherwise
     */
    bool getRunIndicesValid() const {
        return runIndicesValid;
    }

    /**
     * Sets the run indices to valid or invalid
     * @param valid true if the run indices are valid, false otherwise
     */
    void setRunIndicesValid(bool valid) {
        this->runIndicesValid = valid;
    }

    /**
     * Operator overloading, two MoveRanges are equal if their begin, end,
     * beginRun and endRun fields are equal.
     * @param o the MoveRange to compare to this
     * @returns true if this is equal to o
     */
    bool operator==(const MoveRange& o) const {
        return o.getBegin() == begin && o.getEnd() == end &&
               o.getBeginRun() == beginRun && o.getEndRun() == endRun;
    }
    /**
     * Operator overloading, to stream the MoveRange to an output stream.
     * @param os the output stream
     * @param r the MoveRange to print
     */
    friend std::ostream& operator<<(std::ostream& os, const MoveRange& r);
};

/**
 * Operator overloading. Outputs the MoveRange to the output stream
 * @param output The output stream.
 * @param r The MoveRange to print.
 */
std::ostream& operator<<(std::ostream& output, const MoveRange& r);

// ============================================================================
// SARANGE DEFINITION
// ============================================================================

typedef MoveRange SARange; // the range in the suffix array in the move table

// ============================================================================
// CLASS TextOccurrence
// ============================================================================

/**
 * An occurrence in the text. An occurrence is a range in the text and a
 * distance to the mapped read. The distance can be either an edit distance or
 * a hamming distance. The occurrence can also contain a CIGAR string, an
 * assigned sequence, a line in SAM format and a flag representing whether this
 * occurrence is along the forward or reverse complemented strand.
 */
class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)

    std::string stringCIGAR = "*"; // Set to *

    length_t assignedSequenceID = -1; // the ID of the assigned sequence

    bool seqNameChecked = false; // indicates whether the sequence name was
    // found for this occurrence
    SeqNameFound seqNameFound = NOT_FOUND; // indicates whether

    Strand strand = FORWARD_STRAND; // the strand on which this occurrence lies
    PairStatus pairStatus = FIRST_IN_PAIR; // indicates whether this occurrence
    // is the first read in a pair

    std::string outputLine = ""; // the output line for this text Occ

    length_t indexBegin; // the position in the indexed text where the
                         // occurrence starts

    /**
     * Get the flags for the SAM line of this occurrence (single-ended).
     * @param primaryAlignment Indicates whether this is the primary alignment
     * @returns the flags for the SAM line of this occurrence.
     */
    uint16_t getFlagsSE(bool primaryAlignment) const {
        uint16_t result = 0;

        // set rev complement flag
        result |= (isRevCompl() << 4);

        // set secondary alignment flag (256)
        result |= (!primaryAlignment << 8);

        return result;
    }

    /**
     * Calculates the mapping quality based on the score, the number of other
     * matches/pairs and the score of the best match/pair.
     * @param hitNumbers the number of other matches/pairs
     * @param minScore the score of the best match/pair
     */
    int16_t getMapQ(uint32_t hitNumbers, uint32_t minScore) const {
        if (distance != minScore) {
            return 0;
        }
        assert(hitNumbers > 0);
        if (hitNumbers == 1) {
            return MAX_MAPQ;
        } else {
            return round(-10.0 * log10(1 - 1.0 / hitNumbers));
        }
    }

    /**
     * Gets the string that represents this occurrence in the XA tag of the SAM
     * format. XA tag format: rname,pos,CIGAR,NM; The sign of pos indicates the
     * strand
     * @param seqNames the vector with all sequence names
     */
    std::string asXA(const std::vector<std::string>& seqNames) const {
        std::string sign = (isRevCompl()) ? "-" : "+";
        return fmt::format("{},{}{},{},{};", seqNames[assignedSequenceID], sign,
                           range.getBegin() + 1, // SAM is 1-based
                           stringCIGAR, distance);
    }

  public:
    /**
     * Constructor
     * @param range The range of this occurrence in the text.
     * @param distance The (edit or hamming) distance to the mapped read of
     * this occurrence.
     * @param strand The strand on which this occurrence lies.
     * @param pairStatus Indicates whether this occurrence is the first or
     * second read in the pair.
     */
    TextOcc(Range range, length_t distance, Strand strand,
            PairStatus pairStatus)
        : range(range), distance(distance), strand(strand),
          pairStatus(pairStatus), indexBegin(range.getBegin()) {
    }

    /**
     * Constructor
     * @param range The range of this occurrence in the text.
     * @param distance The (edit or hamming) distance to the mapped read of
     * this occurrence.
     * @param CIGAR The CIGAR string of this occurrence.
     * @param strand The strand on which this occurrence lies.
     * @param pairStatus Indicates whether this occurrence is the first or
     * second read in the pair.
     */
    TextOcc(Range range, length_t distance, const std::string& CIGAR,
            Strand strand, PairStatus pairStatus)
        : range(range), distance(distance), stringCIGAR(CIGAR), strand(strand),
          pairStatus(pairStatus), indexBegin(range.getBegin()) {
        assert(CIGAR == "*"); // CIGAR must be *
    }

    /**
     * Constructor for an invalid text occurrence
     */
    TextOcc() : range(0, 0) {
    }

    /**
     * Creates an unmapped occurrence with SAM line for single-end alignment
     * @param bundle The read bundle with info about this sequence
     * @returns the unmapped occurrence.
     */
    static TextOcc createUnmappedSAMOccurrenceSE(const ReadBundle& bundle);

    /**
     * Creates an unmapped RHS record for single-end alignment. The hits column
     *contains an asterisk.
     * @param bundle The read bundle with info about this sequence
     * @returns the unmapped occurrence.
     */
    static TextOcc createUnmappedRHSOccurrenceSE(const ReadBundle& bundle) {
        TextOcc t; // default constructor
        t.outputLine = bundle.getSeqID() + "\t*";
        return t;
    }

    /**
     * Generates the SAM line for this occurrence.
     * @param seqID The sequence ID of the read.
     * @param printSeq The sequence of the read to be printed (* if not
     * the first match).
     * @param printQual The base quality to be printed (* if not the first
     * match).
     * @param nHits The number of other matches.
     * @param minScore The score of the best match.
     * @param primaryAlignment Indicates whether this is the primary alignment
     * @param seqNames the vector with all the reference sequence names
     */
    void generateSAMSingleEnd(const std::string& seqID,
                              const std::string& printSeq,
                              const std::string& printQual, length_t nHits,
                              length_t minScore, bool primaryAlignment,
                              const std::vector<std::string>& seqNames);

    /**
     * Generates the SAM line for this occurrence and put the other occurrences
     * in the XA tag.
     * @param seqID The sequence ID of the read.
     * @param printSeq The sequence of the read (or reverse complement) to be
     * printed
     * @param printQual The base quality to be printed
     * @param nHits The number of matches with the best score (includes this
     * occurrence!)
     * @param otherMatches The other matches.
     * @param seqNames the vector with all the reference sequence names
     */
    void generateSAMSingleEndXA(const std::string& seqID,
                                const std::string& printSeq,
                                const std::string& printQual, length_t nHits,
                                const std::vector<TextOcc>& otherMatches,
                                const std::vector<std::string>& seqNames);

    /**
     * Generates the SAM line for this occurrence. Given that this is the
     * first occurrence.
     * @param bundle the read bundle with info about this sequence
     * @param nHits The number of other matches.
     * @param minScore The score of the best match.
     * @param seqNames the vector with all the reference sequence names
     */
    void generateSAMSingleEndFirst(ReadBundle& bundle, length_t nHits,
                                   length_t minScore,
                                   const std::vector<std::string>& seqNames) {
#ifndef NO_REPORT

        std::string printSeq =
            (isRevCompl()) ? bundle.getRevComp() : bundle.getRead();
        std::string printQual =
            (isRevCompl()) ? bundle.getRevQuality() : bundle.getQual();

        generateSAMSingleEnd(bundle.getSeqID(), printSeq, printQual, nHits,
                             minScore, true, seqNames);
#endif
    }

    /**
     * Generates the SAM line for this occurrence, given that this is the
     * first occurrence. The other alignments are put in the XA tag.
     * @param bundle the read bundle with info about this sequence
     * @param nHits The number of  matches with an equal score (this included!)
     * @param minScore The score of the best match.
     * @param otherMatches The other matches.
     * @param seqNames
     */
    void generateSAMSingleEndXA(ReadBundle& bundle, length_t nHits,
                                length_t minScore,
                                const std::vector<TextOcc>& otherMatches,
                                const std::vector<std::string>& seqNames) {
#ifndef NO_REPORT
        std::string printSeq =
            (isRevCompl()) ? bundle.getRevComp() : bundle.getRead();
        std::string printQual =
            (isRevCompl()) ? bundle.getRevQuality() : bundle.getQual();

        if (printQual.empty()) {
            printQual = "*";
        }
        generateSAMSingleEndXA(bundle.getSeqID(), printSeq, printQual, nHits,
                               otherMatches, seqNames);
#endif
    }

    /**
     * Generates the SAM line for this occurrence. Given that this is not the
     * first occurrence.
     * @param seqID The sequence ID of the read.
     * @param nHits The number of other matches.
     * @param minScore The score of the best match.
     * @param seqNames the vector with all the reference sequence names
     */
    void
    generateSAMSingleEndNotFirst(const std::string& seqID, length_t nHits,
                                 length_t minScore,
                                 const std::vector<std::string>& seqNames) {
#ifndef NO_REPORT
        // print * for sequence and quality
        generateSAMSingleEnd(seqID, "*", "*", nHits, minScore, false, seqNames);
#endif
    }

    /**
     * Generates the SAM line for this occurrence, given that this is the
     * first occurrence. The other alignments are put in the XA tag.
     * @param bundle the read bundle with info about this sequence
     * @param otherMatches The other matches.
     * @param seqNames the vector with all the reference sequence names
     */
    void generateRHSSingleEnd(ReadBundle& bundle,
                              const std::vector<TextOcc>& otherMatches,
                              const std::vector<std::string>& seqNames) {
        outputLine = bundle.getSeqID();
        outputLine += "\t" + getRHS(seqNames);
        for (const auto& other : otherMatches) {
            outputLine += ";" + other.getRHS(seqNames);
        }
    }

    std::string getRHS(const std::vector<std::string>& seqNames) const {
        return "(" + seqNames[assignedSequenceID] + "," +
               fmt::to_string(distance) + ")";
    }

    /**
     * @returns the range of this occurrence
     */
    Range& getRange() {
        return range;
    }

    /**
     * @returns the range of this occurrence
     */
    const Range& getRange() const {
        return range;
    }

    /**
     * @returns the distance of this occurrence
     */
    const length_t getDistance() const {
        return distance;
    }

    /**
     * Sets the distance of this occurrence
     * @param distance the distance to set
     */
    void setDistance(length_t distance) {
        this->distance = distance;
    }

    void setPairStatus(PairStatus pairStatus) {
        this->pairStatus = pairStatus;
    }

    /**
     * @returns the output line of this occurrence
     */
    const std::string& getOutputLine() const {
        return outputLine;
    }

    const length_t getAssignedSequenceID() const {
        return assignedSequenceID;
    }

    const bool hasAssignedSequence() const {
        return seqNameChecked && seqNameFound != NOT_FOUND;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance , their length and finally on the existence of the CIGAR string
     * @param r the occurrence to compare to this (=right hand side of the
     * comparison)
     */
    bool operator<(const TextOcc& r) const {

        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        } else {
            // begin is equal, better ed is smarter
            if (distance != r.getDistance()) {
                return distance < r.getDistance();
            } else if (range.width() != r.getRange().width()) {
                // shorter read is smaller...
                return range.width() < r.getRange().width();
            } else {
                return false;
            }
        }
    }

    /**
     * Operator overloading.  An occurrence is smaller than a value if its
     * distance is smaller than the value.
     */
    bool operator<(length_t value) const {
        return getDistance() < value;
    }

    /**
     * Operator overloading.  Two TextOcc are equal if their ranges and distance
     * are equal.
     * @param r the occurrence to compare to this (=right hand side of the
     * comparison)
     */
    bool operator==(const TextOcc& r) {
        return r.getRange() == range && r.getDistance() == distance;
    }

    /**
     * @returns the width of this occurrence.
     */
    length_t width() const {
        return range.width();
    }

    /**
     * @returns true if this is a valid occurrence
     */
    bool isValid() const {
        return !range.empty();
    }

    /**
     * @returns true if the occurrence is along the reverse complemented strand,
     * false otherwise.
     */
    bool isRevCompl() const {
        return strand == REVERSE_C_STRAND;
    }

    bool isFirstReadInPair() const {
        return pairStatus == FIRST_IN_PAIR;
    }

    /**
     * Sets the CIGAR string of this occurrence.
     * @param cigar the CIGAR string to set
     */
    void setCigar(std::string& cigar) {
        // function does nothing
    }

    /**
     * @returns the CIGAR string of this occurrence.
     */
    std::string& getCigar() {
        assert(stringCIGAR == "*");
        return stringCIGAR;
    }

    /**
     * Sets the sequence to which this occurrence belongs.
     */
    void setAssignedSequence(SeqNameFound snFound, length_t seqID) {

        assignedSequenceID = seqID;
        seqNameChecked = true;
        seqNameFound = snFound;
    }

    void setAssignedSequenceNotFound() {

        assignedSequenceID = -1;
        seqNameChecked = true;
        seqNameFound = NOT_FOUND;
    }

    void removeTrimmingLabel() {
        assert(seqNameFound == FOUND_WITH_TRIMMING);
        seqNameFound = FOUND;
    }

    SeqNameFound getSeqNameFound() const {
        return seqNameFound;
    }

    const bool isSeqNameChecked() const {
        return seqNameChecked;
    }

    const length_t getEnd() const {
        return range.getEnd();
    }
    const length_t getBegin() const {
        return range.getBegin();
    }

    const PairStatus getPairStatus() const {
        return pairStatus;
    }
    const Strand getStrand() const {
        return strand;
    }

    void dropFields() {
        // we no longer need most fields so we can drop those we do not need
        stringCIGAR.clear();
    }

    void setIndexBegin(length_t indexBegin) {
        this->indexBegin = indexBegin;
    }

    const length_t getIndexBegin() const {
        return indexBegin;
    }
};

// ============================================================================
// CLASS SA RANGE PAIR
// ============================================================================

class ToeholdInterface {
  private:
    length_t toehold; // the toehold, which is one occurrence of the current
                      // match in the text
    bool toeholdRepresentsEnd; // indicates whether the toehold represents the
                               // end of the match
    length_t originalDepth;    // the depth of the original match, needed for
                            // offsetting the toehold at the end since the FMPos
                            // depth is decremented with shifts, and there is no
                            // way of knowing how muc

  public:
    /**
     * Constructor. Creates a toehold interface.
     * @param toehold the toehold
     * @param toeholdRepresentsEnd indicates whether the toehold represents the
     * @param originalDepth the depth of the original match (??)
     */
    ToeholdInterface(length_t toehold, bool toeholdRepresentsEnd,
                     length_t originalDepth)
        : toehold(toehold), toeholdRepresentsEnd(toeholdRepresentsEnd),
          originalDepth(originalDepth) {
    }
    /**
     * @returns the toehold
     */
    length_t getToehold() const {
        return toehold;
    }

    /**
     * @brief See if the toehold represents the end of the match
     *
     * @return true if the toehold represents the end of the match
     * @return false otherwise
     */
    bool getToeholdRepresentsEnd() const {
        return toeholdRepresentsEnd;
    }

    /**
     * @returns the original depth of the match
     */
    length_t getOriginalDepth() const {
        return originalDepth;
    }

    bool operator==(const ToeholdInterface& o) const {
        return o.getToehold() == toehold &&
               o.getToeholdRepresentsEnd() == toeholdRepresentsEnd &&
               o.getOriginalDepth() == originalDepth;
    }
};

/**
 * A pair of ranges. The first range is range over the suffix array of the
 * text. The second range is the corresponding range over the suffix array
 * of the reversed text
 */
class SARangePair
    : public ToeholdInterface // include toehold info in b-move
{
  private:
    SARange rangeSA;    // the range over the suffix array
    SARange rangeSARev; // the range over the suffix array of the reversed text

  public:
    /**
     * Default constructor, creates two empty ranges
     */
    SARangePair()
        :
          ToeholdInterface(0, false, 0),
          rangeSA(SARange()), rangeSARev(SARange()) {
    }

    /**
     * Constructor. Creates a pair of ranges.
     * @param rangeSA the range over the suffix array
     * @param rangeSARev the range over the suffix array of the reversed text
     * @param toehold the toehold
     * @param toeholdRepresentsEnd indicates whether the toehold represents the
     */
    SARangePair(SARange rangeSA, SARange rangeSARev, length_t toehold,
                bool toeholdRepresentsEnd, length_t originalDepth)
        : ToeholdInterface(toehold, toeholdRepresentsEnd, originalDepth),
          rangeSA(rangeSA), rangeSARev(rangeSARev) {
    }

    /**
     * @returns the range over the suffix array
     */
    const SARange& getRangeSA() const {
        return rangeSA;
    }

    /**
     * @returns the range over the suffix array of the reversed text
     */
    const SARange& getRangeSARev() const {
        return rangeSARev;
    }

    /**
     * Provides mutable access to the SA range.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the SA Range.
     */
    SARange& getRangeSAMutable() {
        return rangeSA;
    }

    /**
     * Provides mutable access to the reverse SA range.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the reverse SA Range.
     */
    SARange& getRangeSARevMutable() {
        return rangeSARev;
    }

    /**
     * @returns true if the range over the suffix array is empty, false
     * otherwise
     */
    bool empty() const {
        return rangeSA.empty();
    }

    /**
     * @returns the width of the ranges, calculated by using the range over the
     * suffix array
     */
    length_t width() const {
        return rangeSA.width();
    }
    /**
     * Operator overloading
     * @param o the other pair (=right hand side of equation)
     * @returns true if this is equal to rhs
     */
    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.getRangeSA() == rangeSA
               && ToeholdInterface::operator==(o)
            ;
    }

    /**
     * Function for debug purposes. With bidirectional search we expect this to
     * be true. With unidirectional backwards search this should be false.
     * @returns true if the ranges are synchronized, false otherwise
     */
    bool isSynchronized() const {
        return rangeSA.width() == rangeSARev.width();
    }
};

class SARangeWithToehold : public SARange, public ToeholdInterface {
  public:
    SARangeWithToehold() : SARange(), ToeholdInterface(0, false, 0) {
    }

    SARangeWithToehold(SARange range, length_t toehold,
                       bool toeholdRepresentsEnd, length_t originalDepth)
        : SARange(range),
          ToeholdInterface(toehold, toeholdRepresentsEnd, originalDepth) {
    }

    const SARange& getRange() const {
        return *this;
    }

    SARange getRange() {
        return *this;
    }
};

// ============================================================================
// SARANGEBACKWARDS DEFINITION
// ============================================================================

typedef SARangeWithToehold SARangeBackwards; // the range in the suffix array in
                                             // the move table, with toehold

// ============================================================================
// CLASS FMPos
// ============================================================================
/**
 * A position in the bidirectional FM-index.
 */
class FMPos {
  protected:
    SARangePair ranges; // the ranges over the suffix arrays
    length_t depth; // the depth of the prefix of the suffixes of this position

  public:
    /**
     * Default constructor for empty position (= empty ranges and depth of
     * zero)
     */
    FMPos() : ranges(SARangePair()), depth(0) {
    }

    FMPos(SARangePair& ranges, length_t depth) : ranges(ranges), depth(depth) {
    }

    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * Provides mutable access to the ranges.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the ranges.
     */
    SARangePair& getRangesMutable() {
        return ranges;
    }

    const length_t& getDepth() const {
        return depth;
    }

    void setRanges(SARangePair ranges) {
        this->ranges = ranges;
    }

    void setDepth(length_t depth) {
        this->depth = depth;
    }

    /**
     * Operator overloading, two FMPos are equal if their ranges and depth
     * are equal
     * @param rhs the FMPos to compare to this
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMPos& rhs) const {
        return ranges == rhs.getRanges() && depth == rhs.getDepth();
    }
    /**
     * @returns true if the ranges are not empty, false otherwise
     */
    bool isValid() const {
        return !ranges.empty();
    }
};
// ============================================================================
// CLASS FMOcc
// ============================================================================

/**
 * An occurrence in the bidirectional FM-index
 */
class FMOcc {
  private:
    FMPos pos;         // The FM position of this occurrence
    length_t distance; // the distance (hamming or edit)
    length_t shift; // A right-sift to the corresponding positions in the text

    Strand strand = FORWARD_STRAND; // indicates whether this occurrence is
                                    // reverse complemented
    PairStatus pairStatus = FIRST_IN_PAIR; // indicates whether this occurrence
                                           // is the first read in a pair

  public:
    FMOcc() : pos(), distance(0), shift(0) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param ranges the ranges of this approximate match (range in SA and
     * in SA')
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param depth the depth (=length) of this approximate match
     * @param strand the strand on which this occurrence lies, defaults to
     * FORWARD_STRAND
     * @param pairStatus indicates whether this occurrence is the first or
     * second in the pair, defaults to FIRST_IN_PAIR
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(SARangePair ranges, length_t distance, length_t depth,
          Strand strand = FORWARD_STRAND, PairStatus pairStatus = FIRST_IN_PAIR,
          length_t shift = 0)
        : pos(ranges, depth), distance(distance), shift(shift), strand(strand),
          pairStatus(pairStatus) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param pos the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param strand the strand on which this occurrence lies
     * @param pairStatus indicates whether this occurrence is the first or
     * second in the pair
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(FMPos pos, length_t distance, Strand strand, PairStatus pairStatus,
          length_t shift = 0)
        : pos(pos), distance(distance), shift(shift), strand(strand),
          pairStatus(pairStatus) {
    }

    const SARangePair& getRanges() const {
        return pos.getRanges();
    }
    const length_t& getDistance() const {
        return distance;
    }

    const length_t& getDepth() const {
        return pos.getDepth();
    }

    length_t getWidth() const {
        return pos.getRanges().width();
    }

    const length_t& getShift() const {
        return shift;
    }

    void setRanges(SARangePair ranges) {
        pos.setRanges(ranges);
    }

    void setDistance(length_t distance) {
        this->distance = distance;
    }

    void setDepth(length_t depth) {
        pos.setDepth(depth);
    }

    void setStrand(Strand strand) {
        this->strand = strand;
    }

    void setPairStatus(PairStatus pairStatus) {
        this->pairStatus = pairStatus;
    }

    /**
     * @returns true if the position is valid, false otherwise
     */
    bool isValid() const {
        return pos.isValid();
    }

    bool isRevCompl() const {
        return strand == REVERSE_C_STRAND;
    }

    bool isFirstReadInPair() const {
        return pairStatus == FIRST_IN_PAIR;
    }

    Strand getStrand() const {
        return strand;
    }
    PairStatus getPairStatus() const {
        return pairStatus;
    }
    /**
     * Operator overloading to sort FMOcc
     * First the FMOcc are sorted on the begin of the range over the suffix
     * array of their position Then they are sorted on their distance
     * Lastly they are sorted on the width of their ranges
     * @param rhs the FMOcc to compare to this
     * @returns true if this is smaller than rhs
     */
    bool operator<(const FMOcc& rhs) const {
        if (pos.getRanges().getRangeSA().getBegin() !=
            rhs.getRanges().getRangeSA().getBegin()) {
            // begins of ranges are unequal, return smallest begin
            return getRanges().getRangeSA().getBegin() <
                   rhs.getRanges().getRangeSA().getBegin();
        } else if (distance != rhs.getDistance()) {
            // begin is equal, better ed is smarter
            return distance < rhs.getDistance();
        } else if (getRanges().width() != rhs.getRanges().width()) {
            // shorter read is smaller...
            return getRanges().width() < rhs.getRanges().width();
        } else {
            // prefer no shift
            return getShift() < rhs.getShift();
        }
    }
    /**
     * Operator overloading
     * Two FMOcc are equal if their ranges, distance, depth and shift are all
     * equal
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMOcc& rhs) {
        return getRanges() == rhs.getRanges() &&
               distance == rhs.getDistance() && getDepth() == rhs.getDepth() &&
               getShift() == rhs.getShift() && isRevCompl() == rhs.isRevCompl();
    }
    friend std::ostream& operator<<(std::ostream& os, const FMOcc& fmOcc);
};

// ============================================================================
// CLASS FMPosExt
// ============================================================================
/**
 * A single node in the bidirectional FM-index. Its depth is the depth from
 * the start match for a particular phase of a search
 */
class FMPosExt : public FMPos {
  private:
    char c;                // the character of this node
    bool reported = false; // has this particular node already reported?
  public:
    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param ranges the ranges over the suffix and reversed suffix array
     * that go to this node
     * @param row the row of this node in the alignment matrix = depth of
     * this node
     */
    FMPosExt(char character, SARangePair ranges, length_t row)
        : FMPos(ranges, row), c(character), reported(false) {
    }

    /**
     * Default constructor, this Node will have empty ranges
     */
    FMPosExt() : FMPos(), c(char(0)) {
    }

    /**
     * Sets the report flag to true
     */
    void report() {
        reported = true;
    }

    /**
     * Reports the match (with added depth) at this node,
     * @param occ the match will be stored here
     * @param startDepth the depth to add to the match
     * @param EDFound the found edit distance for this node
     * @param noDoubleReports false if this node is allowed to report more
     * than once, defaults to false
     * @param shift right shift of the match, defaults to zero
     */
    void report(FMOcc& occ, const length_t& startDepth, const length_t& EDFound,
                const bool& noDoubleReports = false, length_t shift = 0) {
        if (!reported) {
            Strand strand =
                FORWARD_STRAND; // strand does not matter here-> just choose
            PairStatus pairStatus =
                FIRST_IN_PAIR; // pairStatus does not matter here
            occ = FMOcc(getRanges(), EDFound, depth + startDepth, strand,
                        pairStatus, shift);

            // if finalPiece, report only once
            if (noDoubleReports) {
                report();
            }
        }
    }

    /**
     * Gets the ranges of this node
     * @returns the ranges of this node
     */
    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * Get the character of this node
     * @returns the character of this node
     */
    const char getCharacter() const {
        return c;
    }

    /**
     * Get the row of this node
     * @returns the row of this node
     */
    length_t getRow() const {
        return depth;
    }
};

// ============================================================================
// CLASS CLUSTER
// ============================================================================

/**
 * Cluster class to represent the end of an alignment phase in the bidirectional
 * index using a search scheme. The cluster contains the final column of the
 * alignment matrix with the corresponding distance scores and nodes in the
 * index. The cluster can report the centers of the cluster and the deepest
 * local minimum for efficient switching to a the next phase in the search
 * scheme.
 */
class Cluster {
  private:
    std::vector<uint16_t> eds;   // the edit distances of this cluster
    std::vector<FMPosExt> nodes; // the nodes of this cluster

    length_t lastCell;   // the lastCell of the cluster that was filled in
    uint16_t maxED;      // the maxEd for this cluster
    length_t startDepth; // the startDepth for this cluster (= depth of match
                         // before matrix of this cluster)

    length_t shift; // the right shift of the occurrences in the text
  public:
    /**
     * Constructor
     * @param size the size of the cluster
     * @param maxED the maximal allowed edit distance
     * @param startDepth the depth before this cluster
     * @param shift the right shift of the occurrences in the text
     */
    Cluster(length_t size, length_t maxED, length_t startDepth, length_t shift)
        : eds(size, maxED + 1), nodes(size), lastCell(-1), maxED(maxED),
          startDepth(startDepth), shift(shift) {
    }

    /**
     * Sets the ed and node at index idx to ed and node. Also updates
     * lastCell to be idx
     * @param idx the idx to change
     * @param node the node to set at index idx
     * @param ed the ed to set at index idx
     */
    void setValue(length_t idx, const FMPosExt& node, const length_t& ed) {
        eds[idx] = ed;
        nodes[idx] = node;
        lastCell = idx;
    }

    /**
     * Returns the size of this cluster
     */
    const length_t size() const {
        return eds.size();
    }

    /**
     * @returns vector with all nodes in the cluster that are a centre and
     * under the maximal allowed distance, be aware that if there are
     * multiple centers in the cluster it is very likely that only one of
     * them will be redundant, but the others might eliminate another
     * occurrence
     */
    std::vector<FMOcc> reportCentersAtEnd() {

        std::vector<FMOcc> centers;
        centers.reserve(lastCell + 1);

        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] <= maxED && (i == 0 || eds[i] <= eds[i - 1]) &&
                (i == lastCell || eds[i] <= eds[i + 1])) {
                FMOcc m;
                nodes[i].report(m, startDepth, eds[i], true, shift);
                centers.emplace_back(m);
            }
        }

        return centers;
    }

    /**
     * @returns approximate match that corresponds to the ranges of the
     * deepest global minimum of this cluster, but with the depth of the
     * highest global minimum of this cluster. If the direction is backward
     * a shift will be set such that the occurrence in the text will be as
     * short as possible
     */
    FMOcc reportDeepestMinimum(Direction dir) {
        uint16_t minED = maxED + 1;
        length_t highestBestIdx = -1;
        length_t deepestBestIdx = -1;

        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] < minED) {
                minED = eds[i];
                highestBestIdx = i;
                deepestBestIdx = i;
            }
            if (eds[i] == minED) {
                deepestBestIdx = i;
            }
        }
        FMOcc m;
        if (minED <= maxED) {
            nodes[deepestBestIdx].report(
                m, startDepth - (deepestBestIdx - highestBestIdx), minED, true,
                ((dir == BACKWARD) ? (deepestBestIdx - highestBestIdx) : 0) +
                    shift);
        }
        return m;
    }

    /**
     * This method returns a match that corresponds to the highest cluster
     * centre. Its descendants and the corresponding initialization edit
     * distances are updated. Eds of descendants that are part of a cluster
     * centre which is lower than the lower bound will be updated in the initEds
     * vector
     * @param lowerBound the lower bound for this iteration
     * @param desc the descendants of the highest cluster centre, these
     * will be inserted during the method
     * @param initEds the initialization eds for the next iteration, these
     * correspond to the eds of the highest centre and its descendants,
     * where eds part of a cluster of which the centre is below the
     * lower bound are updated. These values will be inserted during the
     * method
     * @returns The occurrence corresponding to the upper cluster centre
     * which has a valid distance
     */
    FMOcc getClusterCentra(uint16_t lowerBound, std::vector<FMPosExt>& desc,
                           std::vector<uint16_t>& initEds);
};

// ============================================================================
// CLASS COUNTERS
// ============================================================================
/**
 * @brief Class to manage and manipulate various performance counters.
 */
class Counters {
  public:
    enum CounterType {
        NODE_COUNTER, // counts the number of nodes visited in the index
        TOTAL_REPORTED_POSITIONS, // counts the number of matches reported
                                  // (either via in-text verification or
                                  // in-index matching)
        SEARCH_STARTED,   // counts the number of times a search started
        DROPPED_UNIQUE_MATCHES, // counts the number of unique  matches dropped
                                // because of one-line reporting

        // Aggregational counters
        NUMBER_OF_READS,           // counts the number of reads
        TOTAL_UNIQUE_MATCHES,      // counts the number of all unique matches
        MAPPED_READS,              // counts the number of mapped reads
        MAPPED_PAIRS,              // number of mapped pairs
        MAPPED_HALF_PAIRS, // number of pairs for which only one read mapped

        COUNTER_TYPE_MAX // To denote the number of counters
    };

  private:
    std::array<uint64_t, COUNTER_TYPE_MAX> counters; // array with the counters
  public:
    /**
     * @brief Constructor to initialize counters. Resets all counters to 0.
     */
    Counters() {
        resetCounters();
    };

    /**
     * Reset all counters to 0
     */
    void resetCounters() {
        counters.fill(0);
    }

    /**
     * @brief Increments the specified counter by a given amount.
     *
     * @param type The type of counter to increment.
     * @param amount [optional] The amount by which to increment the counter
     * (default is 1).
     */
    void inc(CounterType type, uint64_t amount = 1) {
        counters[type] += amount;
    }

    /**
     * @brief Retrieves the current value of the specified counter.
     *
     * @param type The type of counter to retrieve.
     * @return The current value of the counter.
     */
    uint64_t get(CounterType type) const {
        return counters[type];
    }

    /**
     * @brief Adds counters from another Counters object to this one.
     *
     * @param o The Counters object from which to add counters.
     */
    void addCounters(const Counters& o) {
        for (int i = 0; i < COUNTER_TYPE_MAX; ++i) {
            counters[i] += o.counters[i];
        }
    }

    /**
     * @brief Reports statistics based on the counters and a given sequencing
     * mode.
     *
     * @param sMode The sequencing mode to use for reporting statistics.
     */
    void reportStatistics(const SequencingMode& sMode) const;
};

// ============================================================================
// CLASS Occurrences
// ============================================================================
// This class combines the occurrences in the text and the occurrences
// in the fm index into one data structure
/**
 * @class Occurrences
 * @brief Represents a collection of occurrences in the FM index and the text.
 *
 * The Occurrences class stores the in-index occurrences (FMOcc) and the in-text
 * occurrences (TextOcc). It provides methods to add occurrences, erase
 * duplicates, and retrieve the occurrences.
 */
class Occurrences {
  private:
    std::vector<TextOcc> inTextOcc; // the in-text occurrences
    std::vector<FMOcc> inFMOcc;     // the in-index occurrences

  public:
    /**
     * @brief Constructs an Occurrences object with optional reserve capacity.
     *
     * @param reserve The initial capacity to reserve for the occurrences
     * vectors.
     */
    Occurrences(const length_t reserve = 200) {
        inTextOcc.reserve(reserve);
        inFMOcc.reserve(reserve);
    }

    /**
     * @brief Adds an in-index occurrence to the collection.
     *
     * @param match The in-index occurrence to add.
     */
    void addFMOcc(const FMOcc& match) {
        inFMOcc.emplace_back(match);
    }

    /**
     * @brief Adds an in-index occurrence to the collection.
     *
     * @param ranges The SA ranges of the occurrence.
     * @param score The score of the occurrence.
     * @param depth The depth of the occurrence.
     * @param strand The strand of the occurrence.
     * @param pairStatus The pair status of the read.
     */
    void addFMOcc(const SARangePair& ranges, const length_t& score,
                  const length_t& depth, const Strand strand,
                  const PairStatus pairStatus) {
        inFMOcc.emplace_back(ranges, score, depth, strand, pairStatus);
    }

    /**
     * @brief Adds an in-index occurrence to the collection.
     *
     * @param currentNode The current FM position and extension of the
     * occurrence.
     * @param score The score of the occurrence.
     * @param strand The strand of the occurrence.
     * @param pairStatus The pair status of the read.
     */
    void addFMOcc(const FMPosExt& currentNode, const length_t& score,
                  const Strand strand, const PairStatus pairStatus) {
        inFMOcc.emplace_back(currentNode, score, strand, pairStatus);
    }

    /**
     * @brief Adds an in-text occurrence to the collection.
     *
     * @param occ The in-text occurrence to add.
     */
    void addTextOcc(const TextOcc& occ) {
        inTextOcc.emplace_back(occ);
    }

    /**
     * @brief Adds an in-text occurrence to the collection.
     *
     * @param range The range of the occurrence in the text.
     * @param score The score of the occurrence.
     * @param CIGAR The CIGAR string of the occurrence.
     * @param strand The strand of the occurrence.
     * @param pairStatus The pair status of the read.
     */
    void addTextOcc(const Range& range, const length_t& score,
                    const std::string& CIGAR, Strand strand,
                    PairStatus pairStatus) {
        inTextOcc.emplace_back(range, score, CIGAR, strand, pairStatus);
    }

    /**
     * @brief Adds an in-text occurrence to the collection.
     *
     * @param range The range of the occurrence in the text.
     * @param score The score of the occurrence.
     * @param strand The strand of the occurrence.
     * @param pairStatus The pair status of the read.
     */
    void addTextOcc(const Range& range, const length_t score, Strand strand,
                    PairStatus pairStatus) {
        inTextOcc.emplace_back(range, score, strand, pairStatus);
    }

    /**
     * @brief Adds the given occurrences to the inTextOcc vector.
     *
     * @param occs The occurrences to add.
     */
    void addTextOccs(const std::vector<TextOcc>& occs) {
        inTextOcc.insert(inTextOcc.end(), occs.begin(), occs.end());
    }

    /**
     * @brief Sets the inTextOcc vector to the given occurrences.
     *
     * @param occs The occurrences to set.
     */
    void setTextOccs(std::vector<TextOcc>& occs) {
        inTextOcc = occs;
    }

    /**
     * @brief Sets the inTextOcc vector to the given occurrences.
     * The given occurrences are moved into the inTextOcc vector.
     * @param occs The occurrences to set.
     */
    void moveTextOccs(std::vector<TextOcc>& occs) {
        inTextOcc = std::move(occs);
    }

    /**
     * @brief Erases all duplicate in-index occurrences and sorts the
     * occurrences.
     */
    void eraseDoublesFM() {
#ifdef DEVELOPER_MODE
        stable_sort(inFMOcc.begin(), inFMOcc.end());
#else
        sort(inFMOcc.begin(), inFMOcc.end());
#endif

        inFMOcc.erase(unique(inFMOcc.begin(), inFMOcc.end()), inFMOcc.end());
    }

    /**
     * @brief Erases all duplicate in-text occurrences and sorts the
     * occurrences.
     */
    void eraseDoublesText() {
#ifdef DEVELOPER_MODE
        stable_sort(inTextOcc.begin(), inTextOcc.end());
#else
        sort(inTextOcc.begin(), inTextOcc.end());
#endif
        inTextOcc.erase(unique(inTextOcc.begin(), inTextOcc.end()),
                        inTextOcc.end());
    }

    /**
     * @brief Returns the in-index occurrences.
     *
     * @return The in-index occurrences.
     */
    const std::vector<FMOcc>& getFMOccurrences() const {
        return inFMOcc;
    };

    /**
     * @brief Returns the number of in-text occurrences.
     *
     * @return The number of in-text occurrences.
     */
    size_t textOccSize() const {
        return inTextOcc.size();
    }

    /**
     * @brief Returns the in-text occurrences.
     *
     * @return The in-text occurrences.
     */
    const std::vector<TextOcc>& getTextOccurrences() const {
        return inTextOcc;
    }

    /**
     * @brief Returns the maximum size between in-text and in-index occurrences.
     *
     * @return The maximum size between in-text and in-index occurrences.
     */
    length_t getMaxSize() const {
        return std::max(inTextOcc.size(), inFMOcc.size());
    }

    /**
     * Checks if the Occurrences object is empty.
     *
     * @return true if the Occurrences object is empty, false otherwise.
     */
    bool empty() const {
        return inTextOcc.empty() && inFMOcc.empty();
    }
};

#endif // end of include guard: INDEXHELPERS_H