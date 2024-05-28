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
#ifndef RINDEXHELPERS_H
#define RINDEXHELPERS_H

#include "indexhelpers.h"

#include <sdsl/construct.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/lcp_bitcompressed.hpp>

#include <iostream>
#include <cstdint>
#include <type_traits>


using namespace std;


// ===========================================================================
// CLASS Position
// ===========================================================================

class Position {

    public:

        length_t runIndex;
        length_t posIndex;

        Position(): runIndex(0), posIndex(0){}

        Position(length_t posIndex, length_t runIndex): runIndex(runIndex), posIndex(posIndex) {}

        bool operator==(const Position& p2) const {
            return posIndex == p2.posIndex && runIndex == p2.runIndex;
        }

        static pair<Position, Position> getEmptyRange() {
            pair<Position, Position> result;
            result.first = Position(1, 0);
            result.second = Position(0, 0);

            return result;
        }

        static Position getPositionWithIndex(length_t posIndex, const Position& prevPos) {
            return Position(posIndex, prevPos.runIndex);
        }

};


// ============================================================================
// CLASS PositionRange
// ============================================================================

class PositionRange {
    private:
        Position beginPos; // beginning run index of the range
        Position endPos;   // end run index of the range (inclusive)
        bool runIndicesValid;

    public:
        PositionRange() {
            beginPos = Position(1, 0);
            endPos = Position(0, 0);
            runIndicesValid = false;
        }

        PositionRange(Position begin, Position end): 
        beginPos(begin), endPos(end), runIndicesValid(true) {}

        static PositionRange getEmptyPositionRange() {
            pair<Position, Position> emptyRange = Position::getEmptyRange();

            return PositionRange(emptyRange.first, emptyRange.second);
        }

        const Position& getBeginPos() const {
            return beginPos;
        }

        const Position& getEndPos() const {
            return endPos;
        }

        Position& getBeginPosMutable() {
            return beginPos;
        }

        Position& getEndPosMutable() {
            return endPos;
        }

        /**
         * Check if this range is empty
         * @returns true if the range is empty, false otherwise
         */
        bool empty() const {
            return endPos.posIndex < beginPos.posIndex;
        }

        /**
         * Gets the width of the range (end - begin)
         * @returns the width of this range
         */
        length_t width() const {
            return (empty()) ? 0 : endPos.posIndex - beginPos.posIndex + 1;
        }

        bool getRunIndicesValid() {
            return runIndicesValid;
        }

        void setRunIndicesValid(bool valid) {
            this->runIndicesValid = valid;
        }

        /**
        * Operator overloading, two ranges are equal if their begin and end field
        * are equal
        */
        bool operator==(const PositionRange& o) const {
            return beginPos == o.getBeginPos() && endPos == o.getEndPos();
        }

};


// ============================================================================
// CLASS SAPositionRangePair
// ============================================================================

class SAPositionRangePair{
    private:
        PositionRange rangeSA;
        PositionRange rangeSARev;

    public:
        SAPositionRangePair() {}

        SAPositionRangePair(PositionRange rangeSA, PositionRange rangeSARev): 
        rangeSA(rangeSA), rangeSARev(rangeSARev) {}

        const PositionRange& getRangeSA() const {
            return rangeSA;
        }

        const PositionRange& getRangeSARev() const {
            return rangeSARev;
        }
        /**
        * @returns true if the ranges are empty, false otherwise
        */
        bool empty() const {
            return rangeSA.empty();
        }

        length_t width() const {
            return rangeSA.width();
        }
        /**
        * Operator overloading
        * @returns true if this is equal to rhs
        */
        bool operator==(const SAPositionRangePair& o) const {
            // only the first range matters as the ranges imply each other
            return o.getRangeSA() == rangeSA && o.getRangeSARev() == rangeSARev;
        }

};



// ============================================================================
// CLASS BRSample
// ============================================================================

/**
 * Bidirectional b-move sample
 * Holds a SAPositionRangePair with a range over de SA and reverse SA,
 * and the j, d and len variables, as explained in the br-index paper
*/

class BRSample {
    private:
        SAPositionRangePair ranges;
        // let p be the position of the current sample in the SA 
        length_t textPos;   // position index of the toehold
        length_t offset;    // (= d), offset of textPos to the starting position of the current pattern
        length_t length;    // (= len), length of the current pattern

    public:

        /**
         * Default constructor for empty position 
         * (= empty ranges and all variables are zero)
         */
        BRSample() : ranges(SAPositionRangePair()), textPos(0), offset(0), length(0) {}

        BRSample(SAPositionRangePair & ranges, length_t textPos, length_t offset, length_t length) 
            : ranges(ranges), textPos(textPos), offset(offset), length(length) {}

        const SAPositionRangePair & getRanges() const {
            return ranges;
        }
        
        const length_t& getTextPos() const {
            return textPos;
        }

        const length_t& getOffset() const {
            return offset;
        }

        const length_t& getLength() const {
            return length;
        }

        void setRanges(SAPositionRangePair  ranges) {
            this->ranges = ranges;
        }

        void setTextPos(length_t textPos) {
            this->textPos = textPos;
        }

        void setOffset(length_t offset) {
            this->offset = offset;
        }

        void setLength(length_t length) {
            this->length = length;
        }

        /**
         * @return width of the ranges of the sample
        */
        length_t width() const {
            return ranges.width();
        }

        /**
         * @return true if the ranges are not empty, false otherwise
         */
        bool isValid() const {
            return !ranges.empty();
        }

        /** Operator overloading, two BRsamples are equal if all their attirbutes are equal
         * @param rhs the BRsample to compare to this
         * @return true if this is equal to rhs
         */
        bool operator==(const BRSample& rhs) const {
            return ranges == rhs.getRanges() && textPos == rhs.getTextPos()
                && offset == rhs.getOffset() && length == rhs.getLength();
        }
};

// ============================================================================
// CLASS BRPos
// ============================================================================

/**
 * A position in b-move.
 */

class BRPos {
    protected:
        BRSample sample; // the b-move sample
        length_t depth; // the depth of the prefix of the suffixes of this position

    public:
        /**
         * Default constructor for empty position (= empty sample and depth of
         * zero)
         */
        BRPos() : sample(BRSample()), depth(0) {}

        BRPos(BRSample sample, length_t depth) : sample(sample), depth(depth) {}

        const BRSample& getSample() const {
            return sample;
        }

        const length_t& getDepth() const {
            return depth;
        }

        void setSample(BRSample sample) {
            this->sample = sample;
        }

        void setDepth(length_t depth) {
            this->depth = depth;
        }

        /**
         * Operator overloading, two BRPos are equal if their ranges and depth
         * are equal
         * @param rhs the BRPos to compare to this
         * @returns true if this is equal to rhs
         */
        bool operator==(const BRPos& rhs) const {
            return sample == rhs.getSample() && depth == rhs.getDepth();
        }

        /**
         * @returns true if the sample is valid, false otherwise
         */
        bool isValid() const {
            return sample.isValid();
        }
};

// ============================================================================
// CLASS BROcc
// ============================================================================

/**
 * An occurrence in b-move
 */

class BROcc {
    private:
        BRPos pos;          // The position of this occurrence
        length_t distance;  // the distance (hamming or edit)
        length_t shift;     // A right-sift to the corresponding positions in the text

    public:
        BROcc() : pos(), distance(0), shift(0) {}

        /**
         * Make a bidirectional approximate match in the suffix array
         * @param sample the BRSample of this approximate match 
         * @param distance the (edit or hamming) distance of this approximate
         * match
         * @param depth the depth (=length) of this approximate match
         * @param shift The right shift to the corresponding positions in the
         * text, defaults to zero
         */
        BROcc(BRSample sample, length_t distance, length_t depth,
            length_t shift = 0)
            : pos(sample, depth), distance(distance), shift(shift) {}
        
        /**
         * Make a bidirectional approximate match in the suffix array
         * @param pos the position in the RIndex of this approximate match
         * @param distance the (edit or hamming) distance of this approximate
         * match
         * @param shift The right shift to the corresponding positions in the
         * text, defaults to zero
         */
        BROcc(BRPos pos, length_t distance, length_t shift = 0)
            : pos(pos), distance(distance), shift(shift) {}

        const BRSample& getSample() const {
            return pos.getSample();
        }

        const SAPositionRangePair& getRanges() const {
            return pos.getSample().getRanges();
        }

        const length_t& getDistance() const {
            return distance;
        }

        const length_t& getDepth() const {
            return pos.getDepth();
        }

        length_t getWidth() const {
            return pos.getSample().getRanges().width();
        }

        const length_t& getShift() const {
            return shift;
        }

        void setSample(BRSample sample) {
            pos.setSample(sample);
        }

        void setDistance(length_t distance) {
            this->distance = distance;
        }

        void setDepth(length_t depth) {
            pos.setDepth(depth);
        }
        /**
         * @returns true if the position is valid, false otherwise
         */
        bool isValid() const {
            return pos.isValid();
        }

        /**
         * Operator overloading to sort BROcc
         * First the BROcc are sorted on the begin of the range over the suffix
         * array of their position Then they are sorted on their distance
         * Lastly they are sorted on the width of their ranges
         * @param rhs the BROcc to compare to this
         * @returns true if this is smaller than rhs
         */
        bool operator<(const BROcc& rhs) const {
            if (getRanges().getRangeSA().getBeginPos().posIndex !=
                rhs.getRanges().getRangeSA().getBeginPos().posIndex) {
                // begins of ranges are unequal, return smallest begin
                return getRanges().getRangeSA().getBeginPos().posIndex <
                    rhs.getRanges().getRangeSA().getBeginPos().posIndex;
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
         * Two BROcc are equal if their samples, distance, depth and shift are all
         * equal
         * @param returns true if this is equal to rhs
         */
        bool operator==(const BROcc& rhs) {
            return getSample() == rhs.getSample() && distance == rhs.getDistance() 
                && getDepth() == rhs.getDepth() && getShift() == rhs.getShift();
        }

};

// ============================================================================
// CLASS BRPosExt
// ============================================================================

/**
 * A single node in b-move. Its depth is the depth from
 * the startmatch for a particular phase of a search
 */

class BRPosExt : public BRPos {
    private:
        char c;                // the character of this node
        bool reported = false; // has this particular node already reported?

    public:
        /**
         * Create a node of the search tree
         * @param character the character of this node
         * @param sample the sample of this node
         * @param row the row of this node in thealignment matrix = depth of
         * this node
         */
        BRPosExt(char character, BRSample sample, length_t row)
            : BRPos(sample, row), c(character), reported(false) {}

        /**
         * Default constructor, this Node will have empty ranges
         */
        BRPosExt() : BRPos(), c(char(0)) {}

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
         * @param shift, right shift of the matche, defaults to zero
         */
        void report(BROcc& occ, const length_t& startDepth, const length_t& EDFound,
                    const bool& noDoubleReports = false, length_t shift = 0) {
            if (!reported) {
                occ = BROcc(this->getSample(), EDFound, this->depth + startDepth, shift);

                // if finalPiece, report only once
                if (noDoubleReports) {
                    report();
                }
            }
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
            return this->depth;
        }
};


// ============================================================================
// type definitions
// ============================================================================

typedef Cluster<BROcc, BRPosExt> BRCluster;
typedef Occurrences<BROcc, BRPosExt> BROccurrences;

#endif // RINDEXHELPERS_H