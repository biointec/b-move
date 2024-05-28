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

#include "bitparallelmatrix.h"
#include "substring.h"
#include <sstream>
#include <cstdint>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

#include "wordlength.h"

// ============================================================================
// HELPERS
// ============================================================================

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

#define MAX_MAPQ 60

// ============================================================================
// CLASS RANGE
// ============================================================================

class Range {
  private:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * Constructor
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    Range() : begin(0), end(0) {
    }

    length_t getBegin() const {
        return begin;
    }
    length_t getEnd() const {
        return end;
    }
    /**
     * Check if this range is empty
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
     * are equal
     */
    bool operator==(const Range& o) const {
        return o.getBegin() == begin && o.getEnd() == end;
    }

    /**
     * Operator overloading. Outputs the range as [begin, end) to the output stream
     * @param os the output stream
     * @param r the range to print
     */
    friend std::ostream& operator<<(std::ostream& os, const Range& r) {
        os << "[" << r.begin << ", " << r.end << ")";
        return os;
    }
};

// ============================================================================
// CLASS TextOccurrence
// ============================================================================

class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)

    std::string samLine; // the corresponding output for this occurrence (for
    // now a custom format)

    uint16_t getFlagsSE(bool reverseComplement, bool notPrimary) const {
        uint16_t result = 0;

        // set rev complement flag
        result |= (reverseComplement << 4);
        // set not primary flag
        result |= (notPrimary << 8);

        return result;
    }

    length_t getMapQ(length_t score, length_t hitNumbers,
                     length_t minScore) const {
        if (score != minScore) {
            return 0;
        }
        if (hitNumbers == 1) {
            return MAX_MAPQ;
        } else {
            return round(-10.0 * log10(1 - 1.0 / hitNumbers));
        }
    }

  public:
    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     */
    TextOcc(Range range, length_t distance)
        : range(range), distance(distance), samLine() {
    }

    /**
     * Constructor for an invalid text occurrence
     */
    TextOcc() : range(0, 0) {
    }

    void generateSamSE(bool revCompl, bool notPrimary, const std::string& seqID,
                       const std::string& seq, const std::string& qual,
                       length_t nHits, length_t minscore) {
        std::stringstream s;

        s << seqID << "\t";                            // read name
        s << getFlagsSE(revCompl, notPrimary) << "\t"; // flags
        s << "*\t";                        // reference sequence name
        s << range.getBegin() + 1 << "\t"; // 1-based pos in ref seq
        s << getMapQ(distance, nHits, minscore) << "\t"; // mapping quality
        s << "*\t";       // CIGAR string not yet supported
        s << "*\t";       // mate ref seq name, always set to *
        s << "0\t";       // mate pos, always set to 0
        s << "0\t";       // inferred insert size, always set to 0
        s << seq << "\t"; // read sequence
        s << (qual.empty() ? "*" : qual) << "\t"; // read quality
        s << "AS:i:" << distance << "\t";         // alignment score
        s << "PG:Z:b-move";                      // program name
        samLine = s.str();
    }

    const Range getRange() const {
        return range;
    }
    const length_t getDistance() const {
        return distance;
    }

    const std::string& getSamLine() const {
        return samLine;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance and finally their length
     */
    bool operator<(const TextOcc& r) {

        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        } else {
            // begin is equal, better ed is smarter
            if (distance != r.getDistance()) {
                return distance < r.getDistance();
            } else {
                // shorter read is smaller...
                return range.width() < r.getRange().width();
            }
        }
    }

    bool operator==(const TextOcc& r) {
        return r.getRange() == range && r.getDistance() == distance;
    }

    bool isValid() const {
        return !range.empty();
    }

    length_t width() const {
        return range.width();
    }
};

// ============================================================================
// CLASS CLUSTER
// ============================================================================

template<class IndexOcc, class IndexPosExt>
class Cluster {
  private:
    std::vector<uint> eds;       // the edit distances of this cluster
    std::vector<IndexPosExt> nodes; // the nodes of this cluster

    length_t lastCell;   // the lastCell of the cluster that was filled in
    uint maxED;          // the maxEd for this cluster
    length_t startDepth; // the startdepth for this cluster (= depth of match
                         // before matrix of this cluster)

    length_t shift; // the right shift of the occurrences in the text
  public:
    /**
     * Constructor
     * @param size, the size of the cluster
     * @param maxED, the maximal allowed edit distance
     * @param startDepth, the depth before this cluster
     * @param shift, the right shift of the occurrences in the text
     */
    Cluster(length_t size, length_t maxED, length_t startDepth, length_t shift)
        : eds(size, maxED + 1), nodes(size), lastCell(-1), maxED(maxED),
          startDepth(startDepth), shift(shift) {
    }

    /**
     * Sets the ed and node at index idx to ed and node. Also updates
     * lastCell to be idx
     * @param idx, the idx to change
     * @param node, the node to set at index idx
     * @param ed, the ed to set at index idx
     */
    void setValue(length_t idx, const IndexPosExt& node, const length_t& ed) {
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
    std::vector<IndexOcc> reportCentersAtEnd() {

        std::vector<IndexOcc> centers;
        centers.reserve(lastCell + 1);

        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] <= maxED && (i == 0 || eds[i] <= eds[i - 1]) &&
                (i == lastCell || eds[i] <= eds[i + 1])) {
                IndexOcc m;
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
    IndexOcc reportDeepestMinimum(Direction dir) {
        uint minED = maxED + 1;
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
        IndexOcc m;
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
     * centre. Its descendants and the corresponding initializationeds are
     * updated. Eds of descendants that are part of a cluster centre which
     * is lower than the lower bound will be updated in the initEds vector
     * @param lowerBound, the lower bound for this iteration
     * @param desc, the descendants of the highest cluster centre, these
     * will be inserted during the method
     * @param initEds, the initialization eds for the next iteration, these
     * correspond to the eds of the highest centre and its descendants,
     * where eds part of a cluster of which the centre is below the
     * lower bound are updated. These values will be inserted during the
     * method
     * @returns The occurrence corresponding to the upper cluster centre
     * which has a valid distance
     */
    IndexOcc getClusterCentra(uint lowerBound, std::vector<IndexPosExt>& desc,
                           std::vector<uint>& initEds) {
    
        desc.reserve(eds.size());
        initEds.reserve(eds.size());
        IndexOcc m;
        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] > maxED || eds[i] < lowerBound) {
                continue;
            }
            bool betterThanParent = (i == 0) || eds[i] <= eds[i - 1];
            bool betterThanChild = (i == lastCell) || eds[i] <= eds[i + 1];

            if (betterThanParent && betterThanChild) {
                // this is a valid centre
                nodes[i].report(m, startDepth, eds[i], false, shift);

                // get all the descendants
                initEds.emplace_back(eds[i]);
                for (length_t j = i + 1; j <= lastCell; j++) {
                    desc.emplace_back(nodes[j]);
                    initEds.emplace_back(eds[j]);
                }

                // replace the clusters under the lower bound
                for (length_t k = 1; k < initEds.size(); k++) {
                    if (initEds[k] < lowerBound && initEds[k] <= initEds[k - 1] &&
                        (k == initEds.size() - 1 || initEds[k] <= initEds[k + 1])) {
                        // k is a centre under the lower bound

                        length_t highestPoint = 0;
                        length_t lowestPoint = initEds.size() - 1;
                        // find highest point of this cluster
                        for (length_t l = k; l-- > 0;) {
                            if (initEds[l] != initEds[l + 1] + 1) {
                                highestPoint = l + 1;
                                break;
                            }
                        }
                        // find lowest point of this cluster
                        for (length_t l = k + 1; l < initEds.size(); l++) {
                            if (initEds[l] != initEds[l - 1] + 1) {
                                lowestPoint = l - 1;
                                break;
                            }
                        }

                        // highest and lowest cannot span entire
                        // initEds.size(), otherwise there would not be a
                        // valid cluster centre above the lower bound
                        if (highestPoint != 0 &&
                            lowestPoint != initEds.size() - 1) {
                            // Make /\ with ed values of this cluster
                            // do iE[hp] = ie[hp - 1] + 1 and iE[lp] = iE[lp
                            // + 1] +1 until entire cluster has been
                            // replaced
                            length_t lC = lowestPoint;
                            length_t hC = highestPoint;
                            bool highest = true;
                            // do not go over maxED + 1, to ensure
                            // continuity at the other end
                            while (lC > hC) {
                                if (highest) {
                                    initEds[hC] =
                                        std::min(maxED + 1, initEds[hC - 1] + 1);
                                    hC++;
                                } else {
                                    initEds[lC] =
                                        std::min(maxED + 1, initEds[lC + 1] + 1);
                                    lC--;
                                }
                                highest = !highest;
                            }
                            if (lC == hC) {
                                // change middle element of cluster
                                initEds[lC] = std::min(initEds[lC + 1] + 1,
                                                    initEds[lC - 1] + 1);
                            }

                        } else if (highestPoint == 0 &&
                                lowestPoint != initEds.size() - 1) {
                            // monotonous rise from lowestPoint to
                            // highestPoint
                            for (length_t l = lowestPoint; l-- > 0;) {
                                initEds[l] = initEds[l + 1] + 1;
                            }
                        } else if (highestPoint != 0 &&
                                lowestPoint == initEds.size() - 1) {
                            // monotonous rise from highestPoint to
                            // lowestPoint
                            for (length_t l = highestPoint; l < initEds.size();
                                l++) {
                                initEds[l] = initEds[l - 1] + 1;
                            }
                        }
                    }
                }
                // stop searching
                break;
            }
        }

        return m;
    }
};

// ============================================================================
// STRUCT COUNTERS
// ============================================================================

// A struct of performance counters
struct Counters {
    // performance counters
    uint64_t nodeCounter; // counts the number of nodes visited in the index

    uint64_t totalReportedPositions; // counts the number of matches
                                     // reported

    /**
     * Reset all counters to 0
     */
    void resetCounters() {
        nodeCounter = 0,
        totalReportedPositions = 0;
    }

    Counters() {
        resetCounters();
    }

    void addCounters(const Counters& o) {
        nodeCounter += o.nodeCounter,
            totalReportedPositions += o.totalReportedPositions;
    }
    void incNodeCounter() {
        nodeCounter++;
    }
};

// ============================================================================
// CLASS Occurrences
// ============================================================================

// This class combines the occurrences in the text and the occurrences
// in an index into one datastructure
template<class IndexOcc, class IndexPosExt>
class Occurrences {
  private:
    std::vector<TextOcc> inTextOcc; // the in-text occurrences
    std::vector<IndexOcc> inIndexOcc;     // the in-index occurrences

  public:
    Occurrences(const length_t reserve = 200) {
        inTextOcc.reserve(reserve);
        inIndexOcc.reserve(reserve);
    }

    /**
     * Add an index occurrence
     */
    void addIndexOcc(const IndexOcc& match) {
        inIndexOcc.emplace_back(match);
    }

    /**
     * Add an index occurrence
     */
    template<typename... Args> void addIndexOcc(Args... args) {
        inIndexOcc.emplace_back(args...);
    }

    /**
     * Add an index occurrence
     */
    void addIndexOcc(const IndexPosExt& currentNode, const length_t& score) {
        inIndexOcc.emplace_back(currentNode, score);
    }

    /**
     * Add an in-text occurrence
     */
    void addTextOcc(const TextOcc& textOcc) {
        inTextOcc.emplace_back(textOcc);
    }

    /**
     * Add an in-text occurrence
     */
    void addTextOcc(const Range& range, const length_t score) {
        inTextOcc.emplace_back(range, score);
    }

    /**
     * Erase all double in-index occurrences and sorts the occurrences
     */
    void eraseDoublesIndex() {
        sort(inIndexOcc.begin(), inIndexOcc.end());
        inIndexOcc.erase(unique(inIndexOcc.begin(), inIndexOcc.end()), inIndexOcc.end());
    }

    /**
     * Erase all double in-text occurrences and sorts the occurrences
     */
    void eraseDoublesText() {
        sort(inTextOcc.begin(), inTextOcc.end());
        inTextOcc.erase(unique(inTextOcc.begin(), inTextOcc.end()),
                        inTextOcc.end());
    }

    const std::vector<IndexOcc>& getIndexOccurrences() const {
        return inIndexOcc;
    };

    size_t textOccSize() const {
        return inTextOcc.size();
    }

    const std::vector<TextOcc>& getTextOccurrences() const {
        return inTextOcc;
    }

    length_t getMaxSize() const {
        return std::max(inTextOcc.size(), inIndexOcc.size());
    }
};


#endif // INDEXHELPERS_H