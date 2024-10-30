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

#include "rindexhelpers.h"
#include "logger.h"     // for logger
#include <cstdint>      // for uint16_t, uint64_t
#include <fmt/core.h>   // for format
#include <fmt/format.h> // for to_string

using namespace std;

// ============================================================================
// CLASS RANGE
// ============================================================================

ostream& operator<<(ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

// ============================================================================
// CLASS MOVERANGE
// ============================================================================

ostream& operator<<(ostream& os, const MoveRange& r) {
    os << "[" << r.begin << ", " << r.end << ") in runs [" << r.beginRun << ", "
       << r.endRun << "]";
    return os;
}

// ============================================================================
// CLASS TEXT OCC
// ============================================================================
void TextOcc::generateSAMSingleEnd(const string& seqID, const string& printSeq,
                                   const string& printQual, length_t nHits,
                                   length_t minScore, bool primaryAlignment,
                                   const vector<string>& seqNames) {

#ifndef NO_REPORT

    // Estimate the length of the resulting string and reserve memory
    size_t estimated_length =
        seqID.length() + printSeq.length() + printQual.length() + 100;

    outputLine.reserve(estimated_length);

    // Precompute constant values to avoid repeated calculations
    uint16_t flags = getFlagsSE(primaryAlignment);
    int16_t mapQ = getMapQ(nHits, minScore);
    length_t pos = range.getBegin() + 1; // SAM is 1-based

    // Format the output line
    outputLine =
        fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tNM:i:{}"
                    "\tPG:Z:b-move",
                    seqID,                        // read name
                    flags,                        // sam flags
                    seqNames[assignedSequenceID], // reference sequence name
                    pos,         // 1-based pos in reference sequence
                    mapQ,        // mapping quality
                    stringCIGAR, // CIGAR string
                    printSeq,    // sequence or *
                    printQual,   // quality or *
                    distance,    // AS:i: distance
                    distance     // NM:i: distance
        );
#endif
}

void TextOcc::generateSAMSingleEndXA(const string& seqID,
                                     const string& printSeq,
                                     const string& printQual, length_t nHits,
                                     const vector<TextOcc>& otherMatches,
                                     const vector<string>& seqNames) {

#ifndef NO_REPORT

    generateSAMSingleEnd(seqID, printSeq, printQual, nHits, distance, true,
                         seqNames);
    // add the X0 X1 and XA tag
    length_t x0 = nHits - 1; // number of co-optimal hits (nHits includes this)
    length_t x1 = otherMatches.size() - x0; // number of suboptimal hits

    // append the sam line
    outputLine += "\tX0:i:" + fmt::to_string(x0) +
                  "\tX1:i:" + fmt::to_string(x1) + "\tXA:Z:";
    for (const auto& m : otherMatches) {
        outputLine += m.asXA(seqNames);
    }
#endif
}

TextOcc TextOcc::createUnmappedSAMOccurrenceSE(const ReadBundle& bundle) {
    TextOcc t;
    // Estimate the final size and reserve space to avoid multiple
    // allocations
    size_t estimatedSize = bundle.getSeqID().length() +
                           bundle.getRead().length() +
                           bundle.getQual().length() +
                           50; // 50 is a rough estimate for the fixed parts

    t.outputLine.reserve(estimatedSize);

    // Construct the SAM line directly
    t.outputLine =
        fmt::format("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tPG:Z:b-move",
                    bundle.getSeqID(), // Read name
                    bundle.getRead(),  // Read sequence
                    bundle.getQual()   // Read quality
        );

    return t;
}

// ============================================================================
// CLASS FM OCC
// ============================================================================

ostream& operator<<(ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

// ============================================================================
// CLASS CLUSTER
// ============================================================================

FMOcc Cluster::getClusterCentra(uint16_t lowerBound, vector<FMPosExt>& desc,
                                vector<uint16_t>& initEds) {
    desc.reserve(eds.size());
    initEds.reserve(eds.size());
    FMOcc m;
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
                                    min(maxED + 1, initEds[hC - 1] + 1);
                                hC++;
                            } else {
                                initEds[lC] =
                                    min(maxED + 1, initEds[lC + 1] + 1);
                                lC--;
                            }
                            highest = !highest;
                        }
                        if (lC == hC) {
                            // change middle element of cluster
                            initEds[lC] =
                                min(initEds[lC + 1] + 1, initEds[lC - 1] + 1);
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

// ============================================================================
// CLASS COUNTERS
// ============================================================================

void Counters::reportStatistics(const SequencingMode& sMode) const {
    std::stringstream ss;

    ss << "Average no. nodes: "
       << counters[NODE_COUNTER] / (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. Nodes: " << counters[NODE_COUNTER];
    logger.logDeveloper(ss);

    ss << "Average no. unique matches per read: "
        << counters[TOTAL_UNIQUE_MATCHES] /
                (counters[NUMBER_OF_READS] * 1.0);
    logger.logInfo(ss);

    ss << "Total no. matches: " << counters[TOTAL_UNIQUE_MATCHES];
    logger.logInfo(ss);

    ss << "Average no. matches per read "
        << counters[TOTAL_REPORTED_POSITIONS] /
                (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. reported matches: "
        << counters[TOTAL_REPORTED_POSITIONS];
    logger.logDeveloper(ss);

    ss << "Mapped reads: " << counters[MAPPED_READS];
    logger.logInfo(ss);

    ss << "Number of reads: " << counters[NUMBER_OF_READS];
    logger.logInfo(ss);

    ss << "Percentage reads mapped: "
        << (counters[MAPPED_READS] * 100.0) / counters[NUMBER_OF_READS]
        << "%";
    logger.logInfo(ss);

#ifdef LF_BENCHMARK_FUNCTIONALITY
    ss << "Number of LF queries: " << LF_call_count;
    logger.logInfo(ss);

    ss << "Average number of CPU cycles per LF query: "
        << (elapsed_LF * 1.0) / LF_call_count;
    logger.logInfo(ss);
#endif

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    ss << "Number of phi queries: " << phi_call_count;
    logger.logInfo(ss);

    ss << "Average number of CPU cycles per phi query: "
        << (elapsed_phi * 1.0) / phi_call_count;
    logger.logInfo(ss);

    ss << "Total locate time: " << fixed << elapsed_locate.count() << "s";
    logger.logInfo(ss);

    ss << "Number of locate queries: " << locate_call_count;
    logger.logInfo(ss);

    ss << "Average locate time: "
        << (elapsed_locate.count() * 1000.0) / locate_call_count << "ms";
    logger.logInfo(ss);
#endif
}
