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

#include "rsearchstrategy.h"
#include "logger.h"    // for Logger, logger
#include "substring.h" // for Substring

#include <sstream> // for char_traits, basic_istream, ifstream, strings...
using namespace std;

// ============================================================================
// CLASS SEARCH STRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTOR
// ----------------------------------------------------------------------------

SearchStrategy::SearchStrategy(IndexInterface& argument, PartitionStrategy p,
                               DistanceMetric distanceMetric, MappingMode m,
                               SequencingMode sequencingMode)
    : index(argument), distanceMetric(distanceMetric), partitionStrategy(p),
      mode(m) {

    // set the partition strategy
    switch (p) {
    case UNIFORM:
        partitionPtr = &SearchStrategy::partitionUniform;
        break;
    case DYNAMIC:
        partitionPtr = &SearchStrategy::partitionDynamic;
        break;
    case STATIC:
        partitionPtr = &SearchStrategy::partitionOptimalStatic;
        break;
    default:
        break;
    }

    // set the distance metric
    switch (distanceMetric) {
    case HAMMING:
        startIdxPtr = &SearchStrategy::startIndexHamming;
        filterPtr = &SearchStrategy::filterHamming;
        break;
    case EDIT:
        startIdxPtr = &SearchStrategy::startIndexEdit;
        filterPtr = &SearchStrategy::filterEditWithoutCIGARCalculation;

    default:
        break;
    }

    // set the mode
    setMappingMode(m);
}

// ----------------------------------------------------------------------------
// INFORMATION
// ----------------------------------------------------------------------------
string SearchStrategy::getPartitioningStrategy() const {
    switch (partitionStrategy) {
    case UNIFORM:
        return "UNIFORM";
        break;
    case DYNAMIC:
        return "DYNAMIC";
        break;
    case STATIC:
        return "STATIC";
        break;

    default:
        // should not get here
        return "";
    }
}

string SearchStrategy::getDistanceMetric() const {
    switch (distanceMetric) {
    case HAMMING:
        return "HAMMING";
        break;
    case EDIT:
        return "EDIT";
        break;

    default:
        // should not get here
        return "";
    }
}

string SearchStrategy::getMappingModeString() const {
    switch (mode) {
    case BEST:
        return "BEST";
        break;
    case ALL:
        return "ALL";
        break;

    default:
        // should not get here
        return "";
    }
}

// ----------------------------------------------------------------------------
// PARTITIONING
// ----------------------------------------------------------------------------
void SearchStrategy::partition(const string& pattern, vector<Substring>& parts,
                               const int& numParts, const int& maxScore,
                               vector<SARangePair>& exactMatchRanges,
                               Counters& counters) const {
    assert(exactMatchRanges.size() == (uint32_t)numParts);
    parts.clear();

    if (numParts >= (int)pattern.size() || numParts == 1) {
        // no need of splitting up since all parts would be one
        // character or less or there is only one part
        return;
    }
    parts.reserve(numParts);
    (this->*partitionPtr)(pattern, parts, numParts, maxScore, exactMatchRanges,
                          counters);
}

void SearchStrategy::calculateExactMatchRanges(
    vector<Substring>& parts, vector<SARangePair>& exactMatchRanges,
    Counters& counters) const {

    index.setDirection(FORWARD, false); // Bidirectional ranges
    length_t wordSize = index.getWordSize();

    // Match exact ranges for each part bidirectionally except the last
    for (length_t i = 0; i < parts.size() - 0; ++i) {
        auto& current = parts[i];
        auto size = current.size();
        auto start = current.begin() + ((size >= wordSize) ? wordSize : 0);
        auto initRanges =
            (size >= wordSize)
                ? index.lookUpInKmerTable({current, current.begin(), start})
                : index.getCompleteRange();
        exactMatchRanges[i] = index.matchStringBidirectionally(
            {current, start, current.end()}, initRanges, counters);
    }

    // Last part matched unidirectional backwards
    index.setDirection(BACKWARD, true);
    auto& lastPart = parts.back();
    lastPart.setDirection(BACKWARD);
    auto size = lastPart.size();
    auto end = (size >= wordSize) ? lastPart.end() - wordSize : lastPart.end();
    auto initRanges =
        (size >= wordSize)
            ? index.lookUpInKmerTable({lastPart, end, lastPart.end()})
            : index.getCompleteRange();
    exactMatchRanges.back() = index.matchStringBidirectionally(
        {lastPart, lastPart.begin(), end}, initRanges, counters);
}

// Uniform Partitioning

void SearchStrategy::partitionUniform(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<SARangePair>& exactMatchRanges,
                                      Counters& counters) const {

    for (int i = 0; i < numParts; i++) {
        parts.emplace_back(pattern, (i * 1.0 / numParts) * pattern.size(),
                           ((i + 1) * 1.0 / numParts) * pattern.size());
    }
    // set end of final part correct
    parts.back().setEnd(pattern.size());

    // calculate the exact match ranges
    calculateExactMatchRanges(parts, exactMatchRanges, counters);
}

// Static Partitioning
void SearchStrategy::partitionOptimalStatic(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<SARangePair>& exactMatchRanges,
    Counters& counters) const {

    setParts(pattern, parts, numParts, maxScore);
    calculateExactMatchRanges(parts, exactMatchRanges, counters);
}

void SearchStrategy::setParts(const string& pattern, vector<Substring>& parts,
                              const int& numParts, const int& maxScore) const {
    const vector<double>& begins = getBegins(numParts, maxScore);

    int pSize = pattern.size();
    // set the first part
    parts.emplace_back(pattern, 0, begins[0] * pSize);

    for (unsigned int i = 0; i < begins.size() - 1; i++) {
        parts.emplace_back(pattern, begins[i] * pSize, begins[i + 1] * pSize);
    }
    parts.emplace_back(pattern, begins.back() * pSize, pattern.size());
}

// Dynamic Partitioning

void SearchStrategy::partitionDynamic(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<SARangePair>& exactMatchRanges,
                                      Counters& counters) const {

    int matchedChars =
        seed(pattern, parts, numParts, maxScore, exactMatchRanges);

    int pSize = pattern.size();
    vector<int> weights = getWeights(numParts, maxScore);

    Direction dir = FORWARD;
    int partToExtend = 0;

    // extend the part with the largest range, as to minimize the range
    // for each part do this until all characters are assigned to a
    // part
    for (int j = matchedChars; j < pSize; j++) {

        // find the part with the largest range
        length_t maxRangeWeighted = 0;

        for (int i = 0; i < numParts; i++) {
            bool noLeftExtension =
                (i == 0) || parts[i].begin() == parts[i - 1].end();
            bool noRightExtension =
                (i == numParts - 1) || parts[i].end() == parts[i + 1].begin();
            if (noLeftExtension && noRightExtension) {
                continue;
            }
            if (exactMatchRanges[i].width() * weights[i] > maxRangeWeighted) {
                maxRangeWeighted = exactMatchRanges[i].width() * weights[i];
                partToExtend = i;
                if (noLeftExtension) {
                    // only right extension
                    dir = FORWARD;
                } else if (noRightExtension) {
                    // only left extension
                    dir = BACKWARD;
                } else {
                    // both directions possible, choose direction of
                    // smallest neighbour
                    dir = (exactMatchRanges[i - 1].width() <
                           exactMatchRanges[i + 1].width())
                              ? BACKWARD
                              : FORWARD;
                }
            }
        }

        if (maxRangeWeighted == 0) {
            // no need to keep calculating new range, just extend the
            // parts
            extendParts(pattern, parts);
            return;
        }

        // extend partToExtend in direction
        char c; // the new character
        if (dir == FORWARD) {
            parts[partToExtend].incrementEnd();
            c = pattern[parts[partToExtend].end() - 1];
        } else {
            parts[partToExtend].decrementBegin();
            c = pattern[parts[partToExtend].begin()];
        }

        // match the new character, if this is the last part use unidirectional
        // backward matching
        index.setDirection(dir, partToExtend == numParts - 1);
        index.addChar(c, exactMatchRanges.at(partToExtend), counters);
    }
}

int SearchStrategy::seed(const string& pattern, vector<Substring>& parts,
                         const int& numParts, const int& maxScore,
                         vector<SARangePair>& exactMatchRanges) const {
    auto pSize = pattern.size();
    bool useKmerTable = (numParts * index.getWordSize() < (pSize * 2) / 3);
    int wSize = (useKmerTable) ? index.getWordSize() : 1;

    const auto& seedPercent = getSeedingPositions(numParts, maxScore);

    vector<int> seeds;
    // push the seed for the first part
    seeds.emplace_back(0);

    // push the optimal seeds for the middle parts
    for (int i = 1; i < numParts - 1; i++) {
        seeds.emplace_back((seedPercent[i - 1] * pSize) - (wSize / 2));
    }

    for (int i = 0; i < numParts - 1; i++) {
        parts.emplace_back(pattern, seeds[i], seeds[i] + wSize);
    }

    // push the seeds for the final parts
    parts.emplace_back(pattern, pSize - wSize, pSize);

    exactMatchRanges.resize(numParts);
    for (int i = 0; i < numParts; i++) {

        exactMatchRanges[i] = (useKmerTable)
                                  ? index.lookUpInKmerTable(parts[i])
                                  : index.getRangeOfSingleChar(parts[i][0]);
    }
    return numParts * wSize;
}

void SearchStrategy::extendParts(const string& pattern,
                                 vector<Substring>& parts) const {
    for (length_t i = 0; i < parts.size(); i++) {

        if ((i != parts.size() - 1) &&
            (parts[i].end() != parts[i + 1].begin())) {
            // extend completely to the right
            // it is known that the range will stay [0,0)
            parts[i].setEnd(parts[i + 1].begin());
        }
        if ((i != 0) && (parts[i].begin() != parts[i - 1].end())) {
            // extend completely to the left
            parts[i].setBegin(parts[i - 1].end());
        }
    }
}

// ----------------------------------------------------------------------------
// (APPROXIMATE) MATCHING
// ----------------------------------------------------------------------------

void SearchStrategy::matchWithSearches(const string& seq, const length_t k,
                                       Counters& counters, Occurrences& occs,
                                       const length_t minDistance) {
    // calculate number of parts
    length_t numParts = calculateNumParts(k);

    // create ranges for exact matches of parts
    vector<SARangePair> exactMatchRanges(numParts);

    // create a vector for the parts
    vector<Substring> parts;

    // partition the read
    partition(seq, parts, numParts, k, exactMatchRanges, counters);

    if (parts.empty()) {
        // splitting up was not viable -> just search the entire pattern
        if (name != "Naive backtracking") {
            stringstream ss;
            ss << "Normal bidirectional search was used as "
                  "entered pattern is too short "
               << seq.size();
            logger.logWarning(ss);
        }

        vector<TextOcc> textOccs;

        index.approxMatchesNaive(seq, k, counters, textOccs);
        occs.addTextOccs(textOccs);
        return;
    }

    // B) do the searches from search schemes
    // create the searches
    const vector<Search>& searches = createSearches(k, exactMatchRanges);

    // reserve stacks and matrices for each part
    index.reserveStacks(numParts,
                        seq.length()); // reserve stacks for each part
    // create the bit-parallel alignment matrices
    index.resetMatrices(parts.size()); // reset the alignment matrix that will
                                       // be (possibly) used for each part

    for (const Search& s : searches) {
        doRecSearch(s, parts, occs, exactMatchRanges, counters);
    }
}

void SearchStrategy::matchApproxAllMap(ReadBundle& bundle, length_t maxED,
                                       Counters& counters,
                                       vector<TextOcc>& result) {

    if (maxED == 0) {

        index.setIndexInMode(FORWARD_STRAND, FIRST_IN_PAIR);
        index.exactMatchesOutput(bundle.getRead(), counters, result);
        index.setIndexInMode(REVERSE_C_STRAND, FIRST_IN_PAIR);
        index.exactMatchesOutput(bundle.getRevComp(), counters, result);

        if (generateSAMLines) {
            generateOutputSingleEnd(result, bundle, counters, 0);
        }
        return;
    }

    // The occurrences in the text and index
    Occurrences occ;

    // set the index in forward mode with the correct sequence
    index.setIndexInMode(FORWARD_STRAND);
    // map with the searches
    matchWithSearches(bundle.getRead(), maxED, counters, occ);

    // set the index in reverse complement mode with the correct sequence
    index.setIndexInMode(REVERSE_C_STRAND);
    // map with the searches
    matchWithSearches(bundle.getRevComp(), maxED, counters, occ);

    // C) get all non-redundant matches mapped to the text
    result = (this->*filterPtr)(occ, maxED, counters, bundle);

    // D) generate sam output
    if (generateSAMLines) {
        generateOutputSingleEnd(result, bundle, counters, maxED);
    }
}

void SearchStrategy::checkAlignments(OccVector& occVector, uint32_t& best,
                                     uint32_t l, Counters& counters,
                                     uint32_t cutOff, const string& seq) const {
    std::vector<TextOcc> trimmedOccs = {};
    std::vector<TextOcc> assignedOccs = {};
    for (length_t i = 0; i < occVector[l].second.size(); i++) {
        auto& occ = occVector[l].second[i];
        auto seqFound = assignSequence(occ, counters, cutOff, seq);
        if (seqFound != FOUND) {
            if (seqFound == FOUND_WITH_TRIMMING && occ.getDistance() > l) {
                trimmedOccs.emplace_back(std::move(occ));
            }
        } else {
            assignedOccs.emplace_back(std::move(occ));
            if (l < best) {
                best = l;
            }
        }
    }
    occVector[l].second = std::move(assignedOccs);
    // the occurrences found with trimming should be added
    // to the correct vector according to their new distance
    for (auto& occ : trimmedOccs) {
        occ.removeTrimmingLabel();
        occVector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
}

vector<TextOcc> SearchStrategy::combineOccVectors(OccVector& ovFW,
                                                  OccVector& ovRC,
                                                  length_t best, length_t max) {
    assert(ovRC.size() == ovFW.size());
    assert(max < ovFW.size());

    auto compare = [](const TextOcc& a, const TextOcc& b) {
        return a.getAssignedSequenceID() < b.getAssignedSequenceID() ||
               (a.getAssignedSequenceID() == b.getAssignedSequenceID() &&
                a.getRange().getBegin() < b.getRange().getBegin());
    };
    auto equal = [](const TextOcc& a, const TextOcc& b) {
        return a.getAssignedSequenceID() == b.getAssignedSequenceID() &&
               a.getRange().getBegin() == b.getRange().getBegin();
    };

    std::vector<TextOcc> matches;
    matches.reserve(numElements(ovFW) + numElements(ovRC));

    for (length_t i = best; i <= max; i++) {
// some occurrences might be represented more than once
// sort the occurrences first on assigned sequence name and then on
// begin position in the text
#ifdef DEVELOPER_MODE
        std::stable_sort(ovFW[i].second.begin(), ovFW[i].second.end(), compare);
        std::stable_sort(ovRC[i].second.begin(), ovRC[i].second.end(), compare);
#else
        std::sort(ovFW[i].second.begin(), ovFW[i].second.end(), compare);
        std::sort(ovRC[i].second.begin(), ovRC[i].second.end(), compare);
#endif

        // remove duplicates
        // Remove duplicates from ovFW[i].second
        auto fw_end =
            std::unique(ovFW[i].second.begin(), ovFW[i].second.end(), equal);
        ovFW[i].second.erase(fw_end, ovFW[i].second.end());

        // Remove duplicates from ovRC[i].second
        auto rc_end =
            std::unique(ovRC[i].second.begin(), ovRC[i].second.end(), equal);
        ovRC[i].second.erase(rc_end, ovRC[i].second.end());

        matches.insert(matches.end(),
                       make_move_iterator(ovFW[i].second.begin()),
                       make_move_iterator(ovFW[i].second.end()));
        matches.insert(matches.end(),
                       make_move_iterator(ovRC[i].second.begin()),
                       make_move_iterator(ovRC[i].second.end()));
    }

    return matches;
}

void SearchStrategy::matchApproxBestPlusX(ReadBundle& bundle, length_t x,
                                          Counters& counters,
                                          const length_t minIdentity,
                                          vector<TextOcc>& result) {

    length_t cutOff = getMaxED(minIdentity, bundle.size());
    uint32_t best = cutOff + 1;
    bool bestFound = false;

    // return value
    // create a vector of vectors of TextOcc per edit distance
    OccVector occVectorFW(cutOff + 1), occVectorRC(cutOff + 1);

    // get the read and the revCompl from the bundle
    const auto& read = bundle.getRead();
    const auto& revCompl = bundle.getRevComp();

    // start with exact match, this does not need an in-text verification matrix
    // as an exact match has no insertions or deletions
    index.setIndexInMode(FORWARD_STRAND);
    index.exactMatchesOutput(bundle.getRead(), counters, occVectorFW[0].second);
    index.setIndexInMode(REVERSE_C_STRAND);
    index.exactMatchesOutput(bundle.getRevComp(), counters,
                             occVectorRC[0].second);
    occVectorFW[0].first = true,
    occVectorRC[0].first = true; // processing for 0 finished

    if (!occVectorFW[0].second.empty() || !occVectorRC[0].second.empty()) {
        checkAlignments(occVectorFW, best, 0, counters, cutOff, read);
        checkAlignments(occVectorRC, best, 0, counters, cutOff, revCompl);
        if (best == 0) {
            bestFound = true;
        }
    }

    // moving to approximate pattern matching if necessary
    length_t maxED = (best == 0) ? x : cutOff;
    length_t prevK = 0;

    for (length_t k = std::max(x, (length_t)1); k <= maxED;) {
        bool updated = false;
        // start with forward direction
        updated |= processSeq(read, FIRST_IN_PAIR, FORWARD_STRAND, k,
                              occVectorFW, counters);
        // then reverse complement
        updated |= processSeq(revCompl, FIRST_IN_PAIR, REVERSE_C_STRAND, k,
                              occVectorRC, counters);

        if (updated) {
            // new matches were found, we need to process them to see if
            // they need trimming
            for (length_t l = prevK + 1; l <= std::min(k, best + x); l++) {
                checkAlignments(occVectorFW, best, l, counters, maxED, read);
                checkAlignments(occVectorRC, best, l, counters, maxED,
                                revCompl);
            }
        }
        if (bestFound) {
            break; // this was the last iteration
        }
        if (updated && best < cutOff + 1) {
            bestFound = true;
            if (x == 0) {
                break;
            }
            prevK = k,
            k = std::min(best + x, maxED); // check the final x strata
        } else {
            // continue to next stratum
            if (k == maxED)
                break;
            length_t step = (k < 5) ? 2 : 4;
            prevK = k;
            k = std::min(k + x + step, maxED);
        }
    }

    if (!bestFound) {
        if (unmappedSAM) {
            result.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }

    // find out the number of hits for best
    uint32_t nHits =
        occVectorFW[best].second.size() + occVectorRC[best].second.size();

    result = combineOccVectors(occVectorFW, occVectorRC, best,
                               std::min(best + x, maxED));

    (this->*generateOutputSEPtr)(bundle, nHits, best, result, counters);
}

bool SearchStrategy::processSeq(const string& seq, const PairStatus status,
                                const Strand strand, const length_t maxDist,
                                OccVector& vector, Counters& counters) {
    assert(maxDist < vector.size());
    // if maxDist has been processed then all distances below it have
    // also been processed
    if (!vector[maxDist].first) {

        // find the minimum distance (first distance d for which
        // vector[d].first == false)
        auto it =
            find_if(vector.begin(), vector.end(),
                    [](const BoolAndVector& pair) { return !pair.first; });
        length_t minD = min((length_t)distance(vector.begin(), it), maxDist);

        auto occs = mapRead(seq, maxDist, counters, status, strand, minD);

        // the mapped occurrences can be between the min and max
        // distance they need to be split up according to their
        // distances and added to the correct vector
        for (auto& occ : occs) {
            vector[occ.getDistance()].second.emplace_back(std::move(occ));
        }
        // set distance processed to true
        for (length_t i = minD; i <= maxDist; i++) {
            vector[i].first = true;
        }
    }
    // return true if any of the distances between 0 and maxDist have an
    // occurrence
    bool found =
        any_of(vector.begin(), vector.begin() + maxDist + 1,
               [](const BoolAndVector& pair) { return !pair.second.empty(); });
    return found;
}

void SearchStrategy::handleTrimmedOccs(vector<length_t>& trimmedIds,
                                       const length_t oDist, OccVector& v) {

// sort the vector by id
#ifdef DEVELOPER_MODE
    stable_sort(trimmedIds.begin(), trimmedIds.end());
#else
    sort(trimmedIds.begin(), trimmedIds.end());
#endif
    for (auto idIt = trimmedIds.rbegin(); idIt != trimmedIds.rend(); ++idIt) {
        auto id = *idIt;
        auto& occ = v[oDist].second[id];
        occ.removeTrimmingLabel();
        // insert occ in the correct vector and remove from old vector
        auto& tVector = v[occ.getDistance()].second;
        auto point = lower_bound(tVector.begin(), tVector.end(), occ);
        tVector.insert(point, std::move(occ));
        v[oDist].second.erase(v[oDist].second.begin() + id);
    }
}

void SearchStrategy::addSingleEndedForBest(vector<TextOcc>& matchesSE1,
                                           vector<TextOcc>& matchesSE2,
                                           OccVector& fw1, OccVector& rc1,
                                           OccVector& fw2, OccVector& rc2,
                                           bool read2done) {

    for (length_t i = 0; i < matchesSE1.size(); i++) {
        auto occ = std::move(matchesSE1[i]);
        auto& vector = (occ.getStrand() == FORWARD_STRAND) ? fw1 : rc1;
        vector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
    for (length_t i = 0; i < matchesSE2.size(); i++) {
        auto occ = std::move(matchesSE2[i]);
        occ.setPairStatus(SECOND_IN_PAIR); // set flag for  second read
        auto& vector = (occ.getStrand() == FORWARD_STRAND) ? fw2 : rc2;
        vector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
    for (length_t i = 0; i < fw1.size(); i++) {
        fw1[i].first = rc1[i].first = true;
    }
    for (length_t i = 0; i < fw2.size(); i++) {
        fw2[i].first = rc2[i].first = read2done;
    }
}

void SearchStrategy::doRecSearch(const Search& s, vector<Substring>& parts,
                                 Occurrences& occ,
                                 const vector<SARangePair>& exactMatchRanges,
                                 Counters& counters) const {

    if (s.getUpperBound(0) > 0) {
        // first part is allowed an error so start with an empty match
        s.setDirectionsInParts(parts);

        SARangePair startRange = index.getCompleteRange();
        FMOcc startMatch = FMOcc(startRange, 0, 0);
        (this->*startIdxPtr)(s, startMatch, occ, parts, counters, 0);
        return;
    }

    // first get the bidirectional match of first part
    int first = s.getPart(0);
    SARangePair startRange = exactMatchRanges[first];

    // if this range is bigger than the switch point
    if (startRange.width() > 0) {

        // prepare the parts for this search
        s.setDirectionsInParts(parts);

        // can we continue exact matching according to this search?
        uint16_t partInSearch = 1;
        length_t exactLength = parts[first].size();

        while (s.getUpperBound(partInSearch) == 0) {
            // extend the exact match
            index.setDirection(s.getDirection(partInSearch),
                               s.isUnidirectionalBackwards(partInSearch));
            const auto& part = parts[s.getPart(partInSearch)];

            startRange =
                index.matchStringBidirectionally(part, startRange, counters);
            if (startRange.empty()) {
                return;
            }

            exactLength += part.size();
            partInSearch++;
        }

        // Create a match corresponding to the exact match
        FMOcc startMatch = FMOcc(startRange, 0, exactLength);

        // Start approximate matching in the index
        (this->*startIdxPtr)(s, startMatch, occ, parts, counters, partInSearch);
    }
}

// ----------------------------------------------------------------------------
// PAIRING
// ----------------------------------------------------------------------------

/**
 * @brief Find the first occurrence in a sorted vector of TextOcc that is after
 * a specified position
 * @param occs The vector of TextOcc
 * @param position The position to search for
 * @return An iterator to the first occurrence that is after the specified
 * position
 */
vector<TextOcc>::iterator findFirstOccAfter(vector<TextOcc>& occs,
                                            length_t position) {
    assert(is_sorted(occs.begin(), occs.end()));
    return lower_bound(occs.begin(), occs.end(), position,
                       [](const TextOcc& occ, length_t pos) {
                           return occ.getIndexBegin() < pos;
                       });
}

void assignSequence(TextOcc& occ) {
}

vector<TextOcc> SearchStrategy::findBestMapping(OccVector& fw, OccVector& rc,
                                                const ReadBundle& bundle,
                                                Counters& counters,
                                                PairStatus status) {
    assert(fw.size() == rc.size());
    vector<TextOcc> matches;
    for (length_t i = 0; i < fw.size(); i++) {
        // map the stratum
        mapStratum(fw, i, counters, bundle.getRead(), status, FORWARD_STRAND);
        matches.insert(matches.end(), make_move_iterator(fw[i].second.begin()),
                       make_move_iterator(fw[i].second.end()));
        mapStratum(rc, i, counters, bundle.getRevComp(), status,
                   REVERSE_C_STRAND);
        matches.insert(matches.end(), make_move_iterator(rc[i].second.begin()),
                       make_move_iterator(rc[i].second.end()));
        if (!matches.empty()) {
            return matches;
        }
    }

    return matches;
}
// ----------------------------------------------------------------------------
// POST PROCESSING
// ----------------------------------------------------------------------------

void SearchStrategy::generateOutputSingleEnd(vector<TextOcc>& occs,
                                             ReadBundle& bundle,
                                             Counters& counters,
                                             length_t cutOff) const {

    if (occs.empty()) {
        if (unmappedSAM) {
            occs.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }

    vector<TextOcc> assignedOccs;
    assignedOccs.reserve(occs.size());

    // for each occurrence find the assigned sequence
    for (auto& t : occs) {
        index.setIndexInMode(t.getStrand(), t.getPairStatus());
        const auto& seq = bundle.getSequence(t.getStrand());
        length_t seqID;
        SeqNameFound found =
            index.findSeqName(t, seqID, counters, cutOff, distanceMetric, seq);
        if (found != SeqNameFound::NOT_FOUND) {
            t.setAssignedSequence(found, seqID);
            assignedOccs.emplace_back(std::move(t));
        }
    }

    // keep only the occurrences that have an assigned sequence
    occs = std::move(assignedOccs);

    if (occs.empty()) {
        if (unmappedSAM) {
            occs.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }
#ifdef DEVELOPER_MODE
    // sort the occurrences on score
    // Use stable sort because we expect many elements with same key
    stable_sort(occs.begin(), occs.end(),
                [](const TextOcc& a, const TextOcc& b) {
                    return a.getDistance() < b.getDistance();
                });

    // find the minimal score
    uint32_t minScore = occs.front().getDistance();

    uint32_t nHits = 1;

    // count the number of hits with the minimal score
    // take into account that the occurrences are sorted on score
    // find the index of the first occurrence with a score higher than
    // minScore with a binary search
    auto it = lower_bound(occs.begin(), occs.end(), minScore + 1);
    nHits = distance(occs.begin(), it);

#else
    // Find minimal distance
    auto minElement = std::min_element(
        occs.begin(), occs.end(), [](const TextOcc& a, const TextOcc& b) {
            return a.getDistance() < b.getDistance();
        });

    length_t minScore = minElement->getDistance();

    // Count number of elements with minimal distance
    length_t nHits =
        std::count_if(occs.begin(), occs.end(), [&](const TextOcc& elem) {
            return elem.getDistance() == minScore;
        });

    // Swap element with minimal distance to primary position
    if (minElement != occs.begin()) {
        std::iter_swap(occs.begin(), minElement);
    }
#endif
    (this->*generateOutputSEPtr)(bundle, nHits, minScore, occs, counters);
}

// ============================================================================
// CLASS CUSTOM SEARCH STRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTION
// ----------------------------------------------------------------------------

void CustomSearchStrategy::getSearchSchemeFromFolder(string pathToFolder,
                                                     bool verbose) {

    // get the name of the file
    string line;
    {
        ifstream ifs(pathToFolder + "name.txt");
        if (!ifs) {
            throw runtime_error("Problem reading: " + pathToFolder +
                                "name.txt\nDid you provide a directory to "
                                "a search scheme without a name file?");
        }
        getline(ifs, line);
        name = line;
        ifs.close();
    }

    // get the info per distance score (scores between 1 and MAX_K are
    // looked at)
    for (int i = 1; i <= MAX_K; i++) {
        string name =
            pathToFolder + to_string(i) + PATH_SEPARATOR + "searches.txt";
        ifstream stream_searches(name);
        if (!stream_searches) {
            // this score is not supported
            supportsMaxScore[i - 1] = false;
            continue;
        }

        schemePerED[i - 1] =
            SearchScheme::readScheme(stream_searches, name, i).getSearches();

        if (schemePerED[i - 1].size() > 0) {
            supportsMaxScore[i - 1] = true;
        }
        stream_searches.close();
    }

    // check if the searches are valid
    sanityCheck(verbose);

    // get static positions (if they exist)
    for (int i = 1; i <= MAX_K; i++) {
        ifstream stream_static(pathToFolder + to_string(i) + PATH_SEPARATOR +
                               "static_partitioning.txt");

        if (stream_static) {
            // a file with static partitioning positions exists
            getline(stream_static, line);
            vector<string> positionsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                positionsAsString.push_back(token);
            }

            if (positionsAsString.size() != calculateNumParts(i) - 1) {
                throw runtime_error(
                    "Not enough static positions provided in " + pathToFolder +
                    to_string(i) + PATH_SEPARATOR +
                    "static_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " parts\nProvided: " +
                    to_string(positionsAsString.size()) + " parts");
            }

            for (auto str : positionsAsString) {
                staticPositions[i - 1].push_back(stod(str));
            }

            // check if these positions are valid
            sanityCheckStaticPartitioning(i);
            // if valid set getBeginsPointer to custom
            beginsPointer[i - 1] = &CustomSearchStrategy::getBeginsCustom;

            stream_static.close();
        }
    }

    // get dynamic seeds and weights (if file exists)
    for (int i = 1; i <= MAX_K; i++) {
        ifstream stream_dynamic(pathToFolder + to_string(i) + PATH_SEPARATOR +
                                "dynamic_partitioning.txt");

        if (stream_dynamic) {
            // a file with dynamic partitioning positions exists
            getline(stream_dynamic, line);
            vector<string> seedsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                seedsAsString.push_back(token);
            }

            if (seedsAsString.size() != calculateNumParts(i) - 2) {
                throw runtime_error(
                    "Not enough seeding positions provided in " + pathToFolder +
                    to_string(i) + PATH_SEPARATOR +
                    "dynamic_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " seeds\nProvided: " +
                    to_string(seedsAsString.size()) + " seeds");
            }

            for (auto str : seedsAsString) {
                seedingPositions[i - 1].push_back(stod(str));
            }

            // check if these seeds are valid
            sanityCheckDynamicPartitioning(i);

            // get the weights
            getline(stream_dynamic, line);
            stringstream ss_w(line);
            string stringWeight;
            while (ss_w >> stringWeight) {
                weights[i - 1].push_back(stoi(stringWeight));
            }

            if (weights[i - 1].size() != calculateNumParts(i)) {
                throw runtime_error(
                    "Not enough weights provided for max score " +
                    to_string(i));
            }

            // set the pointers to custom
            seedingPointer[i - 1] =
                &CustomSearchStrategy::getSeedingPositionsCustom;
            weightsPointers[i - 1] = &CustomSearchStrategy::getWeightsCustom;
        }
    }
}

Search CustomSearchStrategy::makeSearch(const string& line,
                                        length_t idx) const {
    stringstream ss(line);

    vector<string> tokens;
    string token;
    while (ss >> token) {
        tokens.push_back(token);
    }

    if (tokens.size() != 3) {
        throw runtime_error("A search should have 3 vectors: order, "
                            "lower bound and upper bound!");
    }

    vector<length_t> order;
    getVector(tokens[0], order);

    vector<length_t> lower_bound;
    getVector(tokens[1], lower_bound);

    vector<length_t> upper_bound;
    getVector(tokens[2], upper_bound);

    return Search::makeSearch(order, lower_bound, upper_bound, idx);
}

void CustomSearchStrategy::getVector(const string& vectorString,
                                     vector<length_t>& vector) const {

    if (vectorString.size() < 2) {
        throw runtime_error(vectorString +
                            " is not a valid vector for a search");
    }
    string bracketsRemoved = vectorString.substr(1, vectorString.size() - 2);

    stringstream ss(bracketsRemoved);
    string token;
    while (getline(ss, token, ',')) {
        vector.emplace_back(stoull(token));
    }
}

// ----------------------------------------------------------------------------
// SANITY CHECKS
// ----------------------------------------------------------------------------

void CustomSearchStrategy::sanityCheckStaticPartitioning(
    const int& maxScore) const {
    const auto& positions = staticPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than
    // 0
    for (unsigned int i = 0; i < positions.size(); i++) {
        if (positions[i] <= 0 || positions[i] >= 1) {
            throw runtime_error("One of the provided static positions for " +
                                to_string(maxScore) +
                                " is not between 0 and 1 (exclusive)");
        }
        if (i < positions.size() - 1 && positions[i] - positions[i + 1] >= 0) {
            throw runtime_error("Provided static positions for " +
                                to_string(maxScore) +
                                " are not strictly increasing");
        }
    }
}
void CustomSearchStrategy::sanityCheckDynamicPartitioning(
    const int& maxScore) const {

    const auto& seeds = seedingPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than
    // 0
    for (unsigned int i = 0; i < seeds.size(); i++) {
        if (seeds[i] <= 0 || seeds[i] >= 1) {
            throw runtime_error("One of the provided static positions for " +
                                to_string(maxScore) +
                                " is not between 0 and 1 (exclusive)!");
        }
        if (i < seeds.size() - 1 && seeds[i] - seeds[i + 1] >= 0) {
            throw runtime_error("Provided seeding positions for " +
                                to_string(maxScore) +
                                " are not strictly increasing");
        }
    }
}

void CustomSearchStrategy::sanityCheck(bool verbose) const {

    // check if for each supported edit distance all error distributions
    // are covered

    for (int K = 1; K <= MAX_K; K++) {

        const auto& scheme = schemePerED[K - 1];
        if (!supportsMaxScore[K - 1]) {
            continue;
        }
        const Search& firstSearch = scheme.front();
        length_t P = firstSearch.getNumParts();
        // check if all searches have same number of parts
        if (any_of(scheme.begin(), scheme.end(),
                   [P](const Search& s) { return s.getNumParts() != P; })) {
            logger.logError("Not all searches for distance " + to_string(K) +
                            " have the same number of "
                            "parts");
            throw runtime_error("Not all searches for "
                                "distance " +
                                to_string(K) +
                                " have the same "
                                "number of parts");
        }

        // check if zero based
        if (any_of(scheme.begin(), scheme.end(),
                   [](const Search& s) { return !s.zeroBased(); })) {
            throw runtime_error(
                "Not all searches are zero based for distance " + to_string(K) +
                "!");
        }

        // check if connectivity satisfied
        if (any_of(scheme.begin(), scheme.end(), [](const Search& s) {
                return !s.connectivitySatisfied();
            })) {
            throw runtime_error("Connectivity property not satisfied "
                                "for all searches with distance " +
                                to_string(K) + "!");
        }

        // check if U and L string are valid
        if (any_of(scheme.begin(), scheme.end(),
                   [](const Search& s) { return !s.validBounds(); })) {
            throw runtime_error("Decreasing lower or upper bounds "
                                "for a search for K  = " +
                                to_string(K));
        }
    }
}
