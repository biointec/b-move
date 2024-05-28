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

#include <sstream> // used for splitting strings

using namespace std;

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTOR
// ----------------------------------------------------------------------------

RSearchStrategy::RSearchStrategy(RIndex& argument, PartitionStrategy p,
                               DistanceMetric distanceMetric)
    : index(argument), partitionStrategy(p), distanceMetric(distanceMetric) {

    // set the partition strategy
    switch (p) {
    case UNIFORM:
        partitionPtr = &RSearchStrategy::partitionUniform;
        break;
    case DYNAMIC:
        partitionPtr = &RSearchStrategy::partitionDynamic;
        break;
    case STATIC:
        partitionPtr = &RSearchStrategy::partitionOptimalStatic;
        break;
    default:
        break;
    }

    // set the distance metric
    switch (distanceMetric) {
    case HAMMING:
        startIdxPtr = &RSearchStrategy::startIndexHamming;
        filterPtr = &RSearchStrategy::filterHamming;
        break;
    case EDITNAIVE:
        startIdxPtr = &RSearchStrategy::startIndexEditNaive;
        filterPtr = &RSearchStrategy::filterEdit;
        break;
    case EDITOPTIMIZED:
        startIdxPtr = &RSearchStrategy::startIndexEditOptimized;
        filterPtr = &RSearchStrategy::filterEdit;
    default:
        break;
    }
}

// ----------------------------------------------------------------------------
// INFORMATION
// ----------------------------------------------------------------------------
string RSearchStrategy::getPartitioningStrategy() const {
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

string RSearchStrategy::getDistanceMetric() const {
    switch (distanceMetric) {
    case HAMMING:
        return "HAMMING";
        break;
    case EDITOPTIMIZED:
        return "(OPTIMIZED) EDIT";
        break;
    case EDITNAIVE:
        return "(NAIVE) EDIT";
        break;

    default:
        // should not get here
        return "";
    }
}
// ----------------------------------------------------------------------------
// SANITY CHECKS
// ----------------------------------------------------------------------------

void RSearchStrategy::genErrorDistributions(
    int P, int K, vector<Distribution>& distributions) {
    Distribution distribution(P, 0);

    for (int i = 0; i < pow(K + 1, P); i++) {
        int sum = accumulate(distribution.begin(), distribution.end(), 0);
        if (sum <= K)
            distributions.push_back(distribution);

        for (int j = 0; j < P; j++) {
            distribution[j]++;
            if (distribution[j] != K + 1)
                break;
            distribution[j] = 0;
        }
    }
}

bool RSearchStrategy::coversDistributions(
    const vector<Distribution>& distributions, const vector<Search>& scheme,
    bool verbose) {
    vector<int> numCover(scheme.size(), 0);

    bool ret = true;

    // and error distribution is simply a vector containing
    // the number of errors in each partition P
    for (const Distribution& distribution : distributions) {
        int numberOfCovers = 0;

        vector<int> searchesThatCover;
        // check if a search covers the distribution
        bool distributionCovered = false;
        for (size_t si = 0; si < scheme.size(); si++) {
            const Search& s = scheme[si];
            bool thisCover = true;
            length_t numErrors = 0;

            // check if search covers distribution
            for (length_t i = 0; i < s.getNumParts(); i++) {
                length_t p = s.getPart(i);
                numErrors += distribution[p];
                if ((numErrors > s.getUpperBound(i)) ||
                    (numErrors < s.getLowerBound(i)))
                    thisCover = false;
            }

            // print the distribution and the search that covers it
            if (thisCover) {
                searchesThatCover.push_back(si + 1);
                numberOfCovers++;
                numCover[si]++;

                if (verbose) {
                    cout << "Pattern: ";
                    for (length_t i = 0; i < distribution.size(); i++)
                        cout << distribution[i];
                    cout << " is covered by search:" << si + 1 << " (";
                    for (length_t i = 0; i < distribution.size(); i++)
                        cout << s.getPart(i);
                    cout << "), (";
                    for (length_t i = 0; i < s.getNumParts(); i++)
                        cout << s.getLowerBound(i);
                    cout << "), (";
                    for (length_t i = 0; i < s.getNumParts(); i++)
                        cout << s.getUpperBound(i);
                    cout << ")";

                    if (distributionCovered)
                        cout << " * "; // * means redundant
                    cout << endl;
                }
                distributionCovered = true;
            }
        }

        if (!distributionCovered && verbose) {
            cout << "Pattern is not covered: ";
            for (length_t i = 0; i < distribution.size(); i++)
                cout << distribution[i];
            cout << endl;
        }

        ret &= distributionCovered;

        if (!ret && !verbose)
            return ret;

        if (verbose) {
            cout << "Pattern ";
            for (length_t i = 0; i < distribution.size(); i++)
                cout << distribution[i];
            cout << " is covered by searches: ";
            for (auto si : searchesThatCover)
                cout << si << ",";
            cout << "\n";
        }
    }

    if (verbose)
        cout << distributions.size() << " distributions covered" << endl;

    // check if all searches cover at least one distribution
    for (size_t i = 0; i < numCover.size(); i++) {
        if (verbose)
            cout << "Search " << i << " is used " << numCover[i] << " times\n";

        if (numCover[i] == 0)
            cout << "Warning: search " << scheme[i]
                 << " covers no error distributions!\n";
    }

    return ret;
}

// ----------------------------------------------------------------------------
// PARTITIONING
// ----------------------------------------------------------------------------
void RSearchStrategy::partition(const string& pattern, vector<Substring>& parts,
                               const int& numParts, const int& maxScore,
                               vector<BRSample>& exactMatchSamples,
                               Counters& counters) const {
    assert(exactMatchSamples.size() == (uint32_t)numParts);
    parts.clear();

    if (numParts >= (int)pattern.size() || numParts == 1) {
        // no need of splitting up since all parts would be one
        // character or less or there is only one part
        return;
    }
    parts.reserve(numParts);
    (this->*partitionPtr)(pattern, parts, numParts, maxScore, exactMatchSamples,
                          counters);
}

// Uniform Partitioning

void RSearchStrategy::partitionUniform(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<BRSample>& exactMatchSamples,
                                      Counters& counters) const {

    for (int i = 0; i < numParts; i++) {
        parts.emplace_back(pattern, (i * 1.0 / numParts) * pattern.size(),
                           ((i + 1) * 1.0 / numParts) * pattern.size());
    }
    // set end of final part correct
    parts.back().setEnd(pattern.size());

    // match the exactSample for each part
    index.setDirection(FORWARD);

    // the size of kmers in the kmer hash table
    length_t wordSize = index.getWordSize();

    for (int i = 0; i < numParts; i++) {
        const auto& current = parts[i];
        BRSample initialSample = index.getInitialSample();

        length_t start = 0;
        if (current.size() >= wordSize) {
            // skip wordsize
            initialSample = index.lookUpInKmerTable(Substring(
                current, current.begin() + 0, current.begin() + wordSize));
            start = wordSize;
        }
        Substring remainingPart(current, current.begin() + start,
                                current.end());
        exactMatchSamples[i] = index.matchStringBidirectionally(
            remainingPart, initialSample, counters);
    }
}

// Static Partitioning
void RSearchStrategy::partitionOptimalStatic(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<BRSample>& exactMatchSamples,
    Counters& counters) const {

    setParts(pattern, parts, numParts, maxScore);

    // match the exactSample for each part
    index.setDirection(FORWARD);

    // the size of kmers in the kmer hash table
    length_t wordSize = index.getWordSize();

    for (int i = 0; i < numParts; i++) {
        const auto& current = parts[i];
        BRSample initialSample = index.getInitialSample();

        length_t start = 0;
        if (current.size() >= wordSize) {
            // skip wordsize
            initialSample = index.lookUpInKmerTable(Substring(
                current, current.begin() + 0, current.begin() + wordSize));
            start = wordSize;
        }
        Substring remainingPart(current, current.begin() + start,
                                current.end());
        exactMatchSamples[i] = index.matchStringBidirectionally(
            remainingPart, initialSample, counters);
    }
}

void RSearchStrategy::setParts(const string& pattern, vector<Substring>& parts,
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

void RSearchStrategy::partitionDynamic(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<BRSample>& exactMatchSamples,
                                      Counters& counters) const {

    int matchedChars =
        seed(pattern, parts, numParts, maxScore, exactMatchSamples);
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
            if (exactMatchSamples[i].width() * weights[i] > maxRangeWeighted) {
                maxRangeWeighted = exactMatchSamples[i].width() * weights[i];
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
                    dir = (exactMatchSamples[i - 1].width() <
                           exactMatchSamples[i + 1].width())
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

        // match the new character
        index.setDirection(dir);
        index.addChar(c, exactMatchSamples.at(partToExtend), counters);
    }
}

int RSearchStrategy::seed(const string& pattern, vector<Substring>& parts,
                         const int& numParts, const int& maxScore,
                         vector<BRSample>& exactMatchSamples) const {
    int pSize = pattern.size();
    bool useKmerTable = (pSize >= 100);
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

    exactMatchSamples.resize(numParts);
    for (int i = 0; i < numParts; i++) {
        exactMatchSamples[i] = (useKmerTable)
                                ? index.lookUpInKmerTable(parts[i])
                                : index.getSampleOfSingleChar(parts[i][0]);
    }
    return numParts * wSize;
}

void RSearchStrategy::extendParts(const string& pattern,
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
#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
size_t RSearchStrategy::matchApprox(const string& pattern,
                                            length_t maxED, Counters& counters,
                                            const string& ID,
                                            const string& qual, bool revCompl) {
#else
vector<TextOcc> RSearchStrategy::matchApprox(const string& pattern,
                                            length_t maxED, Counters& counters,
                                            const string& ID,
                                            const string& qual, bool revCompl, bool report) {
#endif

    if (maxED == 0) {
        index.setDirection(FORWARD);
        auto m = index.exactMatchesOutput(pattern, counters);
#ifndef ENABLE_BENCHMARK_FUNCTIONALITY
        if (report) {
            generateSAM(m, pattern, ID, qual, revCompl);
        }
#endif
        return m;
    }

    // create the parts of the pattern
    vector<Substring> parts;

    // calculate how many parts there will be
    uint numParts = calculateNumParts(maxED);

    // the samples corresponding to the exact match for each part
    vector<BRSample> exactMatchSamples(numParts);

    // partition the read
    partition(pattern, parts, numParts, maxED, exactMatchSamples, counters);

    if (parts.empty() || numParts * maxED >= pattern.size()) {
        // splitting up was not viable -> just search the entire pattern
        cerr << "Warning: Normal bidirectional search was used as "
                "entered pattern is too short "
             << pattern.size() << endl;

        auto m = index.approxMatchesNaive(pattern, maxED, counters);
#ifndef ENABLE_BENCHMARK_FUNCTIONALITY
        if (report) {
            generateSAM(m, pattern, ID, qual, revCompl);
        }
#endif
        return m;
    }

    // The occurrences in the text and index
    BROccurrences occ;

    // END of preprocessing
    // A) do the searches for which the first part occurs more than the switch point

    index.reserveStacks(numParts,
                        pattern.length()); // reserve stacks for each part

    // create the bit-parallel alignment matrices
    index.resetMatrices(parts.size()); // reset the alignment matrix that will
                                       // be (possibly) used for each part

    // create the searches
    const vector<Search>& searches = createSearches(maxED, exactMatchSamples);

    for (const Search& s : searches) {
        doRecSearch(s, parts, occ, exactMatchSamples, counters);
    }

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
    // Iterate over all occurrences in  occ.second and accumulate the range width
    size_t count = 0;
    for (const auto& occ : occ.getIndexOccurrences()) {
            count += occ.getWidth();
    }
    return count;
#else
    // B) get all matches mapped to the text
    auto matches = (this->*filterPtr)(occ, maxED, counters);
    // C) generate sam output
    if (report) {
        generateSAM(matches, pattern, ID, qual, revCompl);
    }

    return matches;
#endif
}

void RSearchStrategy::doRecSearch(const Search& s, vector<Substring>& parts,
                                 BROccurrences& occ,
                                 const vector<BRSample>& exactMatchSamples,
                                 Counters& counters) const {

    if (s.getUpperBound(0) > 0) {
        // first part is allowed an error so start with an empty match
        s.setDirectionsInParts(parts);

        BRSample startSample = index.getInitialSample();
        BROcc startMatch = BROcc(startSample, 0, 0);
        (this->*startIdxPtr)(s, startMatch, occ, parts, counters, 0);
        return;
    }

    // first get the bidirectional match of first part
    int first = s.getPart(0);
    BRSample startSample = exactMatchSamples[first];
    if (!startSample.isValid()) {
        return;
    }

    // prepare the parts for this search
    s.setDirectionsInParts(parts);

    // can we continue exact matching according to this search?
    uint partInSearch = 1;
    length_t exactLength = parts[first].size();

    while (s.getUpperBound(partInSearch) == 0) {
        // extend the exact match
        index.setDirection(s.getDirection(partInSearch - 1));
        const auto& part = parts[s.getPart(partInSearch)];

        startSample =
            index.matchStringBidirectionally(part, startSample, counters);
        if (!startSample.isValid()) {
            return;
        }

        exactLength += part.size();
        partInSearch++;
    }

    // Create a match corresponding to the exact match
    BROcc startMatch = BROcc(startSample, 0, exactLength);
    // Start approximate matching in the index
    (this->*startIdxPtr)(s, startMatch, occ, parts, counters, partInSearch);
}

// ============================================================================
// CLASS CUSTOMSEARCHSTRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTION
// ----------------------------------------------------------------------------

void CustomRSearchStrategy::getSearchSchemeFromFolder(string pathToFolder,
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

    // get the info per distance score (scores between 1 and 4 are looked
    // at)
    for (int i = 1; i <= MAX_K; i++) {
        string name = pathToFolder + to_string(i) + "/searches.txt";
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
        ifstream stream_static(pathToFolder + to_string(i) +
                               "/static_partitioning.txt");

        if (stream_static) {
            // a file with static partitioning positions exists
            getline(stream_static, line);
            vector<string> postionsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                postionsAsString.push_back(token);
            }

            if (postionsAsString.size() != calculateNumParts(i) - 1) {
                throw runtime_error(
                    "Not enough static positions provided in " + pathToFolder +
                    to_string(i) + "/static_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " parts\nProvided: " +
                    to_string(postionsAsString.size()) + " parts");
            }

            for (auto str : postionsAsString) {
                staticPositions[i - 1].push_back(stod(str));
            }

            // check if these positions are valid
            sanityCheckStaticPartitioning(i);
            // if valid set getBeginsPointer to custom
            beginsPointer[i - 1] = &CustomRSearchStrategy::getBeginsCustom;

            stream_static.close();
        }
    }

    // get dynamic seeds and weights (if file exists)
    for (int i = 1; i <= MAX_K; i++) {
        ifstream stream_dynamic(pathToFolder + to_string(i) +
                                "/dynamic_partitioning.txt");

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
                    to_string(i) + "/dynamic_partitioning.txt\nExpected: " +
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
                &CustomRSearchStrategy::getSeedingPositionsCustom;
            weightsPointers[i - 1] = &CustomRSearchStrategy::getWeightsCustom;
        }
    }
}

Search CustomRSearchStrategy::makeSearch(const string& line,
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

void CustomRSearchStrategy::getVector(const string& vectorString,
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

void CustomRSearchStrategy::sanityCheckStaticPartitioning(
    const int& maxScore) const {
    const auto& positions = staticPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than 0
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

void CustomRSearchStrategy::sanityCheckDynamicPartitioning(
    const int& maxScore) const {

    const auto& seeds = seedingPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than 0
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

void CustomRSearchStrategy::sanityCheck(bool verbose) const {

    // check if for each supported edit distance all error distributions  are
    // covered

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
            throw runtime_error("Not all searches for distance " +
                                to_string(K) +
                                " have the same number of parts");
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
