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
#ifndef BUILDRINDEX_HPP
#define BUILDRINDEX_HPP

#include <cstdint>

#include "buildAuxiliary.h"
#include "plcp.h"
#include "buildMove.h"


/**
 * Print usage to stdout
*/
void showUsage();

/**
 * Parse command line arguments
 * @param argc number of arguments
 * @param argv char pointer array of arguments
 * @param [out] baseFN base filename
 * @param blockSize block size
 * @param verbose verbose ouput
 * @return True if succeeded, false otherwise
*/
bool parseArguments(int argc, char* argv[], string& baseFN, length_t& blockSize, bool& verbose);


void writeIntVectorBinary(const std::string& filename, const std::vector<length_t>& array);


/**
 * Builds samplesFirst and samplesLast.
 * @param samplesFirst Vector to fill. First sa samples of each input interval.
 * @param samplesLast Vector to fill. Last sa samples of each input interval.
 * @param SA The suffix array.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
*/
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast, vector<length_t>& SA, const string& BWT);


/**
 * Generate predecessor structures
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 * @param verbose print verbose output to stdout
*/
void generatePredecessors(const vector<length_t>& samplesFirst, 
    const vector<length_t>& samplesLast, SparseBitvec& predFirst, 
    SparseBitvec& predLast, length_t textLength);

/**
 * Generate predToRun arrays
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] firstToRun mapping between rank of ones in predFirst bitvector and run indices
 * @param [out] lastToRun mapping between rank of ones in predLast bitvector and run indices
*/
void generatePredTorun(const vector<length_t>& samplesFirst, const vector<length_t>& samplesLast,
    vector<length_t>& firstToRun, vector<length_t>& lastToRun);

/**
 * Create all datastructures necessary for the R index 
 * and write them to their respective files
 * @param baseFN base filename, all files will be generated as <baseFN>.<ext>
 * @param blockSize block size
 * @param verbose print verbose output to stdout
*/
void createRIndex(const string& baseFN, length_t blockSize, bool verbose);

int main(int argc, char* argv[]);


#endif

