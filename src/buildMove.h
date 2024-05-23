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
#ifndef BUILDMOVE_HPP
#define BUILDMOVE_HPP

#include <string>
#include <vector>
#include <iostream>
#include <cstdint>

#include "alphabet.h"
#include "util.h"
#include "wordlength.h"
#include "move.h"
#include "moveElement.h"
#include "sparseBitvec.h"

using namespace std;

/**
 * Fill dIndex array.
 * @param tIn Balanced tree: stores pairs of I_out with increasing input interval index.
 * @param mappings The vector to fill, mappings[i] = input-output interval i + index j of input interval containing q_i.
 * @param size The amount of input-output intervals.
*/
void fillMappings(map<length_t,pair<uint8_t, length_t>>& tIn, vector<MoveElement>& mappings, length_t size);


/**
 * Create array of characters 'runChar', which contains the character for each input interval.
 * @param runChars The vector to create. Corresponding characters of the input intervals.
 * @param size The amount of input-output intervals.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 * @param sigma The alphabet.
 * @param dPair Vector with pairs of input-output interval start positions.
*/
template <size_t S>
void createRunChars(length_t size,
                    const string& BWT, const Alphabet<S>& sigma,
                    vector<MoveElement>& mappings, string& stringOfRunHeads) {
    // Allocate memory for the array.
    stringOfRunHeads.resize(size);

    // Store the character of each run.
    for (length_t i = 0; i < size; i++) {
        stringOfRunHeads[i] = BWT[mappings[i].inputStart];
    }
}

/**
 * Get run index of a given position.
 * @param position The position for which to find the run index.
 * @param size The amount of input-output intervals.
 * @param dPair Vector with pairs of input-output interval start positions.
 * @returns The run index of the given position.
*/
length_t getRunIndex(length_t position, length_t size, vector<MoveElement>& mappings);


/**
 * Creates the Move structure.
 * @param baseFN the base filename to write the move structure to.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 * @param charCounts Accumulated number of characters in lex order.
 * @param sigma The alphabet.
 * @param SA The suffix array
*/
template <size_t S>
void createMove(const string& baseFN, const string& BWT, vector<length_t>& charCounts, const Alphabet<S>& sigma, vector<length_t>& SA) {

    // Set bwt size
    length_t bwtSize = BWT.size();

    // balanced tree: stores pairs of I_out with increasing input interval index
    map<length_t, pair<uint8_t, length_t>> tIn;

    // Create accumulated charCounts
    length_t charCountsAcc[S];
    length_t total = 0;
    for (size_t i = 0; i < S; i ++) {
        char c = sigma.i2c(i);
        charCountsAcc[i] = total;
        total += charCounts[c];
    }
    
    // fill tIn, tOut
    length_t charsVisisted[S] = {0};
    length_t prevC = S;
    length_t zeroCharRun;
    length_t zeroCharPos;
    for (length_t i = 0; i < bwtSize; i ++) {

        length_t c = sigma.c2i(BWT[i]);
        if (c == 0) {
            zeroCharPos = i;
            zeroCharRun = tIn.size();
        }
        if (prevC != c) {
            length_t lf = charCountsAcc[c] + charsVisisted[c];
            tIn[i] = make_pair(c, lf);
        }

        charsVisisted[c] ++;
        prevC = c;
    }

    // Set dPair and dIndex size
    length_t arraySize = tIn.size();
    length_t size = arraySize;
    cout << "There are " << size << " input-output intervals." << endl;

    // Fill dPair and dIndex
    vector<MoveElement> mappings;
    fillMappings(tIn, mappings, size);

    // Write Move structures to files
    string fileName = baseFN;
    fileName += ".move";
    ofstream ofs(fileName, ios::binary);

    // Write bwtSize, amount of input-output intervals and the alphabet size to the output stream.
    ofs.write((char*)&bwtSize, sizeof(bwtSize));
    ofs.write((char*)&size, sizeof(size));
    ofs.write((char*)&zeroCharPos, sizeof(zeroCharPos));
    ofs.write((char*)&zeroCharRun, sizeof(zeroCharRun));

    // Write mappings to the output stream.
    for (length_t i = 0; i < size; i ++) {
        mappings[i].serialize(ofs);
    }

    ofs.close();
    cout << "Wrote file " << fileName << endl;

    cout << "Number of rows in Move structure: " << size << endl;
    
}

#endif