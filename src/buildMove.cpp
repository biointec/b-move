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
#include "buildMove.h"

using namespace std;


void fillMappings(map<length_t,pair<uint8_t, length_t>>& tIn, vector<MoveElement>& mappings, length_t size) {
    mappings.resize(size);
    
    map<length_t,pair<uint8_t, length_t>>::iterator tInIt = tIn.begin();
    length_t i = 0;
    while (tInIt != tIn.end()) {
        length_t inputStart = tInIt->first;
        uint8_t c = tInIt->second.first;
        length_t outputStart = tInIt->second.second;

        mappings[i] = MoveElement(c, inputStart, outputStart, 0);

        i ++;
        tInIt ++;
    }

    tInIt = tIn.begin();
    i = 0;
    while (tInIt != tIn.end()) {
        length_t outputIndex = getRunIndex(tInIt->second.second, size, mappings);

        mappings[i].mappingIndex = outputIndex;

        i ++;
        tInIt ++;
    }
}


length_t getRunIndex(length_t position, length_t size, vector<MoveElement>& mappings) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the posible range smaller by binary search, untill only 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the test value.
        if (mappings[testIndex].inputStart <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex-1;
        }
    }

    assert(position >= mappings[leftBoundary].inputStart);
    assert(leftBoundary == size-1 || position < mappings[leftBoundary+1].inputStart);

    return leftBoundary;
}
