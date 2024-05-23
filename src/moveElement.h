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
#ifndef MOVEHELPERS_HPP
#define MOVEHELPERS_HPP

#include <string>
#include <vector>
#include <iostream>
#include <cstdint>

#include "wordlength.h"
#include "alphabet.h"
#include "rindexhelpers.h"

#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

class MoveElement {

    public:

        // Run head for this move row. 'c' in the paper.
        uint8_t c;

        // The position index of the start of the input interval in the BWT.
        length_t inputStart;

        // The position index of the mapping of the first element of the input interval (first element of ouput interval).
        length_t outputStart;

        // Index of the input interval containing the mapping.
        length_t mappingIndex;

        MoveElement() {}

        MoveElement(uint8_t c, length_t inputStart, length_t outputStart, length_t mappingIndex): 
            c(c), inputStart(inputStart), outputStart(outputStart), mappingIndex(mappingIndex) {}

        void load(ifstream& ifs) {
            ifs.read((char*)&c, sizeof(c));
            ifs.read((char*)&inputStart, sizeof(inputStart));
            ifs.read((char*)&outputStart, sizeof(outputStart));
            ifs.read((char*)&mappingIndex, sizeof(mappingIndex));
        }

        length_t serialize(ofstream& out) const {
            out.write((char*)&c, sizeof(c));
            out.write((char*)&inputStart, sizeof(inputStart));
            out.write((char*)&outputStart, sizeof(outputStart));
            out.write((char*)&mappingIndex, sizeof(mappingIndex));

            return sizeof(c) + sizeof(inputStart) + sizeof(outputStart) + sizeof(mappingIndex);
        }

};

#endif