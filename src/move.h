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
#ifndef MOVE_HPP
#define MOVE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

#include "alphabet.h"
#include "util.h"
#include "wordlength.h"
#include "moveElement.h"
#include "rindexhelpers.h"
#include "sparseBitvec.h"

using namespace std;

#ifdef ENABLE_BENCHMARK_FUNCTIONALITY
extern length_t LF_call_count;
extern length_t elapsed_LF;
#endif


template <size_t S>
class Move {

    private:

        // Total size of the bwt.
        length_t bwtSize;

        // Total amount of runs.
        length_t runsSize;

        vector<MoveElement> mappings;

        // The run and position index of the zero character
        length_t zeroCharRun;
        length_t zeroCharPos;
 
    public:

        Move() {};

        /**
         * Constructor for the Move structure.
         * @param baseFN the base filename to write the move structure to.
         * @param verbose if True, print extra information during loading the structure.
        */
        Move(const string& baseFN, const bool verbose) {
            load(baseFN, verbose);
        }

        /**
         * Get the run head of the run with given index.
         * @param pos The position for which to get the run head within the block.
         * @return The run head of the given run index.
        */
        uint8_t getRunHead(const Position& pos) const {
            return mappings[pos.runIndex].c;
        }

        /**
         * Check if the run head of a specific run equals the given character.
         * @param runIndex The index of the run for which to check.
         * @param c The character for which to check.
         * @return True if the run head of the given run is c.
         */
        bool runHeadEquals(const length_t runIndex, const length_t c) const {
            // return runHeads.charAtEquals(runIndex, c);
            return mappings[runIndex].c == c;
        }

        /**
         * Perform left extension: match character c in the current range and perform LF operation on the boundaries to construct the new SA range.
         * @param inputRange The previous range in the SA (which match a string x), the range for which we match c.
         * @param c The character to match.
         * @param onlyPosIndex Only compute the position index after lf operation if true.
         * @return The new range in the SA, which contains matches of cx.
        */
        void leftExtension(const PositionRange& inputRange, PositionRange& newRange, length_t c, const bool onlyPosIndex) const;

        /**
         * Perform lf operation on a position.
         * @param input The position to map.
         * @param onlyPosIndex Only compute the position index after lf operation if true.
         * @return The position after lf mapping on the input.
        */
        void lf(const Position& input, Position& output, const bool onlyPosIndex) const;

        /**
         * Get the full BWT range (end inclusive).
         * @return Get the full range of bwt positions with corresponding run indexes.
        */
        PositionRange fullRange() const;

        /**
         * Get the run index of a position.
         * @param position The position of which to find the runIndex.
         * @param possibleRange The range of run indices in which the run index is located (end inclusive).
         * @return The run index of the position.
        */
        void getRunIndex(Position& position, pair<length_t,length_t> possibleRange) const;

        /**
         * Compute run indices of the positions of the given range.
         * @param PositionRange The range of which we will compute the run indices.
         * @return The range with the computed run indices.
        */
        PositionRange computeRunIndices(const PositionRange& positionRange) const;

        bool getNextOccurence(const PositionRange& startPosition,
                              Position& result, const length_t c) const;
        bool getPreviousOccurence(const PositionRange& startPosition,
                                  Position& result, const length_t c) const;

        /**
         * Print dPair and dIndex.
        */
        void printData();

        /**
         * Load the move data struction from an input stream.
         * @param in The input stream to load from.
        */
        void load(const string& baseFN, const bool verbose);


        /**
         * Get the amount of occurences of character c in the given range.
         * @param range The range to scan.
         * @param c count all occurences of this character
         * @returns The amount of occurences of c in the given range.
        */
        length_t getCharCountInRange(const PositionRange& range, const length_t c) const;

       /**
        * Get the amount of characters smaller than a given character c in the given range.
        * @param range The range to scan.
        * @param c count all smaller characters than c (exclusive).
        * @returns The amount of occurences of characters smaller than c in the given range.
       */
        length_t getSmallerCharsCountInRange(const PositionRange& range, const length_t c) const; 

      /**
       * print the allocated memory to stdout.
       * @param out The output stream to which to write the contents.
       * @return The size of the allocated memory.
      */
      size_t printMemSize(ofstream& out);
    
};

#endif