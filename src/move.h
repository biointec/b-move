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
#ifndef MOVELFREPR_H
#define MOVELFREPR_H

#include "definitions.h"  // for length_t
#include "rindexhelpers.h" // for SARange
#include "moveElement.h"         // for MoveRow
#include <memory>            // for allocator_traits<>::value_type
#include <stdint.h>          // for uint8_t
#include <string>            // for string
#include <utility>           // for pair
#include <vector>            // for vector

class MoveLFRepr {

  private:
    // Total size of the bwt.
    length_t bwtSize;

    // Total amount of runs.
    length_t nrOfRuns;

    // Rows of the move table
    std::vector<MoveRowLF> rows;

    // The position index of the zero character
    length_t zeroCharPos;

  public:
    /**
     * @brief Get the total amount of runs.
     *
     * @return length_t
     */
    length_t getNrOfRuns() const {
        return nrOfRuns;
    }

    /**
     * @brief Get the total size of the bwt.
     *
     * @return length_t
     */
    length_t getBwtSize() const {
        return bwtSize;
    }

    /**
     * Get the run head of the run with given index.
     * @param pos The position for which to get the run head within the block.
     * @return The run head of the given run index.
     */
    uint8_t getRunHead(const length_t& runIndex) const {
        return rows[runIndex].c;
    }

    /**
     * @brief Get the Input Start Pos object
     *
     * @param runIndex
     * @return length_t
     */
    length_t getInputStartPos(const length_t& runIndex) const {
        return rows[runIndex].inputStartPos;
    }

    /**
     * Find the run index for a given position.
     * @param position The position for which to find the run index.
     * @param runIndex The run index to be updated.
     * @param possibleRange The possible range for binary search.
     *
     * @return void
     */
    void getRunIndex(const length_t position, length_t& runIndex,
                     std::pair<length_t, length_t>& possibleRange) const;

    /**
     * Compute run indices of the positions of the given range.
     * @param range The range of which we will compute the run indices.
     *
     * @return void
     */
    void computeRunIndices(MoveRange& range) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index nextPos and run index
     * nextRun to indicate the next position in the SA interval corresponding to
     * an occurrence of c in the BWT (starting from the beginning of the
     * interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param nextPos The position in the SA interval that will contain the next
     * occurrence of c
     * @param nextRun The run index that will contain the next occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    bool walkToNextRun(const MoveRange& startRange, length_t& nextPos,
                       length_t& nextRun, const length_t c) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index previousPos and run index
     * previousRun to indicate the previous position in the SA interval
     * corresponding to an occurrence of c in the BWT (starting from the end of
     * the interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param previousPos The position in the SA interval that will contain the
     * previous occurrence of c
     * @param previousRun The run index that will contain the previous
     * occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    void walkToPreviousRun(const MoveRange& startRange, length_t& previousPos,
                           length_t& previousRun, const length_t c) const;

    /**
     * @brief Fast forward the runIndex until it contains the run that contains
     * the positionIndex.
     *
     * @param positionIndex The position index to fast forward to.
     * @param runIndex The run index to fast forward.
     */
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    /**
     * @brief Perform the LF operation on the given positionIndex and runIndex.
     *
     * @param positionIndex The position index to perform the LF operation on.
     * @param runIndex The run index to perform the LF operation on.
     */
    void findLF(length_t& positionIndex, length_t& runIndex) const;

    /**
     * @brief Perform the LF operation on the given positionIndex and runIndex
     * without fast forwarding.
     *
     * @param positionIndex The position index to perform the LF operation on.
     * @param runIndex The run index to perform the LF operation on.
     */
    void findLFWithoutFastForward(length_t& positionIndex,
                                  const length_t& runIndex) const;

    /**
     * @brief Extend the match corresponding to parentRange by prepending
     * character c to the match.
     *
     * @param parentRange The match to extend.
     * @param childRange The extended match.
     * @param c The character to prepend to the match.
     */
    void addChar(const SARange& parentRange, SARange& childRange,
                 const length_t& c) const;

    /**
     * @brief Extend the match corresponding to parentRange by appending
     * character c to the match. Only the count is necessary so no fast
     * forwarding.
     *
     * @param parentRange The match to extend.
     * @param c The character to append to the match.
     * @return length_t - The count of the character c in the extended match.
     */
    length_t countChar(const SARange& parentRange, const length_t& c) const;

    /**
     * @brief Extend the match corresponding to range by appending all
     * characters lexicographically smaller than c to the match. Only the count
     * is necessary so no fast forwarding.
     *
     * @param range The match to extend.
     * @param c The character to append to the match.
     */
    length_t getCumulativeCounts(const SARange& range,
                                 length_t positionInAlphabet) const;

    /**
     * @brief Load the move representation from a file.
     *
     * @param baseFile Base file name
     * @return true if the move representation was loaded successfully
     * @return false otherwise
     */
    bool load(const string& baseFile);
};

class MoveLFReprBP {
  public:
    // Constructor
    MoveLFReprBP(); // Default constructor

    // Initialize function to set up bit packing
    bool initialize(length_t nrOfRuns, length_t bwtSize);

    // Destructor
    ~MoveLFReprBP();

    // Load the move representation from a file
    bool load(const std::string& fileName);

    // Write the move representation to a file
    bool write(const std::string& fileName) const;

    // Get the number of runs
    length_t getNrOfRuns() const {
        return nrOfRuns;
    }

    // Get the total size of the BWT
    length_t getBwtSize() const {
        return textSize;
    }

    // Initialize values for each row
    void setRowValues(length_t rowIndex, uint8_t runChar,
                      length_t inputStartPos, length_t outputStartPos,
                      length_t outputStartRun);

    // Get the run head (character) for a given row
    uint8_t getRunHead(length_t runIndex) const;

    // Get the inputStartPos for a given row
    length_t getInputStartPos(length_t runIndex) const;

    // Get the outputStartPos for a given row
    length_t getOutputStartPos(length_t runIndex) const;

    // Get the outputStartRun for a given row
    length_t getOutputStartRun(length_t runIndex) const;

    // Set the outputStartRun value for a specific row
    void setOutputStartRun(length_t rowIndex, length_t value);

    // Set the zeroCharPos value
    void setZeroCharPos(length_t value) {
        zeroCharPos = value;
    }

    // Perform the LF-mapping step for positionIndex, updating runIndex
    void findLF(length_t& positionIndex, length_t& runIndex) const;

    /**
     * @brief Perform the LF operation on the given positionIndex and runIndex
     * without fast forwarding.
     *
     * @param positionIndex The position index to perform the LF operation on.
     * @param runIndex The run index to perform the LF operation on.
     */
    void findLFWithoutFastForward(length_t& positionIndex,
                                  const length_t& runIndex) const;

    // Fast-forward to the run containing a given position index
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    /**
     * Find the run index for a given position.
     * @param position The position for which to find the run index.
     * @param runIndex The run index to be updated.
     * @param possibleRange The possible range for binary search.
     *
     * @return void
     */
    void getRunIndex(const length_t position, length_t& runIndex,
                     std::pair<length_t, length_t>& possibleRange) const;

    /**
     * Compute run indices of the positions of the given range.
     * @param range The range of which we will compute the run indices.
     *
     * @return void
     */
    void computeRunIndices(MoveRange& range) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index nextPos and run index
     * nextRun to indicate the next position in the SA interval corresponding to
     * an occurrence of c in the BWT (starting from the beginning of the
     * interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param nextPos The position in the SA interval that will contain the next
     * occurrence of c
     * @param nextRun The run index that will contain the next occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    bool walkToNextRun(const MoveRange& startRange, length_t& nextPos,
                       length_t& nextRun, const length_t c) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index previousPos and run index
     * previousRun to indicate the previous position in the SA interval
     * corresponding to an occurrence of c in the BWT (starting from the end of
     * the interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param previousPos The position in the SA interval that will contain the
     * previous occurrence of c
     * @param previousRun The run index that will contain the previous
     * occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    void walkToPreviousRun(const MoveRange& startRange, length_t& previousPos,
                           length_t& previousRun, const length_t c) const;

    /**
     * @brief Extend the match corresponding to parentRange by prepending
     * character c to the match.
     *
     * @param parentRange The match to extend.
     * @param childRange The extended match.
     * @param c The character to prepend to the match.
     */
    void addChar(const SARange& parentRange, SARange& childRange,
                 const length_t& c) const;

    /**
     * @brief Extend the match corresponding to parentRange by appending
     * character c to the match. Only the count is necessary so no fast
     * forwarding.
     *
     * @param parentRange The match to extend.
     * @param c The character to append to the match.
     * @return length_t - The count of the character c in the extended match.
     */
    length_t countChar(const SARange& parentRange, const length_t& c) const;

    /**
     * @brief Extend the match corresponding to range by appending all
     * characters lexicographically smaller than c to the match. Only the count
     * is necessary so no fast forwarding.
     *
     * @param range The match to extend.
     * @param c The character to append to the match.
     */
    length_t getCumulativeCounts(const SARange& range,
                                 length_t positionInAlphabet) const;

  private:
    // The effective values for a given row
    uint8_t* buffer; // Buffer for bit-packed data

    // Number of runs and total size of the BWT
    length_t nrOfRuns; // Number of runs
    length_t textSize; // Total size of the BWT

    // Bit counts for each field in the struct
    uint8_t bitsForC; // Number of bits for the 'c' field
    uint8_t bitsForN; // Number of bits for the position fields
    uint8_t bitsForR; // Number of bits for the run field

    // Total bits and bytes for one packed row
    uint16_t totalBits;  // Total bits for one packed row
    uint16_t totalBytes; // Total bytes for one packed row

    // The position index of the zero character
    length_t zeroCharPos;

    // Helper to extract a value from the bit-packed buffer
    length_t getRowValue(length_t rowIndex, uint16_t bitOffset,
                         uint8_t numBits) const;

    // Helper to set a value in the bit-packed buffer
    void setRowValue(length_t rowIndex, length_t value, uint16_t bitOffset,
                     uint8_t numBits);
};

#endif