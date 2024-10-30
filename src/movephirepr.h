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
#ifndef MOVEPHIREPR_H
#define MOVEPHIREPR_H

#include "definitions.h"  // for length_t
#include "rindexhelpers.h" // for SARange
#include "moveElement.h"         // for MoveRow

class MovePhiRepr {

  private:
    // Total size of the text.
    length_t textSize;

    // Total amount of runs.
    length_t nrOfRuns;

    // Rows of the phi move table
    std::vector<MoveRowPhi> phiRows;

  public:
    /**
     * @brief Load the move representation from a file.
     *
     * @param baseFile Base file name
     * @param verbose Print verbose output
     * @return true if the move representation was loaded successfully
     * @return false otherwise
     */
    bool load(const string& baseFile);

    /**
     * @brief Fast forward the runIndex until it contains the run that contains
     * the positionIndex.
     *
     * @param positionIndex The position index
     * @param runIndex The run index (will be updated)
     */
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    /**
     * @brief Perform the phi operation on the given position index and run
     * index.
     *
     * @param positionIndex The position index (will be updated)
     * @param runIndex The run index (will be updated)
     */
    void phi(length_t& positionIndex, length_t& runIndex) const;
};

class MovePhiReprBP {
  public:
    MovePhiReprBP(); // Default constructor stub

    // Initialize function to handle memory allocation
    bool initialize(length_t nrOfRuns, length_t textSize);

    // Load function to handle initialization
    bool load(const std::string& fileName);

    // Write function to handle serialization
    bool write(const std::string& fileName) const;

    // Destructor
    ~MovePhiReprBP();

    // return the number of runs
    length_t size() const {
        return nrOfRuns;
    }

    // Initialize values for each row
    void setRowValues(length_t rowIndex, length_t inputStartPos,
                      length_t outputStartPos, length_t outputStartRun);

    // Get the inputStartPos value for a specific row
    length_t getInputStartPos(length_t i) const;

    // Get the outputStartPos value for a specific row
    length_t getOutputStartPos(length_t i) const;

    // Get the outputStartRun value for a specific row
    length_t getOutputStartRun(length_t i) const;

    // Set the outputStartRun value for a specific row
    void setOutputStartRun(length_t rowIndex, length_t value);

    // Fast-forward runIndex to the row containing positionIndex
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    // Perform the LF-mapping step for positionIndex, updating the runIndex
    void phi(length_t& positionIndex, length_t& runIndex) const;

  private:
    uint8_t* buffer;     // Buffer holding all packed row data
    length_t nrOfRuns;   // Number of runs
    length_t textSize;   // Size of the text for LF-mapping
    uint8_t bitsForN;    // Number of bits for 'n'
    uint8_t bitsForR;    // Number of bits for 'r'
    uint16_t totalBits;  // Total bits for one packed row
    uint16_t totalBytes; // Total bytes for one packed row

    // Helper to extract a specific value from the buffer
    length_t getRowValue(length_t rowIndex, uint16_t bitOffset,
                         uint8_t numBits) const;

    // Helper to set a specific value in the buffer
    void setRowValue(length_t rowIndex, length_t value, uint16_t bitOffset,
                     uint8_t numBits);
};

#endif