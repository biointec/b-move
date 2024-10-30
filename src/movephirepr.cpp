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

#include "movephirepr.h"
#include "definitions.h"
#include "rindexhelpers.h"
#include "logger.h"
#include "moveElement.h"

#include <cassert>
#include <istream>

using namespace std;

#ifdef PHI_BENCHMARK_FUNCTIONALITY
length_t phi_call_count = 0;
length_t locate_call_count = 0;
length_t elapsed_phi = 0;
chrono::duration<double> elapsed_locate{0.0};
#endif

// ----------------------------------------------------------------------------
//  MovePhiRepr
// ----------------------------------------------------------------------------

bool MovePhiRepr::load(const string& fileName) {
    ifstream ifs(fileName, ios::binary);
    if (!ifs) {
        return false;
    }

    // Load the bwtSize, amount of input intervals and the alphabet size.
    ifs.read((char*)&textSize, sizeof(textSize));
    ifs.read((char*)&nrOfRuns, sizeof(nrOfRuns));

    // Load rows: allocate memory and read.
    phiRows.reserve(nrOfRuns + 1);
    for (length_t i = 0; i < nrOfRuns + 1; i++) {
        phiRows[i] = MoveRowPhi(ifs);
        assert(phiRows[i].outputStartRun <= nrOfRuns);
    }

    ifs.close();
    return true;
}

void MovePhiRepr::fastForward(const length_t& positionIndex,
                              length_t& runIndex) const {
    // Fast forward the runIndex until it contains
    // the run that contains the positionIndex.
    while (phiRows[runIndex].inputStartPos <= positionIndex) {
        runIndex++;
        // Ensure runIndex stays within bounds
        assert(runIndex < nrOfRuns + 1);
    }
    runIndex--;
}
void MovePhiRepr::phi(length_t& positionIndex, length_t& runIndex) const {
#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // Start timing with rdtscp
    unsigned int start_aux;
    uint64_t start = __builtin_ia32_rdtscp(&start_aux);
#endif
    assert(runIndex < nrOfRuns);
    assert(positionIndex < textSize);
    const auto& row = phiRows[runIndex];
    length_t offset = positionIndex - row.inputStartPos;
    positionIndex = row.outputStartPos + offset;
    runIndex = row.outputStartRun;
    // Fast forward to the correct runIndex after LF operation
    fastForward(positionIndex, runIndex);

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // End timing with rdtscp
    unsigned int end_aux;
    uint64_t end = __builtin_ia32_rdtscp(&end_aux);

    // Check if the core has switched
    if (start_aux != end_aux) {
        cerr << "Core switch detected in phi after " << phi_call_count
             << " calls" << endl;

    } else {

        elapsed_phi += end - start;
        phi_call_count++;
    }
#endif
}

// ----------------------------------------------------------------------------
//  MovePhiReprBP
// ----------------------------------------------------------------------------

// Stub default constructor (does nothing)
MovePhiReprBP::MovePhiReprBP()
    : nrOfRuns(0), textSize(0), bitsForN(0), bitsForR(0), totalBits(0),
      totalBytes(0) {
    // No initialization done here
}

bool MovePhiReprBP::initialize(length_t nrOfRuns, length_t textSize) {
    this->nrOfRuns = nrOfRuns;
    this->textSize = textSize;

    // Calculate the number of bits required for n and r
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));

    logger.logInfo("\tInitializing MovePhiReprBP with " +
                   std::to_string(nrOfRuns) + " runs and text size " +
                   std::to_string(textSize) + "...");

    // Calculate total bits for each row and convert to bytes
    totalBits = 2 * bitsForN + bitsForR; // Two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    logger.logInfo("\tBits for n: " + std::to_string(bitsForN));
    logger.logInfo("\tBits for r: " + std::to_string(bitsForR));
    logger.logInfo("\tTotal bits per row: " + std::to_string(totalBits));
    logger.logInfo("\tTotal bytes per row: " + std::to_string(totalBytes));

    // Manually allocate buffer memory
    buffer = new uint8_t[totalBytes * (nrOfRuns + 1)]();

    return true;
}

bool MovePhiReprBP::load(const std::string& fileName) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs) {
        return false; // File could not be opened
    }

    // Load the textSize and nrOfRuns from the file
    ifs.read(reinterpret_cast<char*>(&textSize), sizeof(textSize));
    ifs.read(reinterpret_cast<char*>(&nrOfRuns), sizeof(nrOfRuns));

    // Now calculate the number of bits required for n and r
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));

    // Calculate total bits for each row and convert to bytes
    totalBits = 2 * bitsForN + bitsForR; // Two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    // Manually allocate buffer memory
    buffer = new uint8_t[totalBytes * (nrOfRuns + 1)]();

    // Load the packed rows into the buffer
    ifs.read(reinterpret_cast<char*>(buffer), totalBytes * (nrOfRuns + 1));

    // Ensure the file was properly read
    if (!ifs) {
        return false; // Read failed
    }

    ifs.close();
    return true;
}

bool MovePhiReprBP::write(const std::string& fileName) const {
    std::ofstream ofs(fileName, std::ios::binary);
    if (!ofs) {
        return false; // File could not be opened for writing
    }

    // Write the textSize and nrOfRuns to the file
    ofs.write(reinterpret_cast<const char*>(&textSize), sizeof(textSize));
    ofs.write(reinterpret_cast<const char*>(&nrOfRuns), sizeof(nrOfRuns));

    // Write the buffer (packed rows) to the file
    ofs.write(reinterpret_cast<const char*>(buffer),
              totalBytes * (nrOfRuns + 1));

    // Ensure the file was properly written
    if (!ofs) {
        return false; // Write failed
    }

    ofs.close();
    return true;
}

MovePhiReprBP::~MovePhiReprBP() {
    delete[] buffer; // Deallocate buffer memory
}

void MovePhiReprBP::setRowValues(length_t rowIndex, length_t inputStartPos,
                                 length_t outputStartPos,
                                 length_t outputStartRun) {
    setRowValue(rowIndex, inputStartPos, 0,
                bitsForN); // Pack inputStartPos at bit offset 0
    setRowValue(rowIndex, outputStartPos, bitsForN,
                bitsForN); // Pack outputStartPos at bit offset bitsForN
    setRowValue(rowIndex, outputStartRun, 2 * bitsForN,
                bitsForR); // Pack outputStartRun at bit offset 2 * bitsForN
}

length_t MovePhiReprBP::getRowValue(length_t rowIndex, uint16_t bitOffset,
                                    uint8_t numBits) const {
    length_t mask = (1ULL << numBits) - 1;
    length_t byteIndex = rowIndex * totalBytes + (bitOffset / 8);
    uint8_t bitIndex = bitOffset % 8;

    const __uint128_t* block =
        reinterpret_cast<const __uint128_t*>(&buffer[byteIndex]);
    __uint128_t currentValue = *block;

    return (currentValue >> bitIndex) & mask;
}

length_t MovePhiReprBP::getInputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, 0,
                       bitsForN); // Use getRowValue to get the inputStartPos
}

length_t MovePhiReprBP::getOutputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, bitsForN,
                       bitsForN); // Use getRowValue to get the outputStartPos
}

length_t MovePhiReprBP::getOutputStartRun(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, 2 * bitsForN,
                       bitsForR); // Use getRowValue to get the outputStartRun
}

void MovePhiReprBP::setRowValue(length_t rowIndex, length_t value,
                                uint16_t bitOffset, uint8_t numBits) {
    length_t mask = (1ULL << numBits) - 1;
    value &= mask; // Ensure value fits in the specified number of bits

    length_t byteIndex = rowIndex * totalBytes + (bitOffset / 8);
    uint8_t bitIndex = bitOffset % 8;

    __uint128_t* block = reinterpret_cast<__uint128_t*>(&buffer[byteIndex]);
    __uint128_t currentValue = *block;

    __uint128_t clearMask = ~(__uint128_t(mask) << bitIndex);
    currentValue &= clearMask;
    currentValue |= (__uint128_t(value) << bitIndex);
    *block = currentValue;
}

void MovePhiReprBP::setOutputStartRun(length_t rowIndex, length_t value) {
    assert(rowIndex < nrOfRuns); // Ensure the index is within bounds
    setRowValue(rowIndex, value, 2 * bitsForN,
                bitsForR); // Pack value into the correct bit offset
}

void MovePhiReprBP::fastForward(const length_t& positionIndex,
                                length_t& runIndex) const {
    // Fast forward the runIndex until it contains
    // the run that contains the positionIndex.
    while (getInputStartPos(runIndex) <= positionIndex) {
        runIndex++;
        // Ensure runIndex stays within bounds
        assert(runIndex < nrOfRuns + 1);
    }
    runIndex--;
}

void MovePhiReprBP::phi(length_t& positionIndex, length_t& runIndex) const {
#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // Start timing with rdtscp
    unsigned int start_aux;
    uint64_t start = __builtin_ia32_rdtscp(&start_aux);
#endif
    assert(runIndex < nrOfRuns);
    assert(positionIndex < textSize);
    length_t offset = positionIndex - getInputStartPos(runIndex);
    positionIndex = getOutputStartPos(runIndex) + offset;
    runIndex = getOutputStartRun(runIndex);
    // Fast forward to the correct runIndex after LF operation
    fastForward(positionIndex, runIndex);

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // End timing with rdtscp
    unsigned int end_aux;
    uint64_t end = __builtin_ia32_rdtscp(&end_aux);

    // Check if the core has switched
    if (start_aux != end_aux) {
        cerr << "Core switch detected in phi after " << phi_call_count
             << " calls" << endl;

    } else {

        elapsed_phi += end - start;
        phi_call_count++;
    }
#endif
}
