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

#ifndef BUILDHELPERS_H
#define BUILDHELPERS_H

#include "alphabet.h" // for Alphabet

#include <algorithm>  // for copy, max, transform, max_element, min_...
#include <array>      // for array
#include <cctype>     // for toupper, tolower
#include <exception>  // for exception
#include <fstream>    // for ifstream, ofstream
#include <functional> // Add this header for std::function
#include <iomanip>    // for operator<<, setprecision
#include <iostream>   // for operator<<, ifstream, ofstream, basic_o...
#include <random>     // for minstd_rand, uniform_int_distribution
#include <sdsl/int_vector.hpp>
#include <sstream>
#include <stdexcept> // for runtime_error, invalid_argument
#include <stdlib.h>  // for size_t, exit, EXIT_FAILURE, EXIT_SUCCESS
#include <string>    // for string, operator+, basic_string, allocator
#include <time.h>    // for clock, clock_t
#include <vector>    // for vector

using namespace std;

// Function to replace non-ACGT characters with a random ACGT character
char replaceNonACGT(char original, std::minstd_rand& gen,
                    length_t& replacementCounter);

/**
 * @brief Concatenates and transforms the sequences from a FASTA file.
 *
 * This function reads a FASTA file and concatenates the sequences into a single
 * string. It also replaces non-ACGT characters with a random ACGT character.
 * The start positions of each sequence and the sequence names are stored in
 * separate vectors. The resulting concatenated string is terminated with a '$'
 * character.
 *
 * @param fastaFile The path to the FASTA file.
 * @param concatenation The resulting concatenated string. (output)
 * @param positions The vector to store the start positions of each sequence.
 * (output)
 * @param seqNames The vector to store the sequence names. (output)
 * @param replaceFunc The function to use for replacing non-ACGT characters.
 * @param expectedNumber The expected number of sequences in the FASTA file.
 * Default is 12000.
 */
void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             length_t expectedNumber = 12000);

/**
 * @brief Read the contents of a text file into a string buffer.
 *
 */
void readText(const string& filename, string& buf);

void readSATextMode(const string& filename, vector<length_t>& sa,
                    size_t saSizeHint);

void readSA(const string& filename, vector<length_t>& sa, size_t saSizeHint);

void sanityCheck(const string& T, vector<length_t>& sa);

void createAndWriteHeaderInfo(const string& baseFN,
                              const vector<string>& seqNames,
                              const std::vector<length_t>& positions);

void writePositionsAndSequenceNames(const string& baseFN,
                                    const vector<length_t>& positions,
                                    const vector<string>& seqNames);

void checkTextSize(const size_t tSize);

void countChars(const string& T, vector<length_t>& charCounts);

void writeCharCounts(const string& baseFN, const vector<length_t>& charCounts);

void createAlphabet(const string& T, const vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma);

void createSuffixArray(const string& T, vector<length_t>& SA);

void createRevBWT(const string& baseFN, const string& T,
                  const vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, string& rBWT);

void preprocessFastaFile(const std::string& fastaFile,
                         const std::string& baseFN, string& T,
                         bool noWriting = false);

#endif