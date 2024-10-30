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

#include "alphabet.h"    // for Alphabet
#include "definitions.h" // for length_tr

#include <functional> // for function
#include <random>     // for minstd_rand
#include <stdlib.h>   // for size_t
#include <string>     // for string
#include <vector>     // for vector

// Function to replace non-ACGT characters with a random ACGT character
char replaceNonACGT(char original, std::minstd_rand& gen,
                    const std::string& seed, size_t& seedIndex);

// Function to replace non-ACGT characters with a seeded ACGT character
char replaceNonACGTWithSeed(char original, std::minstd_rand& gen,
                            const std::string& seed, size_t& seedIndex);

/**
 * @brief Concatenates and transforms the sequences from a FASTA file.
 *
 * This function reads a FASTA file and concatenates the sequences into a single
 * std::string. It also replaces non-ACGT characters with a random ACGT
 * character. The start positions of each sequence and the sequence names are
 * stored in separate std::vectors. The resulting concatenated std::string is
 * terminated with a '$' character.
 *
 * @param fastaFile The path to the FASTA file.
 * @param concatenation The resulting concatenated std::string. (output)
 * @param positions The std::vector to store the start positions of each
 * sequence. (output)
 * @param seqNames The std::vector to store the sequence names. (output)
 * @param replaceFunc The function to use for replacing non-ACGT characters.
 * @param expectedNumber The expected number of sequences in the FASTA file.
 * Default is 12000.
 */
void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             std::function<char(char, size_t&)> replaceFunc,
                             length_t expectedNumber = 12000);

/**
 * @brief Read the contents of a text file into a std::string buffer.
 *
 */
void readText(const std::string& filename, std::string& buf);

void readSATextMode(const std::string& filename, std::vector<length_t>& sa,
                    size_t saSizeHint);

void readSA(const std::string& filename, std::vector<length_t>& sa,
            size_t saSizeHint);

void sanityCheck(const std::string& T, std::vector<length_t>& sa);

void createAndWriteHeaderInfo(const std::string& baseFN,
                              const std::vector<std::string>& seqNames,
                              const std::vector<length_t>& positions);

void writePositionsAndSequenceNames(const std::string& baseFN,
                                    const std::vector<length_t>& positions,
                                    const std::vector<std::string>& seqNames);

void checkTextSize(const size_t tSize);

void countChars(const std::string& T, std::vector<length_t>& charCounts);

void writeCharCounts(const std::string& baseFN,
                     const std::vector<length_t>& charCounts);

void createAlphabet(const std::string& T,
                    const std::vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma);

void createSuffixArray(const std::string& T, std::vector<length_t>& SA);

void createRevBWT(const std::string& baseFN, const std::string& T,
                  const std::vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, std::string& rBWT);

void preprocessFastaFile(const std::string& fastaFile,
                         const std::string& baseFN, std::string& T,
                         length_t seedLength, bool noWriting = false);

void writeTagAndCompiledInfo(const std::string& baseFN);

void generateBWT(const std::string& T, const std::vector<length_t>& SA,
                 std::string& BWT);

/**
 * @brief Determine the input mode for the index construction. The input can
 * either be a txt file with created suffix arrays or a fasta file.
 * @param baseFN The base file name of the input file.
 * @param textMode A boolean indicating whether the input is in text mode.
 * (output)
 * @param fastaMode A boolean indicating whether the input is in fasta mode.
 * (output)
 * @param fastaExtension The extension of the fasta file if it is in fastaMode.
 * (output)
 */
void determineInputMode(std::string& baseFN, bool& textMode, bool& fastaMode,
                        std::string& fastaExtension);

void writeCharCountsAndCreateAlphabet(const std::string& baseFN,
                                      const std::string& T,
                                      Alphabet<ALPHABET>& sigma,
                                      std::vector<length_t>& charCounts);

void readOrCreateSAWithSanityCheck(const std::string& baseFN,
                                   std::vector<length_t>& SA,
                                   const std::string& T, bool fromFasta);

void readOrCreateRevSAWithSanityCheck(const std::string& baseFN,
                                      std::vector<length_t>& SA,
                                      const std::string& T, bool fromFasta);

void preprocessTextMode(const std::string& baseFN, std::string& T);
#endif