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

#ifndef BUILDINDEX_H
#define BUILDINDEX_H

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include "wordlength.h"
#include "alphabet.h"

using namespace std;


// Function to replace non-ACGT characters with a random ACGT character
char replaceNonACGT(char original, std::minstd_rand& gen, length_t& replacementCounter) {
    const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        replacementCounter++;
        return validChars[gen() % validChars.length()];
    }
    return original;
}


void processSequence(std::ofstream& ofs,
                     std::stringstream& ss,
                     std::string& sequence,
                     length_t& minRemovalLength,
                     length_t& startPosition,
                     std::vector<length_t>& positions,
                     std::vector<std::string>& seqNames,
                     std::minstd_rand& gen,
                     length_t& totalNsRemoved,
                     length_t& totalNsReplaced,
                     length_t& totalOthersReplaced) {
    ofs << "Processing sequence: " << seqNames.back() << "\n";

    bool inNRun = false;
    length_t nRunStart = 0;

    length_t NsRemoved = 0;
    length_t NsReplaced = 0;
    length_t othersReplaced = 0;

    length_t charCount = 0;
    for (char& c : sequence) {
        if (std::toupper(c) == 'N') {
            if (!inNRun) {
                inNRun = true;
                nRunStart = charCount;
            }
        } else {
            if (inNRun) {
                if (charCount - nRunStart < minRemovalLength) {
                    // Replace N with a random ACGT character
                    for (size_t i = 0; i < charCount - nRunStart; ++i) {
                        ss << replaceNonACGT('N', gen, NsReplaced);
                    }
                    ofs << "Replaced " << charCount - nRunStart << " Ns with random ACGT characters at coordinate " << nRunStart << "\n";
                } else{
                    ofs << "Removed " << charCount - nRunStart << " Ns at coordinate " << nRunStart << "\n";
                    NsRemoved += charCount - nRunStart;
                }
                inNRun = false;
            }
            // Replace non-ACGT characters
            c = replaceNonACGT(std::toupper(c), gen, othersReplaced);
            ss << c;
        }
        charCount++;
    }

    if (inNRun){
        if (charCount - nRunStart < minRemovalLength) {
            // Replace N with a random ACGT character
            for (size_t i = 0; i < charCount - nRunStart; ++i) {
                ss << replaceNonACGT('N', gen, NsReplaced);
            }
            ofs << "Replaced " << charCount - nRunStart << " Ns with random ACGT characters at coordinate " << nRunStart << "\n";
        } else{
            ofs << "Removed " << charCount - nRunStart << " Ns at coordinate " << nRunStart << "\n";
            NsRemoved += charCount - nRunStart;
        }
    }

    positions.emplace_back(startPosition); // push back the start position of this sequence
    startPosition += sequence.length(); // create start position for next sequence
    sequence.clear();  // Clear the sequence for the next one

    ofs << "Removed " << NsRemoved << " Ns\n";
    ofs << "Replaced " << NsReplaced << " Ns with random ACGT characters\n";
    ofs << "Replaced " << othersReplaced << " non-ACGT characters with random ACGT characters\n";
    ofs << "\n";

    totalNsRemoved += NsRemoved;
    totalNsReplaced += NsReplaced;
    totalOthersReplaced += othersReplaced;
}


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
 * @param concatenation The resulting concatenated string.
 * @param positions The vector to store the start positions of each sequence.
 * @param seqNames The vector to store the sequence names.
 * @param expectedNumber The expected number of sequences in the FASTA file.
 * Default is 12000.
 */
void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames, 
                             std::string outputFilename,
                             length_t expectedNumber = 12000) {
    // Open input file
    std::ifstream inputFile(fastaFile);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error opening file: " + fastaFile);
    }

    // Open output file
    std::ofstream ofs(outputFilename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Error opening file: " + outputFilename);
    }

    // Temporary variables for processing
    std::string line;
    length_t startPosition = 0;

    // Seed for random number generation
    std::minstd_rand gen(42);

    // Reserve memory for positions and seqNames vectors
    positions.reserve(expectedNumber);
    seqNames.reserve(expectedNumber);

    // Stringstream to concatenate the sequences
    std::stringstream ss;
    std::string sequence;

    length_t totalNsRemoved = 0;
    length_t totalNsReplaced = 0;
    length_t totalOthersReplaced = 0;

    length_t minRemovalLength = 4;

    while (std::getline(inputFile, line)) {
        if (line.empty())
            continue; // Skip empty lines

        if (line[0] == '>') { // Sequence name line
            // Process the previous sequence
            if (!sequence.empty()) {
                processSequence(ofs, ss, sequence, minRemovalLength,
                                startPosition, positions, seqNames, gen, totalNsRemoved,
                                totalNsReplaced, totalOthersReplaced);
            }
            std::string description = line.substr(1);
            // get the sequence name (before the first space)
            seqNames.emplace_back(description.substr(
                0, description.find(' '))); // Store sequence name
        } else {
            // Append to the current sequence
            sequence += line;
        }
    }

    // Process the last sequence in the file
    if (!sequence.empty()) {
        processSequence(ofs, ss, sequence, minRemovalLength,
                        startPosition, positions, seqNames, gen, totalNsRemoved,
                        totalNsReplaced, totalOthersReplaced);
    }

    // Close input file
    inputFile.close();

    concatenation = ss.str() + "$";

    ofs << "Total removed Ns: " << totalNsRemoved << "\n";
    ofs << "Total replaced Ns: " << totalNsReplaced << "\n";
    ofs << "Total replaced non-ACGT characters: " << totalOthersReplaced << "\n";

    ofs.close();
}


/**
 * Read the text and store it in a string buffer
 * @param filename filename of the text
 * @param [out] buf output string buffer
*/
void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

/**
 * Read the suffix array in text mode
 * @param filename filename of the suffix array file
 * @param [out] sa suffix array
 * @param saSizeHint expected size of the suffix array
*/
void readSATextMode(const string& filename, vector<length_t>& sa,
                    size_t saSizeHint) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    sa.reserve(saSizeHint);
    length_t el;
    while (ifs >> el)
        sa.push_back(el);
}

/**
 * Read the suffix array in binary mode
 * @param filename filename of the suffix array file
 * @param [out] sa suffix array
 * @param saSizeHint expected size of the suffix array
*/
void readSA(const string& filename, vector<length_t>& sa, size_t saSizeHint) {
    ifstream ifs(filename, ios::binary);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    size_t numElements = ifs.tellg() / sizeof(length_t);

    if (numElements == saSizeHint) { // file is likely binary
        sa.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, ios::beg);
        ifs.read((char*)sa.data(), sa.size() * sizeof(length_t));
    } else { // try to read SA in text mode
        readSATextMode(filename, sa, saSizeHint);
    }
}

/**
 * Perform sanity check on the text and the suffix array
 * @param T the text
 * @param sa the suffix array
*/
void sanityCheck(const string& T, vector<length_t>& sa) {
    // check T for correctness
    if (T.back() == '\n')
        throw runtime_error("T should end with a \'$\' character, "
                            "not with a newline");

    if (T.back() != '$')
        throw runtime_error("T should end with a \'$\' character");

    if (sa.size() != T.size())
        throw runtime_error("Text and suffix array contain a "
                            "different number of elements");

    // briefly check the suffix array
    length_t min = *min_element(sa.begin(), sa.end());
    length_t max = *max_element(sa.begin(), sa.end());

    if (min == 1 && max == T.size()) { // rebase to [0..T.size()-1]
        for (auto& el : sa)
            el--;
        min--;
        max--;
    }

    if (min != 0 || max != T.size() - 1)
        throw runtime_error("Suffix array must contain numbers between "
                            "[0 and " +
                            to_string(T.size() - 1) + "]");

    // check if all numbers in the suffix array are present
    vector<bool> bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T
}

/**
 * Write an array (C++ vector) to a binary file
 * @param filename name of the file to write to
 * @param array array to be written
*/
void writeArrayBinary(const string& filename, const vector<length_t>& array) {
    ofstream ofs(filename, ios::binary);
    ofs.write((char*)array.data(), array.size() * sizeof(length_t));
    ofs.close();
}

/**
 * Read the text and generate the alphabet datastructure
 * @param baseFN base filename
 * @param [out] T the text
 * @param [out] charCounts array with length of 256 with number of occurrences of all ASCII characters in the text
 * @param [out] sigma the alphabet datastructure
*/
void generateText(const string& baseFN, string& T) {
    
    // read the text file from disk
    cout << "Reading " << baseFN << ".txt..." << endl;
    readText(baseFN + ".txt", T);
}


/**
 * Read the text and generate the alphabet datastructure
 * @param baseFN base filename
 * @param [out] T the text
 * @param [out] charCounts array with length of 256 with number of occurrences of all ASCII characters in the text
 * @param [out] sigma the alphabet datastructure
*/
void generateAlphabet(string& T, vector<length_t>& charCounts, Alphabet<ALPHABET>& sigma) {

    // count the frequency of each characters in T
    charCounts = vector<length_t>(256, 0);
    for (char c : T) {
        charCounts[(unsigned char)c]++;
    }

    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts) {
        if (count > 0) {
            nUniqueChar++;
        }
    }

    cout << "\tText has length " << T.size() << "\n";
    cout << "\tText has " << nUniqueChar << " unique characters\n";

    if (nUniqueChar > ALPHABET) {
        cerr << "FATAL ERROR: the number of unique characters in the "
             << "text exceeds the alphabet size. Please recompile"
             << "bmove using a higher value for ALPHABET " << endl;
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        cout << "WARNING: the number of unique characters in the "
             << "text is less than the ALPHABET size specified when "
             << "bmove was compiled. Performance may be affected\n";
    }

    sigma = Alphabet<ALPHABET>(charCounts);
}



/**
 * Generate the suffix array and the BWT string
 * @param baseFN base filename
 * @param T the text
 * @param [out] SA the suffix array
 * @param [out] BWT the BWT string
*/
void generateSA(const string& baseFN, const string& T, 
                 vector<length_t>& SA) {   

    // read the suffix array
    cout << "Reading " << baseFN << ".sa..." << endl;
    readSA(baseFN + ".sa", SA, T.size());

    // perform a sanity check on the suffix array
    cout << "\tPerforming sanity checks..." << endl;
    sanityCheck(T, SA);
    cout << "\tSanity checks OK" << endl;
}

/**
 * Generate the suffix array and the BWT string
 * @param baseFN base filename
 * @param T the text
 * @param [out] SA the suffix array
 * @param [out] BWT the BWT string
*/
void generateBWT(const string& T, vector<length_t>& SA, string& BWT) {   

    // build the BWT
    cout << "Generating BWT..." << endl;
    BWT = string(T.size(), '\0');
    for (size_t i = 0; i < SA.size(); i++) {
        if (SA[i] > 0) {
            BWT[i] = T[SA[i] - 1];
        } else {
            BWT[i] = T.back();
        }
    }
}

/**
 * Generate the suffix array and the BWT string of the reversed text
 * @param baseFN base filename
 * @param T the text
 * @param [out] revSA the suffix array of the reversed text
 * @param [out] revBWT the BWT string of the reversed text
*/
void generateRevSA(const string& baseFN, const string& T, 
                    vector<length_t>& revSA) {
    
    // read the reverse suffix array
    cout << "Reading " << baseFN << ".rev.sa..." << endl;
    readSA(baseFN + ".rev.sa", revSA, T.size());

    // perform a sanity check on the suffix array
    cout << "\tPerforming sanity checks..." << endl;
    sanityCheck(T, revSA);
    cout << "\tSanity checks OK" << endl;

}

void generateRevBWT(const string& T, vector<length_t>& revSA, string& revBWT) {

    // build the reverse BWT
    revBWT = string(T.size(), '\0');
    for (size_t i = 0; i < revSA.size(); i++) {
        if (revSA[i] > 0) {
            revBWT[i] = T[T.size() - revSA[i]];
        } else {
            revBWT[i] = T.front();
        }
    }
}

#endif //BUILDINDEX_H