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

#include "buildhelpers.h"
#include "bitvec.h"      // for Bitvec, Bitref
#include "definitions.h" // for length_t, BMOVE_BUILD_INDEX_TAG...
#include "libsais64.h"   // for libsais64
#include "logger.h"      // for Logger, logger

#include <algorithm> // for max, transform, equal, max_element, min_ele...
#include <array>     // for array
#include <cctype>    // for toupper
#include <cstdint>   // for int64_t, uint8_t, uint32_t
#include <iostream>  // for ifstream, ofstream, operator<<, ios, basic_...
#include <limits>    // for numeric_limits
#include <stdexcept> // for runtime_error

using namespace std;

// Random replacement function
char replaceNonACGT(char original, std::minstd_rand& gen,
                    const std::string& seed, size_t& seedIndex) {
    const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        std::uniform_int_distribution<size_t> distribution(
            0, validChars.length() - 1);
        return validChars[distribution(gen)];
    }
    return original;
}

// Seeded replacement function
char replaceNonACGTWithSeed(char original, std::minstd_rand& gen,
                            const std::string& seed, size_t& seedIndex) {
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        char replacement = seed[seedIndex];
        seedIndex = (seedIndex + 1) %
                    seed.length(); // Move to the next character in the seed
        return replacement;
    }
    seedIndex = 0; // Reset seed index
    return original;
}

void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             std::function<char(char, size_t&)> replaceFunc,
                             length_t expectedNumber) {

    logger.logInfo("Reading FASTA file " + fastaFile);

    // Open input file and get its size
    std::ifstream inputFile(fastaFile, std::ios::ate);

    if (!inputFile.is_open()) {
        throw std::runtime_error("Error opening file: " + fastaFile);
    }

    std::streamsize fileSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // Reserve memory for the concatenation string based on file size
    concatenation.reserve(static_cast<size_t>(fileSize));

    // Temporary variables for processing
    std::string line;
    length_t startPosition = 0;

    // Index for seed
    size_t seedIndex = 0;

    // Reserve memory for positions and seqNames vectors
    positions.reserve(expectedNumber);
    seqNames.reserve(expectedNumber);

    std::string sequence;

    while (std::getline(inputFile, line)) {
        if (line.empty())
            continue; // Skip empty lines

        if (line[0] == '>') { // Sequence name line
            // Process the previous sequence
            if (!sequence.empty()) {
                // Process the previous sequence
                for (char& c : sequence) {
                    c = replaceFunc(std::toupper(c), seedIndex);
                    concatenation += c;
                }

                positions.emplace_back(
                    startPosition); // push back the start position of this
                                    // sequence
                startPosition +=
                    sequence
                        .length(); // create start position for next sequence
                sequence.clear();  // Clear the sequence for the next one
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
    for (char& c : sequence) {
        c = replaceFunc(std::toupper(c), seedIndex);
        concatenation += c;
    }
    positions.emplace_back(startPosition);
    // add the end of the text
    positions.emplace_back(startPosition + sequence.length());

    // Close input file
    inputFile.close();

    logger.logInfo("Read " + std::to_string(seqNames.size()) + " sequences");
    logger.logInfo("Concatenating sequences...");

    // Ensure the concatenation ends with a dollar sign
    if (concatenation.back() != '$') {
        concatenation += '$';
    }

    logger.logInfo("Concatenation completed");
}

void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

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

void sanityCheck(const string& T, vector<length_t>& sa) {
    logger.logInfo("Performing sanity checks...");
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
    Bitvec bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T

    logger.logInfo("\tSanity checks OK");
}

void createAndWriteHeaderInfo(const string& baseFN,
                              const vector<string>& seqNames,
                              const vector<length_t>& positions) {
    logger.logInfo("Writing SAM header info for reference text to " + baseFN +
                   ".headerSN.bin...");
    // create the header lines with these reference sequences
    std::ofstream headerStream(baseFN + ".headerSN.bin", ios::binary);
    for (length_t i = 0; i < seqNames.size(); i++) {
        headerStream << "@SQ\tSN:" << seqNames[i]
                     << "\tLN:" << positions[i + 1] - positions[i] << "\n";
    }
    headerStream.close();
}

void writePositionsAndSequenceNames(const string& baseFN,
                                    const vector<length_t>& positions,
                                    const vector<string>& seqNames) {
    logger.logInfo("Write positions and sequence names to " + baseFN +
                   ".pos and " + baseFN + ".sna...");
    // Write the positions to disk
    std::ofstream ofs2(baseFN + ".pos", ios::binary);
    ofs2.write((char*)positions.data(), positions.size() * sizeof(length_t));
    ofs2.close();

    // Write the sequence names to disk
    std::ofstream ofs3(baseFN + ".sna", std::ios::binary);
    for (const auto& str : seqNames) {
        size_t len = str.size();
        ofs3.write(reinterpret_cast<const char*>(&len), sizeof(len));
        ofs3.write(str.c_str(), len);
    }
    ofs3.close();
}

void checkTextSize(const size_t tSize) {
    if (tSize > std::numeric_limits<length_t>::max()) {
        throw std::runtime_error(
            "The size of the text ()" + std::to_string(tSize) +
            ") is too large. Maximum size for the current "
            "compiler option is: " +
            std::to_string(std::numeric_limits<length_t>::max()));
    }
    if ((sizeof(length_t) * 8 == 64) &&
        (tSize <= std::numeric_limits<uint32_t>::max())) {
        logger.logWarning(
            "Program was compiled with 64-bit words, but the text size (" +
            std::to_string(tSize) +
            ") fits in a 32-bit word. Consider recompiling the program to "
            "improve performance and save memory.");
    }
}

void countChars(const string& T, vector<length_t>& charCounts) {
    // initialize the counts to 0 and 256 elements
    charCounts = vector<length_t>(256, 0);
    for (char c : T)
        charCounts[(unsigned char)c]++;
}

void writeCharCounts(const string& baseFN, const vector<length_t>& charCounts) {
    ofstream ofs(baseFN + ".cct", ios::binary);
    ofs.write((char*)charCounts.data(), charCounts.size() * sizeof(length_t));
    ofs.close();
    logger.logInfo("Wrote file " + baseFN + ".cct");
}

void createAlphabet(const string& T, const vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma) {
    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    logger.logInfo("Text has length " + std::to_string(T.size()));
    logger.logInfo("Text has " + std::to_string(nUniqueChar) +
                   " unique characters");

    if (nUniqueChar > ALPHABET) {
        logger.logError("FATAL: the number of unique characters in the "
                        "text exceeds the alphabet size. Please recompile "
                        "b-move using a higher value for ALPHABET");
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        logger.logWarning("the number of unique characters in the "
                          "text is less than the ALPHABET size specified when "
                          "b-move was compiled. Performance may be affected.");
    }

    sigma = Alphabet<ALPHABET>(charCounts);
}

void createSuffixArray(const string& T, vector<length_t>& SA) {
    logger.logInfo("Generating the suffix array using libsais...");

    // create the suffix array with libsais
    // Convert std::string to const uint8_t*
    const uint8_t* tPtr = reinterpret_cast<const uint8_t*>(T.c_str());

    // Length of the input string
    int64_t n = static_cast<int64_t>(T.size());

    // Suffix array, size should be n (or n+fs if you want extra space, but 0 is
    // enough for most cases)
    std::vector<int64_t> suffixArray64(n, 0);
    // Call the libsais64 function
    int64_t result = libsais64(tPtr, suffixArray64.data(), n, 0, nullptr);
    if (result != 0) {
        throw runtime_error("libsais64 failed with error code " +
                            to_string(result));
    }
    // cast int64_t to length_t
    SA = vector<length_t>(suffixArray64.begin(), suffixArray64.end());

    logger.logInfo("Suffix array generated successfully!");
}

void createRevBWT(const string& baseFN, const string& T,
                  const vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, string& rBWT) {

    // build the reverse BWT
    logger.logInfo("Generating BWT of reversed text...");
    rBWT.resize(T.size());
    for (size_t i = 0; i < revSA.size(); i++)
        if (revSA[i] > 0)
            rBWT[i] = T[T.size() - revSA[i]];
        else
            rBWT[i] = T.front();
}

void writeStringToBinaryFile(const std::string& name, const std::string& str){
      std::ofstream outFile(name, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << name << std::endl;
        return;
    }

    // First, write the size of the string
    length_t size = str.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    // Then, write the string data itself
    outFile.write(str.c_str(), size);
    outFile.close();
}

void preprocessFastaFile(const std::string& fastaFile,
                         const std::string& baseFN, string& T,
                         length_t seedLength, bool noWriting) {
    std::vector<length_t> positions; // the start positions of each
    std::vector<string> seqNames;    // the name of each sequence

    logger.logInfo("Preprocessing FASTA file " + fastaFile);

    // Seeded random number generation
    std::minstd_rand gen(42);

    // Generate random seed for seeded replacement
    std::string seed;
    size_t seedIndex = 0;
    for (length_t i = 0; i < seedLength; i++) {
        seed += replaceNonACGT('N', gen, seed, seedIndex);
    }

    // Reset gen to original seed
    gen.seed(42);

    std::function<char(char, size_t&)> replaceFunc;

    if (seedLength == 0) {

        replaceFunc = [&gen, &seed](char c, size_t& seedIndex) -> char {
            return replaceNonACGT(c, gen, seed, seedIndex);
        };

        logger.logInfo(
            "Using random (non-seeded) replacement for non-ACGT characters...");

    } else {

        replaceFunc = [&gen, &seed](char c, size_t& seedIndex) -> char {
            return replaceNonACGTWithSeed(c, gen, seed, seedIndex);
        };

        logger.logInfo("Using seeded replacement for non-ACGT characters, with "
                       "a seed length of " +
                       std::to_string(seedLength));
    }

    concatenateAndTransform(fastaFile, T, positions, seqNames, replaceFunc);
    checkTextSize(T.size());

    logger.logInfo("Converting to uppercase...");
    std::transform(T.begin(), T.end(), T.begin(), ::toupper);

    if (!noWriting) {
        logger.logInfo("Writing concatenated uppercase sequence to disk...");
        writeStringToBinaryFile(baseFN + ".txt.bin", T);
    }

    createAndWriteHeaderInfo(baseFN, seqNames, positions);
    writePositionsAndSequenceNames(baseFN, positions, seqNames);
}

void writeTagAndCompiledInfo(const string& baseFN) {
    ofstream tagFile(baseFN + ".tag");
    tagFile << BMOVE_BUILD_INDEX_TAG;
    tagFile.close();
    // Write the 64 or 32-bit compiled info to a file
    ofstream compiledInfoFile(baseFN + ".comp");
    compiledInfoFile << sizeof(length_t);
    compiledInfoFile.close();
}

void generateBWT(const string& T, const vector<length_t>& SA, string& BWT) {
    // build the BWT
    logger.logInfo("Generating BWT...");
    BWT.resize(T.size());
    for (size_t i = 0; i < SA.size(); i++)
        if (SA[i] > 0)
            BWT[i] = T[SA[i] - 1];
        else
            BWT[i] = T.back();
}

void determineInputMode(string& baseFN, bool& txtMode, bool& fastaMode,
                        std::string& fastaExtension) {

    std::array<string, 4> allowedExtensionsFasta = {".fasta", ".fa", ".FASTA",
                                                    ".FA"};
    std::array<string, 4> requiredExtensionsTxt = {".txt", ".sa", ".rev.sa",
                                                   ".rev.txt"};

    // check which mode to use
    txtMode = true;
    fastaMode = false;

    for (auto& ext : requiredExtensionsTxt) {
        string filename = baseFN + ext;

        // check if file with filename exists
        ifstream ifs(filename);
        if (!ifs) {
            // file does not exist
            txtMode = false;
            break;
        }
    }
    fastaExtension = ".fasta";

    // Check if any file with allowed fasta extensions
    // exists
    if (!fastaMode) {
        for (const auto& ext : allowedExtensionsFasta) {
            std::ifstream ifs(baseFN + ext);
            if (ifs) {
                fastaMode = true;
                fastaExtension = ext;
                logger.logDeveloper("Detected fasta extension: " +
                                    fastaExtension);
                break;
            }
        }
    }

    // Check if the base filename ends with any allowed
    // fasta extensions
    for (const auto& ext : allowedExtensionsFasta) {
        if (baseFN.size() >= ext.size() &&
            std::equal(ext.rbegin(), ext.rend(), baseFN.rbegin())) {
            fastaMode = true;
            fastaExtension = baseFN.substr(baseFN.size() - ext.size());
            logger.logDeveloper("Detected fasta extension: " + fastaExtension);
            baseFN = baseFN.substr(0, baseFN.size() - ext.size());
            break;
        }
    }

    if (txtMode) {
        // txt takes precedence as it has already
        // preprocessed suffix arrays
        fastaMode = false;
    }
}

void writeCharCountsAndCreateAlphabet(const string& baseFN, const string& T,
                                      Alphabet<ALPHABET>& sigma,
                                      vector<length_t>& charCounts) {
    // count the frequency of each characters in T
    countChars(T, charCounts);
    // write the character counts table
    writeCharCounts(baseFN, charCounts);
    // Create the alphabet
    createAlphabet(T, charCounts, sigma);
}

void readOrCreateSAWithSanityCheck(const string& baseFN, vector<length_t>& SA,
                                   const string& T, bool fromFasta) {
    if (fromFasta) {
        // create suffix array from concatenated using radixSA64
        createSuffixArray(T, SA);
    } else {
        logger.logInfo("Reading " + baseFN + ".sa...");
        readSA(baseFN + ".sa", SA, T.size());
    }

    // perform a sanity check on the suffix array
    sanityCheck(T, SA);
}

void readOrCreateRevSAWithSanityCheck(const string& baseFN,
                                      vector<length_t>& revSA, const string& T,
                                      bool fromFasta) {
    if (fromFasta) {
        std::string revT = T;
        std::reverse(revT.begin(), revT.end());

        // create the suffix array of the reverse text
        createSuffixArray(revT, revSA);
        revT.clear();

    } else {
        logger.logInfo("Reading " + baseFN + ".rev.sa...");
        readSA(baseFN + ".rev.sa", revSA, T.size());
    }
    // perform a sanity check on the suffix array
    sanityCheck(T, revSA);
}

void preprocessTextMode(const string& baseFN, string& T) {
    // read the text file from disk
    logger.logInfo("Reading " + baseFN + ".txt...");
    readText(baseFN + ".txt", T);
    if (!T.empty() && T.back() != '$') {
        T += '$';
    }
    checkTextSize(T.size());

    logger.logInfo("Converting text to uppercase...");
    std::transform(T.begin(), T.end(), T.begin(), ::toupper);

    writeStringToBinaryFile(baseFN + ".txt.bin", T);

    // create the sequence names vector, positions vector and the header
    vector<length_t> positions = {0, (length_t)T.size()};
    vector<string> seqNames = {baseFN};

    // create the header lines with these reference sequences
    createAndWriteHeaderInfo(baseFN, seqNames, positions);

    // Write the positions to disk
    writePositionsAndSequenceNames(baseFN, positions, seqNames);
}
