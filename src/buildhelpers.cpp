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
#include "libsais64.h" // for libsais
#include <iostream>

// Random replacement function
char replaceNonACGT(char original, std::minstd_rand& gen,
                    length_t& replacementCounter) {
    const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        replacementCounter++;
        return validChars[gen() % validChars.length()];
    }
    return original;
}

void processSequence(std::string& concatenation, std::string& sequence,
                     length_t& minRemovalLength, length_t& startPosition,
                     std::vector<length_t>& positions,
                     std::vector<std::string>& seqNames, std::minstd_rand& gen,
                     length_t& totalNsRemoved, length_t& totalNsReplaced,
                     length_t& totalOthersReplaced) {

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
                        concatenation += replaceNonACGT('N', gen, NsReplaced);
                    }
                } else {
                    NsRemoved += charCount - nRunStart;
                }
                inNRun = false;
            }
            // Replace non-ACGT characters
            c = replaceNonACGT(std::toupper(c), gen, othersReplaced);
            concatenation += c;
        }
        charCount++;
    }

    if (inNRun) {
        if (charCount - nRunStart < minRemovalLength) {
            // Replace N with a random ACGT character
            for (size_t i = 0; i < charCount - nRunStart; ++i) {
                concatenation += replaceNonACGT('N', gen, NsReplaced);
            }
        } else {
            NsRemoved += charCount - nRunStart;
        }
    }

    positions.emplace_back(
        startPosition); // push back the start position of this sequence
    startPosition +=
        sequence.length(); // create start position for next sequence
    sequence.clear();      // Clear the sequence for the next one

    totalNsRemoved += NsRemoved;
    totalNsReplaced += NsReplaced;
    totalOthersReplaced += othersReplaced;
}

void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             length_t expectedNumber) {

    std::cout << "Reading FASTA file " + fastaFile << std::endl;
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

    // Seed for random number generation
    std::minstd_rand gen(42);

    // Reserve memory for positions and seqNames vectors
    positions.reserve(expectedNumber);
    seqNames.reserve(expectedNumber);

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
                processSequence(concatenation, sequence, minRemovalLength,
                                startPosition, positions, seqNames, gen,
                                totalNsRemoved, totalNsReplaced,
                                totalOthersReplaced);
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
        processSequence(concatenation, sequence, minRemovalLength,
                        startPosition, positions, seqNames, gen, totalNsRemoved,
                        totalNsReplaced, totalOthersReplaced);
    }

    // Close input file
    inputFile.close();

    std::cout << "Read " + std::to_string(seqNames.size()) + " sequences"
              << std::endl;
    std::cout << "Concatenating sequences..." << std::endl;
    // Ensure the concatenation ends with a dollar sign
    if (concatenation.back() != '$') {
        concatenation += '$';
    }

    std::cout << "Concatenation completed" << std::endl;

    cout << "\tTotal removed Ns: " << totalNsRemoved << "\n";
    cout << "\tTotal replaced Ns: " << totalNsReplaced << "\n";
    cout << "\tTotal replaced non-ACGT characters: " << totalOthersReplaced
         << "\n";
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
    sdsl::bit_vector bv(sa.size());
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

void createAndWriteHeaderInfo(const string& baseFN,
                              const vector<string>& seqNames,
                              const std::vector<length_t>& positions) {
    std::cout << "Writing SAM header info for reference text to " + baseFN +
                     ".headerSN.bin..."
              << std::endl;
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
    std::cout << "Write positions and sequence names to " + baseFN +
                     ".pos and " + baseFN + ".sna..."
              << std::endl;
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
        std::cout
            << "Program was compiled with 64-bit words, but the text size (" +
                   std::to_string(tSize) +
                   ") fits in a 32-bit word. Consider recompiling the program "
                   "to "
                   "improve "
                   "performance and save memory."
            << std::endl;
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
    std::cout << "Wrote file " + baseFN + ".cct" << std::endl;
}

void createAlphabet(const string& T, const vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma) {
    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    std::cout << "Text has length " + std::to_string(T.size()) << std::endl;
    cout << "Text has " + std::to_string(nUniqueChar) + " unique characters"
         << endl;

    if (nUniqueChar > ALPHABET) {
        std::cerr << "FATAL: the number of unique characters in the "
                     "text exceeds the alphabet size. Please recompile "
                     "b-move using a higher value for ALPHABET"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        std::cout << "the number of unique characters in the "
                     "text is less than the ALPHABET size specified when "
                     "b-move was compiled. Performance may be affected."
                  << std::endl;
    }

    sigma = Alphabet<ALPHABET>(charCounts);
}

void createSuffixArray(const string& T, vector<length_t>& SA) {
    std::cout << "Generating the suffix array using libsais..." << std::endl;
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

    std::cout << "Suffix array generated successfully!" << std::endl;
}

void createRevBWT(const string& baseFN, const string& T,
                  const vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, string& rBWT) {

    // build the reverse BWT
    std::cout << "Generating BWT of reversed text..." << std::endl;
    rBWT.resize(T.size());
    for (size_t i = 0; i < revSA.size(); i++)
        if (revSA[i] > 0)
            rBWT[i] = T[T.size() - revSA[i]];
        else
            rBWT[i] = T.front();
}

void preprocessFastaFile(const std::string& fastaFile,
                         const std::string& baseFN, string& T, bool noWriting) {
    std::vector<length_t> positions; // the start positions of each
    std::vector<string> seqNames;    // the name of each sequence

    std::cout << "Preprocessing FASTA file " + fastaFile << std::endl;

    concatenateAndTransform(fastaFile, T, positions, seqNames);
    checkTextSize(T.size());

    std::cout << "Converting to uppercase..." << std::endl;
    std::transform(T.begin(), T.end(), T.begin(), ::toupper);

    if (!noWriting) {
        std::cout << "Writing concatenated uppercase sequence to disk..."
                  << std::endl;
        std::ofstream ofs(baseFN + ".txt");
        ofs << T;
        ofs.close();
    }

    createAndWriteHeaderInfo(baseFN, seqNames, positions);
    writePositionsAndSequenceNames(baseFN, positions, seqNames);
}