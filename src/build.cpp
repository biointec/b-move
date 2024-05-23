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
#include "build.h"
#include <iostream>
#include <sstream>

#include "radixSA64/radix.h"


/**
 * Print usage to stdout
*/
void showUsage() {
    cout << "Usage: ./bmove-build <base filename>\n\n";
    cout << " [options]\n";
    cout << "  -v verbose output\n\n";

    cout << "Report bugs to jan.fostier@ugent.be" << endl;
}

/**
 * Parse command line arguments
 * @param argc number of arguments
 * @param argv char pointer array of arguments
 * @param [out] baseFN base filename
 * @param verbose verbose ouput
 * @return True if succeeded, false otherwise
*/
bool parseArguments(int argc, char* argv[], string& baseFN,
                    bool& verbose) {

    if (argc < 2) {
        cout << "Insufficient number of arguments\n";
        return false;
    }

    for (int i = 1; i < argc - 1; i++) {
        const string& arg = argv[i];
        if (arg == "-v") {
            verbose = true;
        } else {
            throw runtime_error("Unknown argument: " + arg);
        }
    }

    baseFN = argv[argc - 1];
    return true;
}

void writeIntVectorBinary(const std::string& filename, const std::vector<length_t>& array) {
    // convert to int_vector
    uint8_t width = (uint8_t)ceil(log2(*max_element(array.begin(), array.end())));
    sdsl::int_vector<> intVector(array.size(), 0, width);
    for (size_t i = 0; i < array.size(); i++) {
        intVector[i] = array[i];
    }
    std::ofstream ofs(filename, std::ios::binary);
    intVector.serialize(ofs);
    ofs.close();
}


/**
 * Builds samplesFirst and samplesLast.
 * @param samplesFirst Vector to fill. First sa samples of each input interval.
 * @param samplesLast Vector to fill. Last sa samples of each input interval.
 * @param SA The suffix array.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
*/
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast, vector<length_t>& SA, const string& BWT) {

    samplesFirst.push_back(SA[0] > 0 ? SA[0] - 1 : BWT.size() - 1);
    for (size_t pos = 0; pos < BWT.size()-1; pos++) {
        if (BWT[pos] != BWT[pos+1]) {
            samplesLast.push_back(SA[pos] > 0 ? SA[pos] - 1 : BWT.size() - 1);
            samplesFirst.push_back(SA[pos + 1] > 0 ? SA[pos + 1] - 1 : BWT.size() - 1);
        }
    }
    samplesLast.push_back(SA[BWT.size() - 1] > 0 ? SA[BWT.size() - 1] - 1 : BWT.size() - 1);

}


/**
 * Generate predecessor structures
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 * @param verbose print verbose output to stdout
*/
void generatePredecessors(const vector<length_t>& samplesFirst, 
    const vector<length_t>& samplesLast, SparseBitvec& predFirst, 
    SparseBitvec& predLast, length_t textLength) {

    vector<bool> predFirstBV(textLength, false);
    vector<bool> predLastBV(textLength, false);

    for (length_t SApos: samplesFirst) {
        predFirstBV[SApos] = true;
    }

    for (length_t SApos: samplesLast) {
        predLastBV[SApos] = true;
    }

    predFirst = SparseBitvec(predFirstBV);
    predLast = SparseBitvec(predLastBV);
}

/**
 * Generate predToRun arrays
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] firstToRun mapping between rank of ones in predFirst bitvector and run indices
 * @param [out] lastToRun mapping between rank of ones in predLast bitvector and run indices
*/
void generatePredTorun(const vector<length_t>& samplesFirst, const vector<length_t>& samplesLast,
    vector<length_t>& firstToRun, vector<length_t>& lastToRun) {

    firstToRun.resize(samplesFirst.size());
    iota(firstToRun.begin(), firstToRun.end(), 0);
    sort(firstToRun.begin(), firstToRun.end(), 
        [&samplesFirst](length_t a, length_t b){return samplesFirst[a] < samplesFirst[b];});
    

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(), 
        [&samplesLast](length_t a, length_t b){return samplesLast[a] < samplesLast[b];});
}

/**
 * Create all datastructures necessary for bmove
 * and write them to their respective files
 * @param baseFN base filename, all files will be generated as <baseFN>.<ext>
 * @param verbose print verbose output to stdout
*/
void createBmoveFromIntermediateFiles(const string& baseFN, bool verbose) {

    string T;
    vector<length_t> charCounts;
    Alphabet<ALPHABET> sigma;
    generateText(baseFN, T);
    generateAlphabet(T, charCounts, sigma);

    // write the character counts table
    writeArrayBinary(baseFN + ".cct", charCounts);
    cout << "Wrote file " << baseFN << ".cct\n";

    vector<length_t> SA;
    string BWT;
    generateSA(baseFN, T, SA);
    generateBWT(T, SA, BWT);

    // generate and write PLCP array
    cout << "Generating PLCP array..." << endl;
    PLCP plcp(T, SA, verbose);
    plcp.write(baseFN + ".plcp");
    cout << "Wrote file " << baseFN << ".plcp\n";

    // Create samplesLast and samplesFirst
    cout << "Generating samplesFirst and samplesLast..." << endl;
    vector<length_t> samplesFirst;
    vector<length_t> samplesLast;
    buildSamples(samplesFirst, samplesLast, SA, BWT);
    SA.clear();
    // write samplesFirst and samplesLast to file
    writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
    cout << "Wrote file " << baseFN << ".smpf\n";
    writeIntVectorBinary(baseFN + ".smpl", samplesLast);
    cout << "Wrote file " << baseFN << ".smpl\n";

    SparseBitvec predFirst;
    SparseBitvec predLast;
    vector<length_t> firstToRun;
    vector<length_t> lastToRun;
    cout << "Generating predecessor structures..." << endl;
    generatePredecessors(samplesFirst, samplesLast, predFirst, predLast, BWT.size());
    cout << "Generating predToRun arrays..." << endl;
    generatePredTorun(samplesFirst, samplesLast, firstToRun, lastToRun);
    // Write predecessor structures
    predFirst.write(baseFN + ".prdf");
    cout << "Wrote file " << baseFN << ".prdf\n";
    predLast.write(baseFN + ".prdl");
    cout << "Wrote file " << baseFN << ".prdl\n";
    writeIntVectorBinary(baseFN + ".ftr", firstToRun);
    cout << "Wrote file " << baseFN << ".ftr\n";
    writeIntVectorBinary(baseFN + ".ltr", lastToRun);
    cout << "Wrote file " << baseFN << ".ltr\n";
    samplesFirst.clear();
    samplesLast.clear();
    firstToRun.clear();
    lastToRun.clear();

    cout << "Creating move structure..." << endl;
    createMove<ALPHABET>(baseFN, BWT, charCounts, sigma, SA);
    cout << "Original text length: " << T.size() << "\n";
    BWT.clear();

    vector<length_t> revSA;
    string revBWT;
    generateRevSA(baseFN, T, revSA);
    generateRevBWT(T, revSA, revBWT);

    // Create samplesLast and samplesFirst
    vector<length_t> revSamplesFirst;
    vector<length_t> revSamplesLast;
    cout << "Generating reverse samplesFirst and reverse samplesLast..." << endl;
    buildSamples(revSamplesFirst, revSamplesLast, revSA, revBWT);
    // write samplesFirst and samplesLast to file
    writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
    cout << "Wrote file " << baseFN << ".rev.smpf\n";
    writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
    cout << "Wrote file " << baseFN << ".rev.smpl\n";
    revSamplesFirst.clear();
    revSamplesLast.clear();

    cout << "Creating reverse move structure..." << endl;
    createMove<ALPHABET>(baseFN + ".rev", revBWT, charCounts, sigma, revSA);
    cout << "Original text length: " << T.size() << "\n";
    
}

void processFastaFileIntoBmove(const std::string& fastafile, const std::string& baseFN) {
    std::string T = "";                   // the (concatenated) text
    std::vector<length_t> positions; // the start positions of each sequence
    std::vector<string> seqNames;    // the name of each sequence

    std::cout << "Processing fasta file " << fastafile << std::endl;
    concatenateAndTransform(fastafile, T, positions, seqNames, baseFN + ".fa.meta");

    std::cout << "Concatenated " << seqNames.size() << " sequences" << endl;

    // Write the concatenated sequence to disk
    std::ofstream ofs(baseFN + ".txt");
    ofs << T;
    ofs.close();

    // Write the positions to disk
    std::ofstream ofs2(baseFN + ".pos", ios::binary);

    ofs2.write((char*)positions.data(), positions.size() * sizeof(length_t));
    ofs2.close();
    // clear the vector
    positions.clear();

    // Write the sequence names to disk

    std::ofstream ofs3(baseFN + ".sna", std::ios::binary);
    for (const auto& str : seqNames) {
        size_t len = str.size();
        ofs3.write(reinterpret_cast<const char*>(&len), sizeof(len));
        ofs3.write(str.c_str(), len);
    }
    ofs3.close();
    // clear the vector
    seqNames.clear();

    // count the frequency of each characters in T
    vector<length_t> charCounts(256, 0);
    for (char c : T)
        charCounts[(unsigned char)c]++;
    
    // write the character counts table
    writeArrayBinary(baseFN + ".cct", charCounts);
    cout << "Wrote file " << baseFN << ".cct\n";

    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

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

    Alphabet<ALPHABET> sigma(charCounts);

    // create suffix array from concatenated using radixSA64
    cout << "Generating the suffix array using radixSA64...\n";
    length_t textLength = T.size();
    clock_t startTime = clock();
    length_t* radixSA = Radix<length_t>((uchar*)&(T)[0], T.size(), 0).build();

    std::vector<length_t> SA(radixSA, radixSA + textLength);

    delete[] radixSA;
    clock_t endTime = clock();
    float mseconds = clock_diff_to_msec(endTime - startTime);
    printf("radixSA64 took [%.2fs]\n", mseconds / 1000.0);

    // perform a sanity check on the suffix array
    cout << "\tPerforming sanity checks..." << endl;
    sanityCheck(T, SA);
    cout << "\tSanity checks OK" << endl;

    // build the BWT
    cout << "Generating BWT..." << endl;
    string BWT;
    generateBWT(T, SA, BWT);

    // generate and write PLCP array
    cout << "Generating PLCP array..." << endl;
    PLCP plcp(T, SA, false);
    plcp.write(baseFN + ".plcp");
    cout << "Wrote file " << baseFN << ".plcp\n";

    // Create samplesLast and samplesFirst
    cout << "Generating samplesFirst and samplesLast..." << endl;
    vector<length_t> samplesFirst;
    vector<length_t> samplesLast;
    buildSamples(samplesFirst, samplesLast, SA, BWT);
    SA.clear();
    // write samplesFirst and samplesLast to file
    writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
    cout << "Wrote file " << baseFN << ".smpf\n";
    writeIntVectorBinary(baseFN + ".smpl", samplesLast);
    cout << "Wrote file " << baseFN << ".smpl\n";

    SparseBitvec predFirst;
    SparseBitvec predLast;
    vector<length_t> firstToRun;
    vector<length_t> lastToRun;
    cout << "Generating predecessor structures..." << endl;
    generatePredecessors(samplesFirst, samplesLast, predFirst, predLast, BWT.size());
    cout << "Generating predToRun arrays..." << endl;
    generatePredTorun(samplesFirst, samplesLast, firstToRun, lastToRun);
    // Write predecessor structures
    predFirst.write(baseFN + ".prdf");
    cout << "Wrote file " << baseFN << ".prdf\n";
    predLast.write(baseFN + ".prdl");
    cout << "Wrote file " << baseFN << ".prdl\n";
    writeIntVectorBinary(baseFN + ".ftr", firstToRun);
    cout << "Wrote file " << baseFN << ".ftr\n";
    writeIntVectorBinary(baseFN + ".ltr", lastToRun);
    cout << "Wrote file " << baseFN << ".ltr\n";
    samplesFirst.clear();
    samplesLast.clear();
    firstToRun.clear();
    lastToRun.clear();

    cout << "Creating move structure..." << endl;
    createMove<ALPHABET>(baseFN, BWT, charCounts, sigma, SA);
    cout << "Original text length: " << T.size() << "\n";
    BWT.clear();
    SA.clear();

    // create the reverse of concatenation
    {
        std::string revT = T;
        std::reverse(revT.begin(), revT.end());

        // write reverse concatenation to disk

        std::ofstream ofs4(baseFN + ".rev.txt");
        ofs4 << revT;
        ofs4.close();

        // create the suffix array of the reverse text
        cout << "Generating the suffix array of the reverse text using "
                "radixSA64...\n";
        clock_t startTime2 = clock();
        length_t* revRadixSA =
            Radix<length_t>((uchar*)&(revT)[0], revT.size(), 0).build();
        std::vector<length_t> revSA(revRadixSA, revRadixSA + textLength);

        delete[] revRadixSA;
        clock_t endTime2 = clock();
        float mseconds = clock_diff_to_msec(endTime2 - startTime2);
        printf("radixSA64 took [%.2fs]\n", mseconds / 1000.0);

        // perform a sanity check on the suffix array
        cout << "\tPerforming sanity checks..." << endl;
        sanityCheck(T, revSA);
        cout << "\tSanity checks OK" << endl;

        // build the reverse BWT
        string rBWT(T.size(), '\0');
        generateRevBWT(T, revSA, rBWT);

        // Create samplesLast and samplesFirst
        vector<length_t> revSamplesFirst;
        vector<length_t> revSamplesLast;
        cout << "Generating reverse samplesFirst and reverse samplesLast..." << endl;
        buildSamples(revSamplesFirst, revSamplesLast, revSA, rBWT);
        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        cout << "Wrote file " << baseFN << ".rev.smpf\n";
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        cout << "Wrote file " << baseFN << ".rev.smpl\n";
        revSamplesFirst.clear();
        revSamplesLast.clear();

        cout << "Creating reverse move structure..." << endl;
        createMove<ALPHABET>(baseFN + ".rev", rBWT, charCounts, sigma, revSA);
        cout << "Original text length: " << T.size() << "\n";
    }
}

int main(int argc, char* argv[]) {
    string baseFN;

    std::array<string, 4> allowedExtensionsFasta = {".fasta", ".fa", ".FASTA",
                                                    ".FA"};
    std::array<string, 4> requiredExtensionsTxt = {".txt", ".sa", ".rev.sa",
                                                   ".rev.txt"};
    bool verbose = false;

    if (!parseArguments(argc, argv, baseFN, verbose)) {
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to the bidirectional move structure construction!\n";
    cout << "Alphabet size is " << ALPHABET - 1 << " + 1\n";
    // check which mode to use
    bool txtMode = true;
    bool fastaMode = true;

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

    // to do set fastamode to true if any of the allowed extensions
    for (auto& ext : allowedExtensionsFasta) {
        string filename = baseFN + ext;
        // check if file with filename exists
        ifstream ifs(filename);
        if (ifs) {
            // file does exist -> use fasta mode
            fastaMode = true;
            break;
        }
    }

    // check if the basefile aready ends with fasta extension

    if (baseFN.size() >= 6) {
        string posExt = baseFN.substr(baseFN.size() - 6, 6);
        for (char& c : posExt) {
            c = std::tolower(c);
        }
        if (posExt == ".fasta") {
            fastaMode = true;
            baseFN = baseFN.substr(0, baseFN.size() - 6);
        }
    }

    if (txtMode) {
        // txt takes precedence as it has already preprocesed suffix arrays
        fastaMode = false;
    }

    try {
        if (txtMode) {
            createBmoveFromIntermediateFiles(baseFN, verbose);
        } else if (fastaMode) {
            processFastaFileIntoBmove(baseFN + ".fasta", baseFN);
        } else {
            throw runtime_error("No valid input files found");
        }

    } catch (const std::exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}