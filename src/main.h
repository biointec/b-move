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
#include "wordlength.h"
#include "nucleotide.h"
#include "bitparallelmatrix.h"
#include "indexhelpers.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <sstream> // used for splitting strings
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>

using namespace std;
vector<string> schemes = {"kuch1", "kuch2",  "kianfar", "manbest", "pigeon",
                          "01*0",  "custom", "naive",   "multiple"};

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

struct ReadRecord {
    string id;
    string read;

    string qual;

    /**
     * Deep copies the strings
     */
    ReadRecord(string id, string read, string qual)
        : id(id), read(read), qual(qual) {
    }
};

vector<ReadRecord> getReads(const string& file) {
    vector<ReadRecord> reads;
    reads.reserve(200000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream ifile(file.c_str());
    if (!ifile) {
        throw runtime_error("Cannot open file " + file);
    }
    if (!fasta && !fastq) {
        // this is a not readable

        throw runtime_error("extension " + extension +
                            " is not a valid extension for the readsfile");
    } else if (fasta) {
        // fasta file
        string read = "";
        string id = "";
        string qual = ""; // empty quality string for fasta
        string line;

        while (getline(ifile, line)) {
            if (line.empty()) {
                continue; // Skip empty lines
            }

            if (line[0] == '>' || line[0] == '@') {
                // This is an ID line
                if (!id.empty()) {
                    // If we already have data, process it and clear
                    reads.emplace_back(id, read, qual);
                    reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
                    id.clear();
                    read.clear();
                }
                id = line.substr(1); // Extract ID (skip '>')
            } else {
                // This is a sequence line
                read += line;
            }
        }

        // Process the last entry if it exists
        if (!id.empty()) {
            reads.emplace_back(id, read, qual);
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string qual = "";
        string plusLine = ""; // Skip the '+' line
        string line;

        while (getline(ifile, id) && getline(ifile, read) &&
               getline(ifile, plusLine) && // Skip the '+' line
               getline(ifile, qual)) {
            if (!id.empty() && id[0] != '@') {
                throw runtime_error("File " + file +
                                    "doesn't appear to be in FastQ format");
            }

            if (id.back() == '\n') {
                id.pop_back();
            }
            if (!read.empty() && read.back() == '\n') {
                read.pop_back();
            }

            assert(id.size() > 1);
            id = (id.substr(1));
            reads.emplace_back(id, read, qual);
            reverse(qual.begin(), qual.end());
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
            id.clear(), read.clear(), qual.clear();
        }
    }

    return reads;
}

void writeToOutput(const string& file, const vector<vector<TextOcc>>& mPerRead,
                   const vector<ReadRecord>& reads, length_t textLength,
                   const std::string& basefile) { 

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "@HD"
       << "\t"
       << "VN:1.6"
       << "\t"
       << "SO:queryname"
       << "\n";
    f2 << "@SQ"
       << "\t"
       << "SN:" << basefile << "\t"
       << "LN:" << textLength << "\n";

    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].id;

        for (auto m : mPerRead[i]) {
            f2 << m.getSamLine() << "\n";
        }

        for (auto m : mPerRead[i + 1]) {
            f2 << m.getSamLine() << "\n";
        }
    }

    f2.close();
}

double findMedian(vector<length_t> a, int n) {

    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(), a.begin() + (n - 1) / 2, a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Value at index (N/2)th
        // is the median
        return (double)a[n / 2];
    }
}
