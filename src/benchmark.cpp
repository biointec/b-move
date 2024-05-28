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
#include "main.h"
#include "move.h"
#include "rindex.h"
#include "rsearchstrategy.h"


void doBench(const vector<ReadRecord>& reads, RIndex& index, 
             RSearchStrategy* strategy, length_t ED,
             const std::string& basefile) {

    size_t totalUniqueMatches = 0, sizes = 0, mappedReads = 0;

    cout << "Benchmarking with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getDistanceMetric() << " distance " << endl;
    cout.precision(2);

    vector<vector<TextOcc>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    std::vector<length_t> numberMatchesPerRead;
    numberMatchesPerRead.reserve(reads.size());

    Counters counters;
    counters.resetCounters();

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 2) {

        const auto& readrecord = reads[i];
        const auto& revReadrecord = reads[i + 1];

        const string& id = readrecord.id;
        const string& read = readrecord.read;
        const string& revCompl = revReadrecord.read;
        const string& qual = readrecord.qual;
        const string& revQual = revReadrecord.qual;

        if (((i >> 1) - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i / 2 << "/" << reads.size() / 2 << "\r";
            cout.flush();
        }

        sizes += read.size();

        auto matches =
            strategy->matchApprox(read, ED, counters, id, qual, false);
        totalUniqueMatches += matches;

        // do the same for the reverse complement
        auto matchesRevCompl =
            strategy->matchApprox(revCompl, ED, counters, id, revQual, true);
        totalUniqueMatches += matchesRevCompl;

        // keep track of the number of mapped reads
        mappedReads += !((matchesRevCompl==0) && (matches==0));

        matchesPerRead.emplace_back(matches);
        matchesPerRead.emplace_back(matchesRevCompl);

        numberMatchesPerRead.emplace_back(matches + matchesRevCompl);
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";
    cout << "Average no. nodes: " << counters.nodeCounter / (reads.size() / 2.0)
         << endl;
    cout << "Total no. Nodes: " << counters.nodeCounter << "\n";

    cout << "Average no. unique matches: "
         << totalUniqueMatches / (reads.size() / 2.0) << endl;
    cout << "Total no. unique matches: " << totalUniqueMatches << "\n";
    cout << "Average no. reported matches "
         << counters.totalReportedPositions / (reads.size() / 2.0) << endl;
    cout << "Total no. reported matches: " << counters.totalReportedPositions
         << "\n";
    cout << "Mapped reads: " << mappedReads << endl;
    cout << "Median number of occurrences per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;

    cout << "Average size of reads: " << sizes / (reads.size() / 2.0) << endl;
    
    cout << "Number of LF queries: " << LF_call_count << endl;
    cout << "Average number of CPU cycles per LF query: " << elapsed_LF*1.0 / LF_call_count << endl;
}

void showUsage() { 
    cout << "Usage: ./bmove-benchmarkCharExt [options] basefilename readfile.[ext]\n\n"
            "[options]\n"
            "  -e  --max-errors       maximum edit distance [default = 0]\n"
            "  -p  --partitioning     Add flag to do uniform/static/dynamic \n"
            "                          partitioning [default = dynamic]\n"
            "  -m  --metric           Add flag to set distance metric \n"
            "                          (edit/hamming) [default = edit]\n"
            "  -ks --kmer-size        The size of the seeds for dynamic \n"
            "                          partitioning [default = 10]\n"
            "  -ss --search-scheme    Choose the search scheme [default = kuch1]\n"
            "                          options:\n"
            "                            kuch1    Kucherov k + 1\n"
            "                            kuch2    Kucherov k + 2\n"
            "                            kianfar  Optimal Kianfar scheme\n"
            "                            manbest  Manual best improvement for \n"
            "                                     kianfar scheme (only for ed = 4)\n"
            "                            pigeon   Pigeon hole scheme\n"
            "                            01*0     01*0 search scheme\n"
            "                            custom   Custom search scheme, the next \n"
            "                                     parameter should be a path to \n"
            "                                     the folder containing this search \n"
            "                                     scheme\n"
            "                            multiple Multiple search scheme, the next \n"
            "                                     parameter should be a path to \n"
            "                                     the folder containing the different \n"
            "                                     search schemes to choose from with \n"
            "                                     dynamic selection.\n\n"
            "[ext]\n"
            "  One of the following: fq, fastq, FASTA, fasta, fa\n\n"
            "Following input files are required:\n"
            "  <base filename>.cct:      character counts table\n"
            "  <base filename>.move:     move table of T\n"
            "  <base filename>.rev.move: move table of the reverse of T\n\n";

    cout << "Report bugs to jan.fostier@ugent.be" << endl;
}

int main(int argc, char* argv[]) {

    std::cout << "Using " << LENGTH_TYPE_NAME << std::endl;

    int requiredArguments = 2; // baseFile of files and file containing reads

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments" << endl;
        showUsage();
        return EXIT_FAILURE;
    }
    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        showUsage();
        return EXIT_SUCCESS;
    }

    cout << "Welcome to the bidirectional move structure!\n";

    string maxED = "0";
    string searchscheme = "kuch1";
    string customFile = "";

    string kmerSize = "10";

    PartitionStrategy pStrat = DYNAMIC;
    DistanceMetric metric = EDITOPTIMIZED;

    bool printMemUsage = false;

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-p" || arg == "--partitioning") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "uniform") {
                    pStrat = UNIFORM;
                } else if (s == "dynamic") {
                    pStrat = DYNAMIC;
                } else if (s == "static") {
                    pStrat = STATIC;
                } else {
                    throw runtime_error(
                        s + " is not a partitioning option\nOptions are: "
                            "uniform, static, dynamic");
                }
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-e" || arg == "--max-errors") {
            if (i + 1 < argc) {
                maxED = argv[++i];
            }
        } else if (arg == "-ss" || arg == "--search-scheme") {
            if (i + 1 < argc) {
                searchscheme = argv[++i];
                if (find(schemes.begin(), schemes.end(), searchscheme) ==
                    schemes.end()) {
                    throw runtime_error(searchscheme +
                                        " is not on option as search scheme");
                }
                if (searchscheme == "custom" || searchscheme == "multiple") {
                    if (i + 1 < argc) {
                        customFile = argv[++i];
                    } else {
                        throw runtime_error("custom/multiple search scheme "
                                            "takes a folder as argument");
                    }
                }
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-m" || arg == "--metric") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "edit") {
                    metric = EDITOPTIMIZED;
                } else if (s == "hamming") {
                    metric = HAMMING;
                } else {
                    throw runtime_error(s +
                                        " is not a metric option\nOptions are: "
                                        "edit, hamming");
                }
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-ks" || arg == "--kmer-size") {
            if (i + 1 < argc) {
                kmerSize = argv[++i];
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } // else if (arg == "-mu" || arg == "--mem-usage") {
        //     printMemUsage = true;
        // }
        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return false;
        }
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > MAX_K) {
        cerr << ed << " is not allowed as maxED should be in [0, " << MAX_K
             << "]" << endl;

        return EXIT_FAILURE;
    }

    length_t kMerSize = stoi(kmerSize);

    if (ed > 4 && searchscheme != "custom" && searchscheme != "multiple" &&
        searchscheme != "naive") {
        throw runtime_error(
            "Hard-coded search schemes are only available for "
            "up to 4 errors. Use a custom search scheme instead.");
    }

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string baseFile = argv[argc - 2];
    string readsFile = argv[argc - 1];

    RIndex index = RIndex(baseFile, true, kMerSize);

    if (printMemUsage) {
        index.printMemSize();
    }

    cout << "Reading in reads from " << readsFile << endl;
    vector<ReadRecord> reads;
    try {
        reads = getReads(readsFile);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }

    RSearchStrategy* strategy;
    if (searchscheme == "kuch1") {
        strategy = new RKucherovKplus1(index, pStrat, metric);
    } else if (searchscheme == "kuch2") {
        strategy = new RKucherovKplus2(index, pStrat, metric);
    } else if (searchscheme == "kianfar") {
        strategy = new ROptimalKianfar(index, pStrat, metric);
    } else if (searchscheme == "manbest") {
        strategy = new RManBestStrategy(index, pStrat, metric);
    } else if (searchscheme == "01*0") {
        strategy = new RO1StarSearchStrategy(index, pStrat, metric);
    } else if (searchscheme == "pigeon") {
        strategy = new RPigeonHoleSearchStrategy(index, pStrat, metric);
    } else if (searchscheme == "custom") {
        strategy = new CustomRSearchStrategy(index, customFile, pStrat, metric);
    } else if (searchscheme == "multiple") {
        strategy =
            new RMultipleSchemesStrategy(index, customFile, pStrat, metric);
    } else if (searchscheme == "naive") {
        strategy = new RNaiveBackTrackingStrategy(index, pStrat, metric);
    } else {
        // should not get here
        throw runtime_error(searchscheme + " is not on option as search scheme");
    }

    doBench(reads, index, strategy, ed, baseFile);
    delete strategy;
    cout << "Bye...\n";
}