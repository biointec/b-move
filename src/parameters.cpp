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

#include "parameters.h"
#include "logger.h"         // for Logger, logger
#include "rsearchstrategy.h" // for SearchStrategy, ColumbaSearchStrategy

#include <algorithm>     // for find, max
#include <cassert>       // for assert
#include <iomanip>       // for setw
#include <iostream>      // for operator<<, basic_ostream, cout, ostream
#include <memory>        // for shared_ptr
#include <stddef.h>      // for size_t
#include <stdexcept>     // for runtime_error
#include <stdint.h>      // for uint16_t
#include <stdlib.h>      // for exit
#include <string>        // for string, operator==, operator+, char_tr...
#include <thread>        // for thread
#include <unordered_map> // for unordered_map
#include <unordered_set> // for unordered_set
#include <vector>        // for vector

const std::vector<std::string> schemes = {
    "kuch1",  "kuch2", "kianfar",  "pigeon", "01*0",
    "custom", "naive", "multiple", "minU",   "columba"};

/**
 * Option for the command line arguments considering the partitioning strategy.
 */
class PartitioningOption : public Option {
  public:
    PartitioningOption() : Option("p", "partitioning", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "uniform") {
            params.pStrategy = UNIFORM;
        } else if (arg == "dynamic") {
            params.pStrategy = DYNAMIC;
        } else if (arg == "static") {
            params.pStrategy = STATIC;
        } else {
            logger.logWarning(arg +
                              " is not a partitioning option\nOptions "
                              "are: uniform, static, dynamic" +
                              ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Partitioning strategy to use. Options are: uniform, "
               "static, "
               "dynamic. Default is dynamic.";
    }
};

/**
 * Option for the command line arguments considering the search scheme to be
 * used.
 */
class SearchSchemeOption : public Option {
  public:
    SearchSchemeOption()
        : Option("S", "search-scheme", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.searchScheme = arg;
        if (find(schemes.begin(), schemes.end(), params.searchScheme) ==
            schemes.end()) {
            logger.logWarning(params.searchScheme +
                              " is not an option as a search scheme." +
                              ignoreMessage());
            params.searchScheme = "columba";
        }
    }

    std::string getDescription() const override {
        std::string help = "Search scheme to use. Options are: ";
        help += schemes.front();
        for (size_t i = 1; i < schemes.size(); ++i) {
            help += ", " + schemes[i];
        }
        help += ". Default is columba.";
        return help;
    }
};

/**
 * Option for the command line arguments considering the custom search scheme to
 * be used.
 */
class CustomOption : public Option {
  public:
    CustomOption() : Option("c", "custom", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.custom = arg;
    }

    std::string getDescription() const override {
        return "Path to custom search scheme (overrides default "
               "search scheme).";
    }
};

/**
 * Option for the command line arguments considering the dynamic selection
 * custom collection of search schemes option to be used.
 */
class DynamicSelectionOption : public Option {
  public:
    DynamicSelectionOption()
        : Option("d", "dynamic-selection", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.dynamicSelectionPath = arg;
    }

    std::string getDescription() const override {
        return "Path to custom search scheme with dynamic selection "
               "(overrides "
               "default search scheme).";
    }
};

/**
 * Option for the command line arguments considering the distance metric to be
 * used.
 */
class MetricOption : public Option {
  public:
    MetricOption() : Option("m", "metric", true, STRING, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "edit") {
            params.metric = EDIT;
        } else if (arg == "hamming") {
            params.metric = HAMMING;
        } else {
            logger.logWarning(arg +
                              " is not a distance metric option\nOptions are: "
                              "edit, hamming" +
                              ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Distance metric to use. Options are: edit, hamming. Default is "
               "edit.";
    }
};

/**
 * Option for the command line arguments considering the log file to be used.
 */
class LogFileOption : public Option {
  public:
    LogFileOption() : Option("l", "log-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.logFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the log file. Default is stdout.";
    }
};

/**
 * Option for the command line arguments considering the first reads file to be
 * used.
 */
class FirstReadsOption : public Option {
  public:
    FirstReadsOption()
        : Option("f", "reads-file", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.firstReadsFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the reads file.";
    }
};

/**
 * Option for the command line arguments considering the reference to be
 * used.
 */
class ReferenceOption : public Option {
  public:
    ReferenceOption() : Option("r", "reference", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.base = arg;
    }

    std::string getDescription() const override {
        return "Path to the basename of the index.";
    }
};

/**
 * Option for the command line arguments considering the output file to be
 * used.
 */
class OutputFileOption : public Option {
  public:
    OutputFileOption() : Option("o", "output-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.outputFile = arg;
        if (params.outputFile.size() > 3 &&
            params.outputFile.substr(params.outputFile.size() - 4) != ".sam") {
            params.outputIsSAM = false;
        }
    }

    std::string getDescription() const override {
        return "Path to the output file. Should be .sam or .rhs. Default is "
               "bmoveOutput.sam.";
    }
};

/**
 * Option for the command line arguments considering the maximal distance to be
 * used.
 */
class MaxDistanceOption : public Option {
  public:
    MaxDistanceOption() : Option("e", "max-distance", true, STRING, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.maxDistance = std::stoi(arg);
        } catch (...) {
            logger.logWarning("Max Distance should be a positive integer" +
                              ignoreMessage());
        }
        if (params.maxDistance < 0 || params.maxDistance > MAX_K) {
            logger.logWarning("Max Distance should be in [0, " +
                              std::to_string(MAX_K) + "]" + ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "The maximum allowed distance. Default is 0.";
    }
};

/**
 * Option for the command line arguments considering the k-mer size to be
 * used.
 */
class KmerSizeOption : public Option {
  public:
    KmerSizeOption() : Option("K", "kmer-size", true, INTEGER, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultKmerSize = params.kmerSize;
        try {
            params.kmerSize = std::stoi(arg);
        } catch (...) {
            logger.logWarning("kmerSkipSize should be an integer." +
                              ignoreMessage());
        }
        if (params.kmerSize < 0 || params.kmerSize > 15) {
            logger.logWarning("kmerSkipSize should be in [0, 15]." +
                              ignoreMessage());
            params.kmerSize = defaultKmerSize;
        }
    }

    std::string getDescription() const override {
        return "The size of k-mers in the hash table (used as seeds during "
               "partitioning). Default is 10.";
    }
};

/**
 * Option for the command line arguments to turn of output of unmapped reads.
 */
class NoUnmappedOption : public Option {
  public:
    NoUnmappedOption() : Option("U", "no-unmapped", false, NONE, OUTPUT) {
    }
    void process(const std::string& arg, Parameters& params) const override {
        params.noUnmappedRecord = true;
    }

    std::string getDescription() const override {
        return "Do not output unmapped reads.";
    }
};

/**
 * Option for the command line arguments to turn on XA tags for multiple
 * alignments.
 */
class XATagOption : public Option {
  public:
    XATagOption() : Option("T", "XA-tag", false, NONE, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.XATag = true;
    }

    std::string getDescription() const override {
        return "Output secondary alignments in XA tag for SAM format.";
    }
};

/**
 * Option for the command line arguments to trigger the help.
 */
class HelpOption : public Option {
  public:
    HelpOption() : Option("h", "help", false, NONE, HELP) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        Parameters::printHelp();
        exit(0);
    }

    std::string getDescription() const override {
        return "Print this help message.";
    }
};

const std::vector<std::shared_ptr<Option>> Parameters::options = {
    std::make_shared<PartitioningOption>(),
    std::make_shared<MetricOption>(),
    std::make_shared<LogFileOption>(),
    std::make_shared<FirstReadsOption>(),
    std::make_shared<ReferenceOption>(),
    std::make_shared<OutputFileOption>(),
    std::make_shared<KmerSizeOption>(),
    std::make_shared<MaxDistanceOption>(),
    std::make_shared<NoUnmappedOption>(),

    std::make_shared<XATagOption>(),
    std::make_shared<SearchSchemeOption>(),
    std::make_shared<CustomOption>(),
    std::make_shared<DynamicSelectionOption>(),
    std::make_shared<HelpOption>()
};

void addOption(Option* option,
               std::unordered_map<std::string, Option*>& options) {
    // assert not in map yet
    assert(options.find(option->getShortOpt()) == options.end());
    assert(options.find(option->getLongOpt()) == options.end());
    options[option->getShortOpt()] = option;
    options[option->getLongOpt()] = option;
}

Parameters Parameters::processOptionalArguments(int argc, char** argv) {

    Parameters params;
    params.command = getCommand(argc, argv);

    std::unordered_map<std::string, Option*> options;
    std::unordered_set<Option*> requiredOptions;
    for (const auto& option : Parameters::options) {
        addOption(option.get(), options);
        if (option->isRequired()) {
            requiredOptions.insert(option.get());
        }
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.size() < 2 || arg[0] != '-') {
            logger.logWarning("Invalid option: " + arg);
            continue;
        }

        std::string optionName = arg[1] == '-' ? arg.substr(2) : arg.substr(1);
        auto optionIt = options.find(optionName);
        if (optionIt == options.end()) {
            logger.logWarning("Invalid option: " + optionName);
            continue;
        }

        Option* option = optionIt->second;

        // If the option requires an argument
        if (option->argumentRequired()) {
            if (i + 1 >= argc) {
                logger.logWarning("Missing argument for option: " + arg + ". " +
                                  option->ignoreMessage());
                continue;
            }
            option->process(argv[++i], params);
        } else {
            // If the option does not require an argument but can take one
            // optionally
            if (option->hasArg() && i + 1 < argc && argv[i + 1][0] != '-') {
                option->process(argv[++i], params);
            } else {
                option->process("", params);
            }
        }

        if (option->isRequired() &&
            requiredOptions.find(option) != requiredOptions.end()) {
            requiredOptions.erase(option);
        }
    }

    if (!requiredOptions.empty()) {
        std::string missingOptions = "";
        for (const auto& o : requiredOptions) {
            missingOptions +=
                "\n\t" + o->getOptionTag() + " " + o->getOptionDescription();
        }
        throw std::runtime_error("Missing required options:" + missingOptions);
    }

    if (params.maxDistance > 4 && params.custom.empty() &&
        params.dynamicSelectionPath.empty() && params.searchScheme != "naive") {
        if (params.searchScheme != "minU" && params.searchScheme != "columba") {

            throw std::runtime_error(
                "Hard-coded search schemes are only available for "
                "up to 4 errors. Use a custom search scheme instead.");
        } else if (params.searchScheme == "minU" && params.maxDistance > 7) {
            throw std::runtime_error(
                "MinU search scheme is only available for up to 7 errors.");
        }
    }

    if (params.metric == EDIT && params.maxDistance > MAX_K_EDIT) {
        throw std::runtime_error("Edit distance only supports up to " +
                                 std::to_string(MAX_K_EDIT) + " errors");
    }

    if (params.maxDistance != 4 && params.searchScheme == "manbest") {
        throw std::runtime_error("manbest only supports 4 allowed errors");
    }

    if (params.maxDistance > MAX_K) {
        throw std::runtime_error("Max distance cannot be higher than " +
                                 std::to_string(MAX_K));
    }

    return params;
}

void printOption(const Option* opt, size_t tagWidth, size_t typeWidth,
                 size_t descWidth) {
    std::cout << std::left << std::setw(tagWidth) << opt->getOptionTag()
              << std::setw(typeWidth) << opt->getOptionArgumentType();

    // Handle description to wrap within maxWidth
    std::string description = opt->getOptionDescription();

    // remove leading spaces
    while (description.size() > 0 && description[0] == ' ') {
        description = description.substr(1);
    }

    while (description.size() > 0) {
        if (description.size() <= descWidth) {
            std::cout << description << std::endl;
            return;
        }

        // find last space before end
        size_t end = std::min(description.size(), descWidth);
        size_t lastSpace = description.find_last_of(' ', end);
        if (lastSpace == std::string::npos) {
            lastSpace = end;
        }

        // split into line and remaining
        std::string line = description.substr(0, lastSpace);
        description = description.substr(lastSpace);

        std::cout << line << std::endl;

        // remove leading spaces
        while (description.size() > 0 && description[0] == ' ') {
            description = description.substr(1);
        }
        if (description.size() > 0) {
            std::cout << std::setw(tagWidth + typeWidth) << "";
        }
    }
}

std::string optionTypeToString(OptionType type) {
    switch (type) {
    case REQUIRED:
        return "Required";
    case ALIGNMENT:
        return "Alignment";
    case OUTPUT:
        return "Output";
    case ADVANCED:
        return "Advanced";
    case HELP:
        return "Help";
    default:
        return "Unknown type";
    }
}

void Parameters::printHelp() {
    std::cout << "Usage: b-move [OPTIONS]\n\n";

    const size_t maxWidth = 100;

    // get the maximum value of sizeCol1 of all options
    size_t tagWidth = 0, typeWidth = 0;
    for (const auto& option : options) {
        tagWidth = std::max(tagWidth, option->sizeOptionTag());
        typeWidth = std::max(typeWidth, option->sizeOptionArgumentType());
    }

    tagWidth += 3;  // padding
    typeWidth += 3; // padding

    // Determine the maximum substring length that fits within
    // remaining space
    size_t descWidth = maxWidth - tagWidth - typeWidth;

    // print each option by option type
    for (uint16_t i = 0; i < OptionType::NUM_OPTION_TYPES; i++) {
        std::cout << optionTypeToString((OptionType)i) << " options:\n";
        for (const auto& opt : options) {
            if (opt->getOptionType() == (OptionType)i) {
                printOption(opt.get(), tagWidth, typeWidth, descWidth);
            }
        }
        std::cout << "\n";
    }
}

std::unique_ptr<SearchStrategy>
Parameters::createStrategy(IndexInterface& index) const {
    std::unique_ptr<SearchStrategy> strategy;

    if (!dynamicSelectionPath.empty()) {
        strategy.reset(new MultipleSchemesStrategy(
            index, dynamicSelectionPath, pStrategy, metric, mMode, sMode));
    } else if (!custom.empty()) {
        DynamicCustomStrategy* customStrategy =
            DynamicCustomStrategy::createDynamicCustomSearchStrategy(
                index, custom, pStrategy, metric, mMode, sMode);

        strategy.reset(customStrategy);

        if (mMode == ALL) {
            if (!strategy->supportsDistanceScore(maxDistance)) {
                throw std::runtime_error(
                    "Custom search scheme with path " + custom +
                    " does not support the given distance score: " +
                    std::to_string(maxDistance));
            }

        } else if (mMode == BEST) {
            int maxForBest = 0;
            if (!strategy->supportsBestMapping(maxForBest)) {
                throw std::runtime_error("Custom search scheme with path " +
                                         custom +
                                         " does not support best mapping");
            }
            if (maxForBest < MAX_K) {
                logger.logWarning(
                    "Custom search scheme with path " + custom +
                    " does not support best mapping for more than " +
                    std::to_string(maxForBest) + " errors");
            }
        }

    } else {
        // Handle other search schemes
        if (searchScheme == "kuch1") {
            strategy.reset(
                new KucherovKPlus1(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "kuch2") {
            strategy.reset(
                new KucherovKPlus2(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "kianfar") {
            strategy.reset(
                new OptimalKianfar(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "01*0") {
            strategy.reset(new O1StarSearchStrategy(index, pStrategy, metric,
                                                    mMode, sMode));
        } else if (searchScheme == "pigeon") {
            strategy.reset(new PigeonHoleSearchStrategy(index, pStrategy,
                                                        metric, mMode, sMode));
        } else if (searchScheme == "naive") {
            strategy.reset(new NaiveBackTrackingStrategy(index, pStrategy,
                                                         metric, mMode, sMode));
        } else if (searchScheme == "minU") {
            strategy.reset(
                new MinUSearchStrategy(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "columba") {
            DynamicColumbaStrategy* columbaStrategy =
                DynamicColumbaStrategy::createDynamicColumbaStrategy(
                    index, pStrategy, metric, mMode, sMode);

            // Assuming DynamicCustomStrategy is a subclass of SearchStrategy
            strategy.reset(columbaStrategy);
        } else {
            // should not get here
            throw std::runtime_error(searchScheme +
                                     " is not an option as search scheme");
        }
    }

    // Set common parameters for all strategies
    strategy->setUnmappedSam(!noUnmappedRecord);
    strategy->setXATag(XATag);
    strategy->setSamOutput(outputIsSAM);

    return strategy;
}
