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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "definitions.h" // for length_t, BEST, DYNAMIC

#include <memory>   // for shared_ptr
#include <stddef.h> // for size_t
#include <string>   // for string
#include <vector>   // for vector

class IndexInterface;
class SearchStrategy;

// The hardcoded search schemes that can be used
extern const std::vector<std::string> schemes;

// forward declaration
struct Parameters;

/**
 * Enum to represent the different types of arguments
 */
enum ArgumentType { STRING, INTEGER, NONE };

/**
 * Enum to represent the different types of options
 */
enum OptionType {
    REQUIRED,
    ALIGNMENT,
    OUTPUT,
    ADVANCED,
    HELP,
    NUM_OPTION_TYPES // helper to get the number of option types
};

/**
 * Function to convert an option type to a string for the help function
 */
std::string optionTypeToString(OptionType type);

/**
 * A class to represent a command line option.
 */
class Option {

  protected:
    std::string shortOpt; ///<  The short flag for the option
    std::string longOpt;  ///< The long flag for the option

    bool hasArgument;          ///< If the option has an argument
    ArgumentType argumentType; ///< The type of the argument (NONE if no)

    OptionType optionType; ///< The type of the option

    /**
     * @brief Get the description of the option
     */
    virtual std::string getDescription() const = 0;

    /**
     * @brief Convert an argument type to a string
     */
    static std::string argumentTypeToString(ArgumentType argumentType) {
        switch (argumentType) {
        case STRING:
            return "STR";
        case INTEGER:
            return "INT";
        case NONE:
            return "   ";
        default:
            return "   ";
        }
    }

  public:
    /**
     * @brief Process the argument of the option
     *
     * @param arg the argument
     * @param params the parameters struct, if the argument is valid the
     * relevant parameters are updated
     */
    virtual void process(const std::string& arg, Parameters& params) const = 0;

    /**
     * @brief Construct a new Option object
     *
     * @param shortOpt the short flag for the option
     * @param longOpt the long flag for the option
     * @param hasArgument if the option has an argument
     * @param argumentType the type of the argument (NONE if no argument)
     * @param optionType the type of the option
     */
    Option(const std::string& shortOpt, const std::string& longOpt,
           bool hasArgument, ArgumentType argumentType, OptionType optionType)
        : shortOpt(shortOpt), longOpt(longOpt), hasArgument(hasArgument),
          argumentType(argumentType), optionType(optionType) {
    }

    virtual ~Option() {
    } // Virtual destructor

    /**
     * @brief Check if the option has an argument
     * @return true if the option has an argument
     */
    bool hasArg() const {
        return hasArgument;
    }

    /**
     * @brief Check if the option requires an argument
     * @return true if the option requires an argument
     */
    virtual bool argumentRequired() const {
        return hasArgument;
    }

    /**
     * @brief return the message to print when the option is ignored due to
     * invalid configuration.
     * @return the message to print when the option is ignored
     */
    std::string ignoreMessage() const {
        return "Ignoring -" + shortOpt + "/--" + longOpt + " option";
    }

    /**
     * @brief Get the short flag for the option
     * @return the short flag for the option
     */
    std::string getShortOpt() const {
        return shortOpt;
    }

    /**
     * @brief Get the long flag for the option
     * @return the long flag for the option
     */
    std::string getLongOpt() const {
        return longOpt;
    }

    /**
     * @brief Check if the option is required
     * @return true if the option is required
     */
    bool isRequired() const {
        return optionType == REQUIRED;
    }

    /**
     * @brief get the size of the option tag (for pretty printing)
     */
    size_t sizeOptionTag() const {
        return 1 + shortOpt.size() + 4 + longOpt.size();
    }

    /**
     * @brief get the size of the option argument type (for pretty printing)
     */
    size_t sizeOptionArgumentType() const {
        return argumentTypeToString(argumentType).size();
    }

    /**
     * @brief get the size of the option description (for pretty printing)
     */
    size_t sizeOptionDescription() const {
        return getDescription().size();
    }

    /**
     * @brief get option tag (for pretty printing)
     */
    std::string getOptionTag() const {
        return "  -" + shortOpt + ", --" + longOpt;
    }

    /**
     * @brief get option argument type (for pretty printing)
     */
    std::string getOptionArgumentType() const {
        return argumentTypeToString(argumentType);
    }

    /**
     * @brief get option description (for pretty printing)
     */
    std::string getOptionDescription() const {
        return getDescription();
    }

    /**
     * @brief get the option type
     */
    OptionType getOptionType() const {
        return optionType;
    }
};

/**
 * Struct Parameters contains all the parameters that can be set by the user
 * Provides a method to process the optional arguments.
 * Provides a method to create a search strategy based on the parameters.
 * Provides a method to get the maximum distance or identity based on the mode.
 */
struct Parameters {

    // Parameters for the search strategy mode
    PartitionStrategy pStrategy = DYNAMIC; // default dynamic partitioning
    DistanceMetric metric = EDIT;          // default edit distance (optimized)
    MappingMode mMode = ALL;              // Default mode is best-map
    SequencingMode sMode = SINGLE_END;     // Default mode is single end reads

    // Parameters for the search scheme to use
    std::string searchScheme = "columba"; // default search scheme
    std::string custom =
        ""; // path to custom search scheme (overrides default search scheme)
    std::string dynamicSelectionPath =
        ""; // path to custom search scheme with dynamic selection (overrides
            // default search scheme)

    // Input output parameters
    std::string firstReadsFile = "";  // path to first reads file
    std::string secondReadsFile = ""; // path to second reads file (optional)
    std::string base = "";            // path to the basename of the index
    std::string outputFile = "bmoveOutput.sam"; // path to the output file
    bool outputIsSAM = true;                      // output in SAM format
    std::string logFile = "";                     // Path to the log file
    std::string command = ""; // The command that was used to run the program
    bool noUnmappedRecord = false; // Do not output the unmapped reads
    bool XATag = false;            // Add secondary alignments via XA tag

    // Numerical parameters
    length_t maxDistance = 0; // the maximum allowed distance (for ALL mode)
    length_t kmerSize = 10;   // The size of k-mers in the hash table (used as
                              // seeds during partitioning)

    /**
     * Get the command that was used to run the program
     */
    static std::string getCommand(int argc, char** argv) {
        std::string command;
        for (int i = 0; i < argc; ++i) {
            command += argv[i];
            if (i < argc - 1) {
                command += " ";
            }
        }
        return command;
    }

    /**
     * Process the optional arguments given by the user
     * @param argc The number of arguments
     * @param argv The arguments
     * @return The parameters struct
     */
    static Parameters processOptionalArguments(int argc, char** argv);

    /**
     * Create a search strategy based on the parameters
     * @param index Pointer to he index to search in
     * @return The search strategy
     */
    std::unique_ptr<SearchStrategy> createStrategy(IndexInterface& index) const;

    /**
     * Get the maximum distance or identity based on the mode
     * @return The maximum distance or minimum identity
     */
    length_t getMaxOrIdentity() const {
        return maxDistance;
    }

    static void printHelp();

  private:
    /// @brief  The options that can be set by the user
    const static std::vector<std::shared_ptr<Option>> options;
};

#endif // PARAMETERS_H