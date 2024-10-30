#!/bin/bash

# Script: bmove-build-pfp.sh
# Description: This script performs the b-move build process for PFP.
# Author: Lore Depuydt - lore.depuydt@ugent.be

# Capture start time
start_time=$(date +%s)

# We assume that this script is run from the build folder
bmove_build_exe="./bmove-build"
big_bwt_exe="./../external/Big-BWT/bigbwt"

# Default seed length
seedLength=100

# Function to show usage
showUsage() {
    echo "Usage: $0 [-l <seedLength>] <input_file>"
    echo
    echo "Optional arguments:"
    echo "  -l <seedLength>  Seed length for replacing non-ACGT characters (default: $seedLength). 0 means that no seed is used."
}

# Function to check if a file has a valid FASTA extension
# Usage: isValidFastaFile <filename>
isValidFastaFile() {
    local filename="$1"
    case "$filename" in
        *.fasta|*.fa|*.FASTA|*.FA)
            return 0  # Valid FASTA file extension
            ;;
        *)
            return 1  # Invalid FASTA file extension
            ;;
    esac
}

# Function to run a command with /usr/bin/time -v and extract time and memory usage
# Usage: runCommandWithTime <command> [<args>...]
runCommandWithTime() {
    local command="$1"
    shift
    # Run the command and capture output, while measuring time and memory usage
    (/usr/bin/time -v "$command" "$@") || {
        local status=$?
        echo "Error: Command '$command $@' failed with exit status $status." >&2
        exit $status
    }
}

# Function to parse command-line options
parseOptions() {
    # Parse command-line options
    while getopts ":l:" opt; do
        case $opt in
            l)
                seedLength=$OPTARG
                ;;
            \?)
                echo "Invalid option: -$OPTARG" >&2
                showUsage
                exit 1
                ;;
            :)
                echo "Option -$OPTARG requires an argument." >&2
                showUsage
                exit 1
                ;;
        esac
    done
    # Shift off the options and optional --
    shift $((OPTIND-1))

    # Check if input file is provided
    if [ $# -ne 1 ]; then
        showUsage
        exit 1
    fi

    input_file=$1
}

# Main script logic

# Parse command-line options
parseOptions "$@"

echo "Welcome to the b-move build process with prefix-free parsing!"
echo "-------------------------------------------------------------"
echo "Input file: $input_file"
echo "Seed length: $seedLength"
echo "-------------------------------------------------------------"

# Check if the input file is a FASTA file
if ! isValidFastaFile "$input_file"; then
    echo "Error: Input file must have a valid FASTA extension (.fasta, .fa)." >&2
    exit 1
fi

# Start the preprocessing
echo "Start preprocessing the fasta file with b-move..."
runCommandWithTime "$bmove_build_exe" --preprocess -l "$seedLength" "$input_file"
echo "Preprocessing done!"
echo "-------------------------------------------------------------"

base="${input_file%.*}"

# Start the prefix-free parsing
echo "Start prefix-free parsing for the original string..."
runCommandWithTime "$big_bwt_exe" -e -s -v "$base"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

echo "Start prefix-free parsing for the reverse string..."
runCommandWithTime "$big_bwt_exe" -e -s -v "${base}.rev"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

# Start building the b-move index
echo "Start building the b-move index..."
runCommandWithTime "$bmove_build_exe" --pfp "$base"
echo "b-move index built!"
echo "-------------------------------------------------------------"

# Remove the temporary files
echo "Remove temporary files..."
rm "${base}.bwt"
rm "${base}.rev.bwt"
rm "${base}.ssa"
rm "${base}.rev.ssa"
rm "${base}.esa"
rm "${base}.rev.esa"
rm "${base}.log"
rm "${base}.rev.log"
rm "${base}"
rm "${base}.rev"
echo "Temporary files removed!"
echo "-------------------------------------------------------------"

# Capture end time
end_time=$(date +%s)

# Calculate total elapsed time
total_time=$((end_time - start_time))

echo "The b-move build process with prefix-free parsing is finished!"
echo "Total time elapsed: $total_time seconds."
