#!/bin/bash

# Directory to store log files
log_dir="Chrom19/logs/align"
output_dir="${log_dir}/results"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir" || handle_error "Failed to create output directory"

# Function to handle errors
handle_error() {
    local error_message="$1"
    echo "Error: $error_message" >&2
    exit 1
}

# Check if run_param is passed as an argument, otherwise default to "10"
run_param="${1:-10}"

# List of X and ED values
X_values=(2 4 8 16 32 64 128 256 512)
ED_values=(0 1 2 3 4 5 6 7)

count_csv="${output_dir}/count_duration_seconds_chrom19_${run_param}.csv"
LF_csv="${output_dir}/LF_duration_seconds_chrom19_${run_param}.csv"
number_of_LF_calls_csv="${output_dir}/number_of_LF_calls_chrom19_${run_param}.csv"

> "$count_csv" && echo "ED,bmove_BP_hamming,bmove_full_hamming,rindex_hamming," > "$count_csv"
> "$LF_csv" && echo "ED,bmove_BP_hamming,bmove_full_hamming,rindex_hamming," > "$LF_csv"
> "$number_of_LF_calls_csv" && echo "ED,bmove_BP_hamming,bmove_full_hamming,rindex_hamming," > "$number_of_LF_calls_csv"

# Loop through each ED value
for ED in "${ED_values[@]}"; do
    # Define log file paths
    bmove_BP_log_BP="${log_dir}/columba_seq_benchmark_character_extensions_BP_X_512_ED_${ED}_run_${run_param}.log"
    bmove_BP_log_full="${log_dir}/columba_seq_benchmark_character_extensions_full_X_512_ED_${ED}_run_${run_param}.log"
    rindex_log="${log_dir}/original_brindex_benchmark_LF_X_512_ED_${ED}_run_${run_param}.log"

    # Append information to CSV files
    echo -n "$ED," >> "$count_csv"
    echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$bmove_BP_log_BP")," >> "$count_csv"
    echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$bmove_BP_log_full")," >> "$count_csv"
    echo -n "$(grep -oP 'Total time : \K[0-9.]+' "$rindex_log" | awk '{printf "%.6f", $1/1000}')," >> "$count_csv"

    echo -n "$ED," >> "$LF_csv"
    echo -n "$(grep -oP 'Average number of CPU cycles per LF query: \K[0-9.]+' "$bmove_BP_log_BP")," >> "$LF_csv"
    echo -n "$(grep -oP 'Average number of CPU cycles per LF query: \K[0-9.]+' "$bmove_BP_log_full")," >> "$LF_csv"
    echo -n "$(grep -oP 'Average number of CPU cycles per LF query: \K[0-9.]+' "$rindex_log")," >> "$LF_csv"

    echo -n "$ED," >> "$number_of_LF_calls_csv"
    echo -n "$(grep -oP 'Number of LF queries: \K[0-9]+' "$bmove_BP_log_BP")," >> "$number_of_LF_calls_csv"
    echo -n "$(grep -oP 'Number of LF queries: \K[0-9]+' "$bmove_BP_log_full")," >> "$number_of_LF_calls_csv"
    echo -n "$(grep -oP 'Number of LF queries: \K[0-9]+' "$rindex_log")," >> "$number_of_LF_calls_csv"

    # Add newline at the end of each line
    echo "" >> "$count_csv"
    echo "" >> "$LF_csv"
    echo "" >> "$number_of_LF_calls_csv"

    # Initialize CSV files for each metric
    duration_csv="${output_dir}/duration_seconds_ED_${ED}_chrom19_${run_param}.csv"
    resident_set_size_csv="${output_dir}/resident_set_size_kbytes_ED_${ED}_chrom19_${run_param}.csv"
    unique_matches_csv="${output_dir}/unique_matches_ED_${ED}_chrom19_${run_param}.csv"
    elapsed_time_csv="${output_dir}/elapsed_time_ED_${ED}_chrom19_${run_param}.csv"

    # Clear and set headers for CSV files
    > "$duration_csv" && echo "num_genomes,columba,bmove_report,bmove_OG,bmove_BP,rindex_hamming," > "$duration_csv"
    > "$resident_set_size_csv" && echo "num_genomes,columba,bmove_report,bmove_OG,bmove_BP,rindex_hamming," > "$resident_set_size_csv"
    > "$unique_matches_csv" && echo "num_genomes,columba,bmove_report,bmove_OG,bmove_BP,rindex_hamming," > "$unique_matches_csv"
    > "$elapsed_time_csv" && echo "num_genomes,columba,bmove_report,bmove_OG,bmove_BP,rindex_hamming," > "$elapsed_time_csv"

    # Loop through each X value
    for X in "${X_values[@]}"; do
        # Define log file paths
        columba_log="${log_dir}/columba_align_X_${X}_ED_${ED}_run_${run_param}.log"
        bmove_report_log="${log_dir}/bmove_alignWithReport_X_${X}_ED_${ED}_run_${run_param}.log"
        bmove_OG_log="${log_dir}/bmove_OG_align_X_${X}_ED_${ED}_run_${run_param}.log"
        bmove_BP_log="${log_dir}/bmove_BP_align_X_${X}_ED_${ED}_run_${run_param}.log"
        rindex_log="${log_dir}/original_brindex_align_X_${X}_ED_${ED}_run_${run_param}.log"

        # Append information to CSV files
        echo -n "$X," >> "$duration_csv"
        echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$columba_log")," >> "$duration_csv"
        echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$bmove_report_log")," >> "$duration_csv"
        echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$bmove_OG_log")," >> "$duration_csv"
        echo -n "$(grep -oP 'Total duration: \K[0-9.]+' "$bmove_BP_log")," >> "$duration_csv"
        echo -n "$(grep -oP 'Total time     : \K[0-9.]+' "$rindex_log" | awk '{printf "%.6f", $1/1000000}')," >> "$duration_csv"

        echo -n "$X," >> "$resident_set_size_csv"
        echo -n "$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$columba_log")," >> "$resident_set_size_csv"
        echo -n "$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$bmove_report_log")," >> "$resident_set_size_csv"
        echo -n "$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$bmove_OG_log")," >> "$resident_set_size_csv"
        echo -n "$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$bmove_BP_log")," >> "$resident_set_size_csv"
        echo -n "$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$rindex_log")," >> "$resident_set_size_csv"

        echo -n "$X," >> "$unique_matches_csv"
        echo -n "$(grep -oP 'Total no. unique matches: \K[0-9]+' "$columba_log")," >> "$unique_matches_csv"
        echo -n "$(grep -oP 'Total no. unique matches: \K[0-9]+' "$bmove_report_log")," >> "$unique_matches_csv"
        echo -n "$(grep -oP 'Total no. unique matches: \K[0-9]+' "$bmove_OG_log")," >> "$unique_matches_csv"
        echo -n "$(grep -oP 'Total no. unique matches: \K[0-9]+' "$bmove_BP_log")," >> "$unique_matches_csv"
        echo -n "$(grep -oP 'Total number of occurrences  occ = \K[0-9]+' "$rindex_log")," >> "$unique_matches_csv"

        echo -n "$X," >> "$elapsed_time_csv"
        echo -n "$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K[0-9:]+' "$columba_log")," >> "$elapsed_time_csv"
        echo -n "$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K[0-9:]+' "$bmove_report_log")," >> "$elapsed_time_csv"
        echo -n "$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K[0-9:]+' "$bmove_OG_log")," >> "$elapsed_time_csv"
        echo -n "$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K[0-9:]+' "$bmove_BP_log")," >> "$elapsed_time_csv"
        echo -n "$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K[0-9:]+' "$rindex_log")," >> "$elapsed_time_csv"

        # Add newline at the end of each line
        echo "" >> "$duration_csv"
        echo "" >> "$resident_set_size_csv"
        echo "" >> "$unique_matches_csv"
        echo "" >> "$elapsed_time_csv"
    done
done

# Completion message
echo "Script completed successfully"
