#!/bin/bash

# Script: align_sequences.sh
# Description: This script aligns sequences using various tools on Chrom19 and EColi data.
# Usage: ./align_sequences.sh <counter>

# Custom functions
check_or_create_dir() {
    local dir="$1"
    [ -d "$dir" ] || mkdir -p "$dir" || handle_error "Failed to create directory $dir"
}

check_input_file() {
    local file="$1"
    [ -f "$file" ] || handle_error "Input file $file not found"
}

check_executable() {
    local exe="$1"
    [ -x "$exe" ] || handle_error "Executable $exe not found or not executable"
}

handle_error() {
    local error_message="$1"
    echo "Error: $error_message" >&2
    exit 1
}

execute_tool() {
    local tool="$1"
    local log_file="$2"
    local extra_args="$3"
    
    > "$log_file"
    
    /usr/bin/time -v "${executable_dir}/${tool}" $extra_args >> "$log_file" 2>&1 || handle_error "Failed executing $tool"
    
    echo "Executed $tool. Log saved to $log_file"
}

align_for_tool_and_dataset() {
    local tool="$1"
    local dataset="$2"
    local X="$3"
    local ED="$4"
    local read_file="$5"
    
    local output_dir="output_dir_${dataset}"
    local log_dir="log_dir_${dataset}"
    
    output_subdir="${!output_dir}/X_${X}_${tool}"
    input_file="$output_subdir/output_modified"
    log_file="${!log_dir}/${tool}_align_X_${X}_ED_${ED}_run_${counter}.log"
    
    local args=""
    if [ "$tool" == "columba" ]; then
        args="-s 8 -a all -e ${ED} -i 0 -r $input_file -f $read_file"
        execute_tool "${tool}_seq" "$log_file" "$args" || handle_error "Failed executing $tool for X=$X and ED=$ED"

    elif [ "$tool" == "original_brindex" ]; then
        input_file="${input_file}.txt.bri"
        args="-m ${ED} $input_file $read_file"
        execute_tool "${tool}_seq" "$log_file" "$args" || handle_error "Failed executing $tool for X=$X and ED=$ED"

    elif [ "$tool" == "bmove" ]; then 
        output_subdir="${!output_dir}/X_${X}_${tool}2"
        input_file="$output_subdir/output_modified"
        args="-a all -e ${ED} -r $input_file -f $read_file"

        log_file="${!log_dir}/${tool}_BP_align_X_${X}_ED_${ED}_run_${counter}.log"
        execute_tool "${tool}_seq_BP" "$log_file" "$args" || handle_error "Failed executing $tool BP for X=$X and ED=$ED"
        
        log_file="${!log_dir}/${tool}_OG_align_X_${X}_ED_${ED}_run_${counter}.log"
        execute_tool "${tool}_seq_OG" "$log_file" "$args" || handle_error "Failed executing $tool OG for X=$X and ED=$ED"

        log_file="${!log_dir}/${tool}_alignWithReport_X_${X}_ED_${ED}_run_${counter}.log"
        execute_tool "${tool}_seq" "$log_file" "$args" || handle_error "Failed executing $tool for X=$X and ED=$ED"
    fi
}

# Main script
output_dir_chrom19="Chrom19"
output_dir_EColi="EColi"
log_dir_chrom19="$output_dir_chrom19/logs/align"
log_dir_EColi="$output_dir_EColi/logs/align"
read_file_Chrom19="$output_dir_chrom19/reads/SRR17981962_sampled.fastq"
read_file_EColi="$output_dir_EColi/reads/SRR28249370_sampled.fastq"
executable_dir="Executables2"

counter="$1"
if [ -z "$counter" ]; then
    handle_error "Counter argument is required"
fi

check_or_create_dir "$output_dir_chrom19"
check_or_create_dir "$log_dir_chrom19"
check_or_create_dir "$output_dir_EColi"
check_or_create_dir "$log_dir_EColi"

check_input_file "$read_file_Chrom19"
check_input_file "$read_file_EColi"

for tool in "original_brindex_benchmark_LF" "original_brindex_benchmark_phi" "columba_seq_benchmark_character_extensions_BP" "columba_seq_benchmark_character_extensions_full" "columba_seq_benchmark_phi_BP" "columba_seq_benchmark_phi_full" "columba_seq_benchmark_phi_OG"; do
    check_executable "${executable_dir}/${tool}"
done

# search_scheme_opt="search_schemes/multiple_opt/individual_schemes/scheme1/"
# search_scheme_pigeon="search_schemes/pigeon/"

Chrom19_values=(2 4 8 16 32 64 128 256 512)
EColi_values=(1024 2048 3264) 
ED_values=(0 1 2 3 4)

for ED in "${ED_values[@]}"; do
    # Benchmark LF on chrom19
    output_subdir="$output_dir_chrom19/X_512_bmove2"
    input_file="$output_subdir/output_modified"
        
    log_file="$log_dir_chrom19/columba_seq_benchmark_character_extensions_BP_X_512_ED_${ED}_run_${counter}.log"
    execute_tool "columba_seq_benchmark_character_extensions_BP" "$log_file" "-r $input_file -f $read_file_Chrom19 -p uniform -a all -e ${ED} -S pigeon -m hamming -K 1" || handle_error "Failed executing columba_seq_benchmark_character_extensions_BP for X=512 and ED=$ED"
        
    log_file="$log_dir_chrom19/columba_seq_benchmark_character_extensions_full_X_512_ED_${ED}_run_${counter}.log"
    execute_tool "columba_seq_benchmark_character_extensions_full" "$log_file" "-r $input_file -f $read_file_Chrom19 -p uniform -a all -e ${ED} -S pigeon -m hamming -K 1" || handle_error "Failed executing columba_seq_benchmark_character_extensions_full for X=512 and ED=$ED"
    
    output_subdir="$output_dir_chrom19/X_512_original_brindex"
    input_file="$output_subdir/output_modified"

    log_file="$log_dir_chrom19/original_brindex_benchmark_LF_X_512_ED_${ED}_run_${counter}.log"
    execute_tool "original_brindex_benchmark_LF" "$log_file" "-m ${ED} ${input_file}.txt.bri $read_file_Chrom19" || handle_error "Failed executing original_brindex_benchmark_LF for X=512 and ED=$ED"
    
    # Benchmark phi on E. coli
    output_subdir="$output_dir_EColi/X_3264_bmove2"
    input_file="$output_subdir/output_modified"
        
    log_file="$log_dir_EColi/columba_seq_benchmark_phi_BP_X_3264_ED_${ED}_run_${counter}.log"
    execute_tool "columba_seq_benchmark_phi_BP" "$log_file" "-r $input_file -f $read_file_EColi -p uniform -a all -e ${ED} -S pigeon -m hamming -K 1" || handle_error "Failed executing columba_seq_benchmark_phi_BP for X=3264 and ED=$ED"
        
    log_file="$log_dir_EColi/columba_seq_benchmark_phi_full_X_3264_ED_${ED}_run_${counter}.log"
    execute_tool "columba_seq_benchmark_phi_full" "$log_file" "-r $input_file -f $read_file_EColi -p uniform -a all -e ${ED} -S pigeon -m hamming -K 1" || handle_error "Failed executing columba_seq_benchmark_phi_full for X=3264 and ED=$ED"
        
    log_file="$log_dir_EColi/columba_seq_benchmark_phi_OG_X_3264_ED_${ED}_run_${counter}.log"
    execute_tool "columba_seq_benchmark_phi_OG" "$log_file" "-r $input_file -f $read_file_EColi -p uniform -a all -e ${ED} -S pigeon -m hamming -K 1" || handle_error "Failed executing columba_seq_benchmark_phi_OG for X=3264 and ED=$ED"
    
    output_subdir="$output_dir_EColi/X_3264_original_brindex"
    input_file="$output_subdir/output_modified"

    log_file="$log_dir_EColi/original_brindex_benchmark_phi_X_3264_ED_${ED}_run_${counter}.log"
    execute_tool "original_brindex_benchmark_phi" "$log_file" "-m ${ED} ${input_file}.txt.bri $read_file_EColi" || handle_error "Failed executing original_brindex_benchmark_phi for X=3264 and ED=$ED"

    for X in "${Chrom19_values[@]}"; do
        for tool in "bmove" "columba" "original_brindex"; do
            align_for_tool_and_dataset "$tool" "chrom19" "$X" "$ED" "$read_file_Chrom19"
        done
        for tool in "bmove" "columba" "original_brindex"; do
            align_for_tool_and_dataset "$tool" "EColi" "$X" "$ED" "$read_file_EColi"
        done
    done
    for X in "${EColi_values[@]}"; do
        for tool in "bmove" "columba" "original_brindex"; do
            align_for_tool_and_dataset "$tool" "EColi" "$X" "$ED" "$read_file_EColi"
        done
    done
done

echo "Script completed successfully"
