# Data for WABI 2024 Submission

Information on the data used for the experiments with human chromosome 19 haplotypes can be accessed [here](HumanChromosome19/).

Information on the data used for the experiments with E. coli strains can be accessed [here](EColi/).

# Experiments for WABI 2024 Submission

This section contains an overview of the commands used for benchmarking on the datasets detailed above.

## Overview

The benchmarking process involves aligning sequence data using three different tools: [`Columba`](https://github.com/biointec/columba), [`br-index`](https://github.com/U-Ar/br-index), and [`b-move`](https://github.com/biointec/b-move). Each tool is executed with certain parameters to assess their performance.

## Tools and Commands

1. **Columba**:
   - Version: [v1.2](https://github.com/biointec/columba/releases/tag/v1.2)
   - Command: `columba`
   - Parameters:
     - `-e <errors>`: Set the maximum number of errors to `<errors>`.
     - `-s 8`: Set suffix array sparseness to 8.
     - `-p dynamic`: Use dynamic partitioning.
     - `-m editopt`: Use the optimized edit distance metric.
     - `-ks 10`: Set the seed k-mer size to 10.
     - `-i 0`: No in-text verification.
     - `-ss custom search_schemes/multiple_opt/individual_schemes/scheme1/`: Use the most optimized custom search scheme.

2. **br-index**:
   - Commands: `bri-locate`, `bri-count`
   - Parameters:
     - `-m <errors>`: Set the maximum number of errors to `<errors>` (Hamming distance only).

3. **b-move**:
   - Version: [v1.0.0](https://github.com/biointec/b-move/releases/tag/v1.0.0)
   - Commands: `bmove-locate`, `bmove-benchmarkCharExt`
   - Parameters:
     - `-e <errors>`: Set the maximum number of errors to `<errors>`.
     - `-p dynamic`: Use dynamic partitioning.
     - `-m edit`: Use the edit distance metric.
     - `-ks 10`: Set the seed k-mer size to 10.
     - `-ss custom search_schemes/multiple_opt/individual_schemes/scheme1/`: Use the most optimized custom search scheme.
     - For Hamming distance alignment: `-m hamming`
     - If the output SAM file is required for `bmove-locate`: `-o <filename>.sam`

