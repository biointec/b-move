# Data for AMB special issue Submission

Information on the data used for the experiments with human chromosome 19 haplotypes can be accessed [here](HumanChromosome19/).

Information on the data used for the experiments with E. coli strains can be accessed [here](EColi/).

# Experiments for AMB special issue Submission

This section contains an overview of the commands used for benchmarking on the datasets detailed above.

## Overview

The benchmarking process involves aligning sequence data using three different tools: [`Columba`](https://github.com/biointec/columba), [`br-index`](https://github.com/U-Ar/br-index), and [`b-move`](https://github.com/biointec/b-move). Each tool is executed with certain parameters to assess their performance.

## Tools and Commands

1. **Columba**:
   - Version: [v2.0-beta](https://github.com/biointec/columba/tree/columba-2.0-beta)
   - Command: `columba` (section 4.4)
   - Parameters:
     - `-a all`: report all alignments up to the maximum number of errors.
     - `-e <errors>`: Set the maximum number of errors to `<errors>`.
     - `-s 8`: Set suffix array sparseness to 8.
     - `-i 0`: No in-text verification.

2. **br-index**:
   - Commands: `bri-locate` (section 4.3, 4.4), `bri-count` (section 4.2)
   - Parameters:
     - `-m <errors>`: Set the maximum number of errors to `<errors>` (Hamming distance only).

3. **b-move**:
   - Version: [v1.2.0](https://github.com/biointec/b-move/releases/tag/v1.2.0)
   - Commands: `bmove-locate`, `bmove-locate-no-report` (section 4.4), `bmove-benchmarkPhi` (section 4.3), `bmove-benchmarkCharExt` (section 4.2)
   - Parameters:
     - `-e <errors>`: Set the maximum number of errors to `<errors>`.
     - `-p dynamic` (section 4.4) or `-p uniform` (section 4.2, 4.3): Use dynamic or uniform partitioning.
     - `-m edit` (section 4.4) or `-m hamming` (section 4.2, 4.3): Use the edit or Hamming distance metric.
     - `-K 10` (section 4.4) or `-K 1` (section 4.2, 4.3): Set the seed k-mer size to 10 or 1 (no seeds).
     - `-S columba` (section 4.4) or `-S pigeon` (section 4.2, 4.3): Use the most optimized custom search scheme or the pigeonhole principle.
     - If the output SAM file is required for `bmove-locate`: `-o <filename>.sam`

