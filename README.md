# b-move

b-move is the first tool that constructs **bidirectional move structures** and supports **lossless approximate pattern matching** against them.

The [move structure](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2021.101) is a run-length compressed index structure that can represent a search text in O(r) space, where r is the number of runs in the Burrows-Wheeler Transform (BWT) of the text. Unlike the r-index, the move structure supports very cache-efficient character extensions during pattern matching. The bidirectional move structure is an extension of the move structure that also indexes the reverse of the string. This allows for efficient pattern matching in both directions.

## Prerequisites

b-move requires a number of packages to be installed on your system. 

Required: 
* CMake (recommended minimum version 3.10)
* gcc (recommended minimum version 7.5)
* Google's Sparsehash
* [SDSL-lite](https://github.com/simongog/sdsl-lite)

b-move internally integrates code from the following external repositories:
* [Columba](https://github.com/biointec/columba): the base of the search schemes implementation
* [radixSA64](https://github.com/mariusmni/radixSA64): for building suffix arrays
* [br-index](https://github.com/U-Ar/br-index): for building the permuted longest common prefix array

<!-- How to install these packages:

As a root, execute the following commands:

on Redhat / Fedora distributions
```bash
yum install cmake
yum install sparsehash-devel
``` 

on Ubuntu / Debian distributions
```bash
apt-get install cmake
apt-get install libsparsehash-dev
```   -->

## Installation

Clone b-move from the GitHub address

    git clone "https://github.com/biointec/b-move.git"

From this `b-move` directory, run the following commands:
```bash
mkdir build
cd build
cmake ..
make 
```

## Usage
b-move can align patterns to a bidirectional move structure. To do this, you need to build the bidirectional move structure, based on the input data. Currently, we only support input data with an alphabet of length 5 (for DNA: A, C, G, T + sentinel character: $). b-move can build its indexes from a fasta file directly. The following commands are available:
* bmove-build for building the bidirectional move structure
* bmove-locate for lossless approximate pattern matching
* bmove-benchmarkCharExt for benchmarking character extension during pattern matching


### bmove-build

To build a b-move index, only the search text (the pan-genome) is required: run the following command in the `build` folder. The index can be built immediately from fasta files. The basefile parameter also determines where the index will be stored. 

```bash
./bmove-build <base filename>
```


Details:

```
Usage: ./bmove-build <base filename>

 [options]
  -v verbose output
```
 
### bmove-locate
bmove-locate executable supports lossless approximate pattern matching using search schemes (see [Columba](https://github.com/biointec/columba)). b-move can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. Alignments are written to SAM format (`.sam`).To align your reads, use the following command: 

```bash
./bmove-locate [options] basefilename readfile.[ext]
```

Note that currently in most cases, the following search scheme option is most efficient: `-ss custom b-move/search_schemes/multiple_opt/individual_schemes/scheme1/` (see [here](https://doi.org/10.1007/978-1-0716-3989-4_11)).

Details:
```
[options]
  -e  --max-errors       maximum edit distance [default = 0]
  -p  --partitioning     Add flag to do uniform/static/dynamic 
                            partitioning [default = dynamic]
  -m  --metric           Add flag to set distance metric 
                            (edit/hamming) [default = edit]
  -ks --kmer-size        The size of the seeds for dynamic 
                            partitioning [default = 10]
  -o  --output           The name of the outputfile if writing 
                            out the matches is desired. This file 
                            must be in .sam format. [default = "", 
                            no output file for benchmarking purposes]
  -ss --search-scheme    Choose the search scheme [default = kuch1]
                            options:
                            kuch1    Kucherov k + 1
                            kuch2    Kucherov k + 2
                            kianfar  Optimal Kianfar scheme
                            manbest  Manual best improvement for 
                                     kianfar scheme (only for ed = 4)
                            pigeon   Pigeon hole scheme
                            01*0     01*0 search scheme
                            custom   Custom search scheme, the next 
                                     parameter should be a path to 
                                     the folder containing this search 
                                     scheme
                            multiple Multiple search scheme, the next 
                                     parameter should be a path to 
                                     the folder containing the different 
                                     search schemes to choose from with 
                                     dynamic selection.

[ext]
  One of the following: fq, fastq, FASTA, fasta, fa
Following input files are required:
  <base filename>.cct:      character counts table
  <base filename>.move:     move table of T
  <base filename>.rev.move: move table of the reverse of T
  <base filename>.smpf:     first-of-run samples of the SA
  <base filename>.rev.smpf: first-of-run samples of the SA of the 
                              reverse of T
  <base filename>.smpl:     last-of-run samples of the SA
  <base filename>.rev.smpl: last-of-run samples of the SA of the 
                              reverse of T
  <base filename>.prdf:     predecessor structure of first-of-run samples
  <base filename>.prdl:     predecessor structure of last-of-run samples
  <base filename>.ftr:      mapping between first-of-run predecessor 
                              structure and run indices
  <base filename>.ltr:      mapping between last-of-run predecessor 
                              structure and run indices
  <base filename>.plcp:     PLCP array
```
 
### bmove-benchmarkCharExt
bmove-benchmarkCharExt executable benchmarks bidirectional character extension during pattern matching. It outputs how many cycles one bidirectional character extension takes on average. To run an experiment, use the following command:

```bash
./bmove-benchmarkCharExt [options] basefilename readfile.[ext]
```

Details:
```
[options]
  -e  --max-errors       maximum edit distance [default = 0]
  -p  --partitioning     Add flag to do uniform/static/dynamic 
                          partitioning [default = dynamic]
  -m  --metric           Add flag to set distance metric 
                          (edit/hamming) [default = edit]
  -ks --kmer-size        The size of the seeds for dynamic 
                          partitioning [default = 10]
  -ss --search-scheme    Choose the search scheme [default = kuch1]
                          options:
                            kuch1    Kucherov k + 1
                            kuch2    Kucherov k + 2
                            kianfar  Optimal Kianfar scheme
                            manbest  Manual best improvement for 
                                     kianfar scheme (only for ed = 4)
                            pigeon   Pigeon hole scheme
                            01*0     01*0 search scheme
                            custom   Custom search scheme, the next 
                                     parameter should be a path to 
                                     the folder containing this search 
                                     scheme
                            multiple Multiple search scheme, the next 
                                     parameter should be a path to 
                                     the folder containing the different 
                                     search schemes to choose from with 
                                     dynamic selection.

[ext]
  One of the following: fq, fastq, FASTA, fasta, fa
Following input files are required:
  <base filename>.cct:      character counts table
  <base filename>.move:     move table of T
  <base filename>.rev.move: move table of the reverse of T
```

 ## Contact

Questions and suggestions can be directed to: 
* lore.depuydt@ugent.be
* jan.fostier@ugent.be
