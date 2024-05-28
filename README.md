# b-move

**b-move** is the first tool to construct **bidirectional move structures** and support **lossless approximate pattern matching** against them.

The [move structure](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2021.101), introduced by Nishimoto and Tabei, is a run-length compressed index structure that represents a search text in $O(r)$ space, where $r$ is the number of runs in the Burrows-Wheeler Transform (BWT) of the text. Unlike the [r-index](https://github.com/nicolaprezza/r-index), the move structure supports very cache-efficient character extensions during pattern matching.

**b-move** extends this concept by constructing a **bidirectional move structure**. This extension indexes both the forward and reverse of the string, enabling efficient bidirectional character extensions. This allows **b-move** to extend a pattern $P$ in both directions, i.e., to both $cP$ and $Pc$. This is particularly useful for lossless approximate pattern matching, where a pattern is matched to the index to find all occurrences in the text.

### Key Features

* **Efficient Bidirectional Search**: Supports fast, cache-efficient bidirectional character extensions in run-length compressed space.
* **High Performance**: Achieves bidirectional character extensions up to 8 times faster than the [br-index](https://github.com/U-Ar/br-index) (bidirectional r-index), significantly narrowing the performance gap with FM-index-based tools such as [Columba](https://github.com/biointec/columba).
* **Scalability**: Maintains the favorable memory characteristics of the br-index, allowing it to handle large genome collections efficiently. For example, it can store the index of over 3000 complete E. coli genomes from NCBI's RefSeq collection within the RAM of a typical laptop.
* **Lossless Approximate Pattern Matching**: Provides a robust implementation for approximate pattern matching, ensuring no loss of alignments during searches.

## Prerequisites

**b-move** requires several packages to be installed on your system:
* CMake (recommended minimum version 3.10)
* gcc (recommended minimum version 7.5)
* Google's Sparsehash
* [SDSL-lite](https://github.com/simongog/sdsl-lite)

**b-move** internally integrates code from the following external repositories:
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

To install **b-move**, follow these steps:

```bash
git clone https://github.com/biointec/b-move.git
cd b-move
mkdir build
cd build
cmake ..
make 
```

## Usage
**b-move** can align patterns to a bidirectional move structure. To do this, you need to build the bidirectional move structure, based on the input data. Currently, we only support input data with an alphabet of length 5 (for DNA: A, C, G, T + sentinel character: $). **b-move** can build its indexes from a fasta file directly. The following commands are available:
* `bmove-build`: for building the bidirectional move structure
* `bmove-locate`: for lossless approximate pattern matching
* `bmove-benchmarkCharExt`: for benchmarking character extension during pattern matching


### bmove-build

To build a **b-move** index, only the search text (the pan-genome) is required. Run the following command in the `build` folder. The index can be built directly from fasta files. The basefile parameter also determines where the index will be stored:

```bash
./bmove-build <base filename>
```

Options:

```
  -v verbose output
```
 
### bmove-locate
The `bmove-locate` executable supports lossless approximate pattern matching using search schemes (see [Columba](https://github.com/biointec/columba)). **b-move** can align single end reads in fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. Alignments are written to SAM format (`.sam`). To align your reads, use the following command: 

```bash
./bmove-locate [options] basefilename readfile.[ext]
```

Currently, the following search scheme option is most efficient: `-ss custom b-move/search_schemes/multiple_opt/individual_schemes/scheme1/` (see [here](https://doi.org/10.1007/978-1-0716-3989-4_11)).

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
The `bmove-benchmarkCharExt` executable benchmarks bidirectional character extension during pattern matching. It outputs how many cycles one bidirectional character extension takes on average. To run an experiment, use the following command:

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

## License

**b-move** is released under the AGPL-3.0 license. See the [LICENSE](LICENSE) file for details.

## Contact

Questions and suggestions can be directed to: 
* lore.depuydt@ugent.be
* jan.fostier@ugent.be
