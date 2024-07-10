# b-move

**b-move** is the first tool to construct **bidirectional move structures** and support **lossless approximate pattern matching** against them.

The [move structure](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2021.101) [1], introduced by Nishimoto and Tabei, is a run-length compressed index structure that represents a search text in $O(r)$ space, where $r$ is the number of runs in the Burrows-Wheeler Transform (BWT) [2] of the text. Unlike the [r-index](https://github.com/nicolaprezza/r-index) [3,4], the move structure supports very cache-efficient character extensions during pattern matching.

**b-move** extends this concept by constructing a **bidirectional move structure**. This extension indexes both the forward and reverse of the string, enabling efficient bidirectional character extensions. This allows **b-move** to extend a pattern $P$ in both directions, i.e., to both $cP$ and $Pc$. This is particularly useful for lossless approximate pattern matching, where a pattern is matched to the index to find all occurrences in the text.

### Key Features

* **Efficient Bidirectional Search**: Supports fast, cache-efficient bidirectional character extensions in run-length compressed space.
* **High Performance**: Achieves bidirectional character extensions up to 8 times faster than the [br-index](https://github.com/U-Ar/br-index) (bidirectional r-index) [5], significantly narrowing the performance gap with FM-index-based tools such as [Columba](https://github.com/biointec/columba) [6,7,8].
* **Scalability**: Maintains the favorable memory characteristics of the br-index, allowing it to handle large genome collections efficiently. For example, it can store the index of over 3000 complete E. coli genomes from NCBI's RefSeq collection within the RAM of a typical laptop.
* **Lossless Approximate Pattern Matching**: Provides a robust implementation for approximate pattern matching, ensuring no loss of alignments during searches.

## Prerequisites

**b-move** requires several packages to be installed on your system:
* CMake (minimum version 3.14)
* gcc (recommended minimum version 7.5)
* [SDSL-lite](https://github.com/simongog/sdsl-lite) [9]

**b-move** internally integrates code from the following external repositories:
* [Columba](https://github.com/biointec/columba) [6,7,8]: the base of the search schemes implementation
* [libsais](https://github.com/IlyaGrebnov/libsais) (by Ilya Grebnov): for building suffix arrays
* [br-index](https://github.com/U-Ar/br-index) [5]: for building the permuted longest common prefix array during the in-memory construction option
* [full-br-index](https://github.com/U-Ar/full-br-index) (by Yuma Arakawa): for constructing the index after prefix-free parsing
* [Big-BWT](https://gitlab.com/manzai/Big-BWT) [10]: for building the Burrows-Wheeler Transform and the suffix array samples using prefix-free parsing

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
git clone --recurse-submodules https://github.com/biointec/b-move.git
cd b-move
mkdir build
cd build
cmake ..
make 
```

## Index construction

Before **b-move** can align patterns to the bidirectional move structure, you need to build the bidirectional move structure, based on the input data. Currently, we only support input data with an alphabet of length 5 (for DNA: A, C, G, T + sentinel character: $). **b-move** can build its indexes from a fasta file directly. Stretches of non-ACGT are then either removed or replaced by a random ACGT-character.

**Important:** To use **b-move** v1.1.0 or higher instead of v1.0.0, the index must be rebuilt.

### 1. In-memory

In-memory constructions requires storing the entire suffix array in RAM using [libsais](https://github.com/IlyaGrebnov/libsais). The in-memory method is fast but requires a substantial amount of RAM for large pan-genomes. For example, constructing the index for a pan-genome of 512 human chromosome 19 haplotypes takes 2 hours and 550 GB.

To build a **b-move** index in-memory, only the search text (the pan-genome) is required. Run the following command in the `build` folder. The index can be built directly from fasta files. The base filename of the input file also determines where the index will be stored:

```bash
./bmove-build <input fasta>
```

### 2. Prefix-free parsing based
The prefix-free parsing based construction option builds the BWTs and suffix array samples using [Big-BWT](https://gitlab.com/manzai/Big-BWT). This prefix-free parsing based method is slower but uses much less memory. For example, constructing the index for a pan-genome of 512 human chromosome 19 haplotypes takes 5 hours and 64 GB, suitable for most workstations.

To build a **b-move** index using prefix-free parsing, only the search text (the pan-genome) is required. Run the following command in the `build` folder. The index must be built directly from fasta files. The base filename of the input file also determines where the index will be stored:

```bash
bash build/bmove-build-pfp.sh <input fasta>
```
 
## Alignment

The following commands are available:
* `bmove-locate`: for lossless approximate pattern matching
* `bmove-benchmarkCharExt`: for benchmarking character extension during pattern matching
 
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

## Citing b-move
The **b-move** paper has been accepted at [WABI 2024](https://doi.org/10.4230/LIPIcs.WABI.2024.10), and a preprint is available on [bioRxiv](https://doi.org/10.1101/2024.05.30.596587). When referencing **b-move** in your research, please cite:

```
Lore Depuydt, Luca Renders, Simon Van de Vyver, Lennart Veys, Travis Gagie, and Jan Fostier. b-move: faster bidirectional character extensions in a run-length compressed index. In 24th International Workshop on Algorithms in Bioinformatics (WABI 2024). Leibniz International Proceedings in Informatics (LIPIcs), Volume [TBD], pp. 10:1–10:18, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2024)
https://doi.org/10.4230/LIPIcs.WABI.2024.10
```

If you use **b-move** in your research, please also cite Nishimoto and Tabei's [paper](https://doi.org/10.4230/LIPIcs.ICALP.2021.101), which introduced the move structure. Cite as:
    
```
Takaaki Nishimoto and Yasuo Tabei. Optimal-Time Queries on BWT-Runs Compressed Indexes. In 48th International Colloquium on Automata, Languages, and Programming (ICALP 2021). Leibniz International Proceedings in Informatics (LIPIcs), Volume 198, pp. 101:1-101:15, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2021)
https://doi.org/10.4230/LIPIcs.ICALP.2021.101
```

## License

**b-move** is released under the AGPL-3.0 license. See the [LICENSE](LICENSE) file for details.

## Contact

Questions and suggestions can be directed to: 
* lore.depuydt@ugent.be
* jan.fostier@ugent.be

## References
[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-Time Queries on BWT-Runs Compressed Indexes. In 48th International Colloquium on Automata, Languages, and Programming (ICALP 2021). Leibniz International Proceedings in Informatics (LIPIcs), Volume 198, pp. 101:1-101:15, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2021)

[2] Michael Burrows and David Wheeler. A Block-Sorting Lossless Data Compression Algorithm. Research Report 124, Digital Equipment Corporation Systems Research Center, 130 Lytton Avenue, Palo Alto, California 94301, May 1994.

[3] Travis Gagie, Gonzalo Navarro, and Nicola Prezza. Optimal-time text indexing in bwtruns bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, SODA 2018, New Orleans, LA, USA, January 7-10, 2018, pages 1459–1477. SIAM, 2018.

[4] Travis Gagie, Gonzalo Navarro, and Nicola Prezza. Fully Functional Suffix Trees and Optimal Text Searching in BWT-Runs Bounded Space. J. ACM, 67(1), jan 2020.

[5] Yuma Arakawa, Gonzalo Navarro, and Kunihiko Sadakane. Bi-Directional r-Indexes. In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022), volume 223 of Leibniz International Proceedings in Informatics (LIPIcs), pages 11:1–11:14, Dagstuhl, Germany, 2022. Schloss Dagstuhl – Leibniz-Zentrum für Informatik.

[6] Luca Renders, Kathleen Marchal, and Jan Fostier. Dynamic partitioning of search patterns for approximate pattern matching using search schemes. iScience, 24(7):102687, 2021.

[7] Luca Renders, Lore Depuydt, and Jan Fostier. Approximate Pattern Matching Using Search Schemes and In-Text Verification. In Bioinformatics and Biomedical Engineering, pages 419–435, Cham, 2022. Springer International Publishing. 

[8] Luca Renders, Lore Depuydt, Sven Rahmann, and Jan Fostier. Automated design of efficient search schemes for lossless approximate pattern matching. In Research in Computational Molecular Biology, pages 164–184, Cham, 2024. Springer Nature Switzerland.

[9] Gog, S., Beller, T., Moffat, A., Petri, M. (2014). From Theory to Practice: Plug and Play with Succinct Data Structures. In: Gudmundsson, J., Katajainen, J. (eds) Experimental Algorithms. SEA 2014. Lecture Notes in Computer Science, vol 8504. Springer, Cham.

[10] Christina Boucher, Travis Gagie, Alan Kuhnle, and Giovanni Manzini. Prefix-Free Parsing for Building Big BWTs. In Laxmi Parida and Esko Ukkonen, editors, 18th International Workshop on Algorithms in Bioinformatics (WABI 2018), volume 113 of Leibniz International Proceedings in Informatics (LIPIcs), pages 2:1–2:16, Dagstuhl, Germany, 2018. Schloss Dagstuhl – LeibnizZentrum für Informatik.