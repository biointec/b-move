# b-move

**b-move** is the first tool to construct the **bidirectional move structure** and support **lossless approximate pattern matching** against it.

The [move structure](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2021.101) [1], introduced by Nishimoto and Tabei, is a run-length compressed index structure that represents a search text in $O(r)$ space, where $r$ is the number of runs in the Burrows-Wheeler Transform (BWT) [2] of the text. Unlike the [r-index](https://github.com/nicolaprezza/r-index) [3,4], the move structure supports very cache-efficient character extensions, $\phi$ and $\phi^{-1}$ operations during pattern matching.

**b-move** extends this concept by constructing a **bidirectional move structure**. This extension indexes both the forward and reverse of the string, enabling efficient bidirectional character extensions. This allows **b-move** to extend a pattern $P$ in both directions, i.e., to both $cP$ and $Pc$. This is particularly useful for lossless approximate pattern matching, where a pattern is matched to the index to find all occurrences in the text.

### Key Features

* **Efficient Bidirectional Search**: Supports fast, cache-efficient bidirectional character extensions in run-length compressed space.
* **High Performance**: Achieves bidirectional character extensions as well as $\phi$ and $\phi^{-1}$ operations up to 7 times faster than the [br-index](https://github.com/U-Ar/br-index) (bidirectional r-index) [5], significantly narrowing the performance gap with FM-index-based tools such as [Columba](https://github.com/biointec/columba) [6,7,8].
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
* [{fmt}](https://github.com/fmtlib/fmt) (by Victor Zverovich): for formatting output
* [The Parallel Hashmap](https://github.com/greg7mdp/parallel-hashmap) (by Gregory Popovitch): for building the hash table that contains all wordSize-mers
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
git clone https://github.com/biointec/b-move.git
cd b-move
mkdir build
cd build
cmake ..
make 
```

### Compiler Options
This project includes configurable compiler options to adjust the how the index is queried. The following options are available:

#### LF move table bit-packing
- **Option**: `LF_MOVE_BIT_PACKED` (Default: `ON`)
- **Description**: Enables bit-packing for LF move tables to reduce memory usage. Bit-packing the move tables is always recommended, as it significantly reduces the memory footprint of the index and does not significantly affect the performance.
- **Usage**: Set to `OFF` to disable bit-packing (`-DLF_MOVE_BIT_PACKED`).

#### $\phi$ and $\phi^{-1}$ move tables
- **Option**: `PHI_MOVE` (Default: `ON`)
- **Description**: Enables move tables specifically for $\phi$ and $\phi^{-1}$ operations. This option increases the memory usage of the index but also improves the performance of $\phi$ and $\phi^{-1}$ operations. Disabling this option results in a smaller index but slower $\phi$ and $\phi^{-1}$ operations. This tradeoff is to be considered by the user.
- **Usage**: Set `PHI_MOVE` to `OFF` to disable $\phi$ and $\phi^{-1}$ move tables.
- **Sub-options**:
  - **$\phi$ and $\phi^{-1}$ move balancing**: `PHI_MOVE_BALANCED` (Default: `ON`)
    - Enables balancing the $\phi$ and $\phi^{-1}$ move tables. Balancing is always recommended, as it improves the performance of $\phi$ and $\phi^{-1}$ operations significantly.
  - **$\phi$ and $\phi^{-1}$ move bit-packing**: `PHI_MOVE_BIT_PACKED` (Default: `ON`)
    - Enables bit-packing to optimize memory use. Bit-packing the $\phi$ and $\phi^{-1}$ move tables is always recommended, as it significantly reduces the memory footprint of the index and does not significantly affect the performance.


## Index construction

Before **b-move** can align patterns to the bidirectional move structure, you need to build the bidirectional move structure, based on the input data. Currently, we only support input data with an alphabet of length 5 (for DNA: A, C, G, T + sentinel character: $). **b-move** can build its indexes from a fasta file directly. Stretches of non-ACGT characters are then replaced by a random ACGT-character.

**Important:** **b-move**'s alignment mode is not guaranteed to be compatible with indexes that were built using an older version. Consider rebuilding your index.

### 1. In-memory

In-memory construction requires storing the entire suffix array in RAM using [libsais](https://github.com/IlyaGrebnov/libsais). The in-memory method is faster but requires a substantial amount of RAM for large pan-genomes. For example, constructing the index for a pan-genome of 512 human chromosome 19 haplotypes takes 2 hours and 520 GB.

To build a **b-move** index in-memory, only the search text (the pan-genome) is required. Run the following command in the `build` folder. The index can be built directly from fasta files. The base filename of the input file also determines where the index will be stored:

```bash
./bmove-build [-l <seedLength>] [--pfp] [--preprocess] <fasta file>

Required arguments:
         fasta file: Path to the FASTA file containing the reference sequences.

Optional arguments:
        -l <seedLength> (default: 100): Seed length for replacing non-ACGT characters. Seed length 0 means that no seed is used.
```

The seed length determines how non-ACGT characters are replaced. A non-zero seed length replaces stretches of non-ACGT characters by repeating a random seed of length `seedLength`. A seed length of 0 refers to non-seeded replacement. The latter option has a negative influence on the number of runs in the BWT. The default seed length of 100 is recommended.

### 2. Prefix-free parsing based
The prefix-free parsing based construction option builds the BWTs and suffix array samples using [Big-BWT](https://gitlab.com/manzai/Big-BWT). This prefix-free parsing based method is slower but uses much less memory. For example, constructing the index for a pan-genome of 512 human chromosome 19 haplotypes takes 7 hours and 84 GB, suitable for most workstations.

To build a **b-move** index using prefix-free parsing, only the search text (the pan-genome) is required. Run the following command in the `build` folder. The index must be built directly from fasta files. The base filename of the input file also determines where the index will be stored:

```bash
bash bmove-build-pfp.sh [-l <seedLength>] <input_file>

Optional arguments:
  -l <seedLength>  Seed length for replacing non-ACGT characters (default: 100). 0 means that no seed is used.
```
 
## Alignment

The following commands are available:
* `bmove-locate`: for lossless approximate pattern matching
* `bmove-locate-no-report`: for lossless approximate pattern matching without reporting the alignments
* `bmove-benchmarkCharExt`: for benchmarking character extension during pattern matching
* `bmove-benchmarkPhi`: for benchmarking the $\phi$ and $\phi^{-1}$ operations
 
### bmove-locate
The `bmove-locate` executable supports lossless approximate pattern matching using search schemes (see [Columba](https://github.com/biointec/columba)). **b-move** can align single end reads in fasta (`.fasta`, `.fa`, `.fna`) or fastq (`.fq`, `.fastq`) format. Alignments are written to SAM format (`.sam`). To align your reads, use the following command: 

```bash
./bmove-locate [options] -r basefilename -f readfile.[ext]
```

Details:
```
Required options:
  -f, --reads-file        STR   Path to the reads file.
  -r, --reference         STR   Path to the basename of the index.

Alignment options:
  -m, --metric            STR   Distance metric to use. Options are: edit, hamming. Default is edit.
  -e, --max-distance      STR   The maximum allowed distance. Default is 0.

Output options:
  -l, --log-file          STR   Path to the log file. Default is stdout.
  -o, --output-file       STR   Path to the output file. Should be .sam or .rhs. Default is
                                bmoveOutput.sam.
  -U, --no-unmapped             Do not output unmapped reads.
  -T, --XA-tag                  Output secondary alignments in XA tag for SAM format.

Advanced options:
  -p, --partitioning      STR   Partitioning strategy to use. Options are: uniform, static, dynamic.
                                Default is dynamic.
  -K, --kmer-size         INT   The size of k-mers in the hash table (used as seeds during
                                partitioning). Default is 10.
  -S, --search-scheme     STR   Search scheme to use. Options are: kuch1, kuch2, kianfar, pigeon,
                                01*0, custom, naive, multiple, minU, columba. Default is columba.
  -c, --custom            STR   Path to custom search scheme (overrides default search scheme).
  -d, --dynamic-selection STR   Path to custom search scheme with dynamic selection (overrides
                                default search scheme).

Help options:
  -h, --help                    Print this help message.
```
See [Columba](https://github.com/biointec/columba) for more information on the options.

### bmove-locate-no-report
The `bmove-locate-no-report` executable supports lossless approximate pattern matching using search schemes in an analogous way to `bmove-locate`. The difference is that it does not report the alignments. This executable can be used to benchmark alignment performance without the overhead of writing the alignments to disk. 
 
### bmove-benchmarkCharExt
The `bmove-benchmarkCharExt` executable benchmarks bidirectional character extension during pattern matching. It outputs how many cycles one bidirectional character extension takes on average. No locating is performed. 

### bmove-benchmarkPhi
The `bmove-benchmarkPhi` executable benchmarks the $\phi$ and $\phi^{-1}$ operations during pattern matching. It outputs how many cycles one $\phi$/$\phi^{-1}$ operation takes on average. It also outputs how long one locating operation (for one SA interval) takes on average.

## Citing b-move
The **b-move** paper has been accepted at [WABI 2024](https://doi.org/10.4230/LIPIcs.WABI.2024.10). When referencing **b-move** in your research, please cite:

```
Lore Depuydt, Luca Renders, Simon Van de Vyver, Lennart Veys, Travis Gagie, and Jan Fostier. b-move: Faster Bidirectional Character Extensions in a Run-Length Compressed Index. In 24th International Workshop on Algorithms in Bioinformatics (WABI 2024). Leibniz International Proceedings in Informatics (LIPIcs), Volume 312, pp. 10:1-10:18, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2024)
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