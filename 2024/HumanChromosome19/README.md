# Human Chromosome 19 Data

## References

We utilized 512 haplotypes of human chromosome 19 from the [1000 Genomes Project](https://www.internationalgenome.org/data-portal/data-collection/phase-3). The data can also be accessed directly [here](http://dolomit.cs.tu-dortmund.de/tudocomp/chr19.1000.fa.xz). This file was originally created for [PHONI: Streamed Matching Statistics with Multi-Genome References](https://doi.org/10.1109/DCC50243.2021.00027). We selected the first 512 haplotypes, generating the list of sample identifiers available [here](SampleIDs.txt).

To create a pan-genome with X number of human chromosome 19 haplotypes, we selected the first X haplotypes from this file.

## Reads

We randomly sampled 100,000 single-end 151bp reads from the read set with accession [SRR17981962](https://www.ebi.ac.uk/ena/browser/view/SRR17981962). The samples are accessible [here](SRR17981962_sampled.fastq).
