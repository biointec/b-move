# E. coli Data

## References

We obtained all complete E. coli genomes from the NCBI RefSeq database using the following command:

```bash
ncbi-genome-download -F fasta -l complete -g Escherichia bacteria
```

The resulting files were shuffled, generating the list of NCBI reference sequence identifiers available [here](EColi_NCBIReferenceSequenceIDs.txt).

To create a pan-genome with X number of E. coli strains, we selected the first X strains from this accession list.

## Reads

We randomly sampled 100,000 single-end 151bp reads from the read set with accession [SRR28249370](https://www.ebi.ac.uk/ena/browser/view/SRR28249370). The samples are accessible [here](SRR28249370_sampled.fastq).