## Introduction

Nucleic codon switching to new codon coding same amino acid with max probability value. The codon probability values are given as a two column data where the first column is nucleic codon and the second column is the probability of occurrence.

If the input nucleic sequence has introns and you would like to skip the introns, the list of introns and its start and end position can be provided. No codon switching will be done in the intron region.

## Usage:

python3 scripts/codon_replace.py codon_probability_table.txt input.fasta

python3 scripts/codon_replace.py codon_probability_table.txt input.fasta intron_list.txt


## Codon probability table sample

```
UUU     15.4
UCU      7.7
UAU      6.4
UGU      5.1
UUC     51.5
UCC     16.7

```

## Intron list sample

```
seq1    10      18
seq1    25      27
```
