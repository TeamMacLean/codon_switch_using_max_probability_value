#!/usr/bin/env python3

# script to replace codon with another codon with higher probability value (as provided in the input table) without changing the amino acid sequence

import os, sys
from Bio import SeqIO



translate_table={
					"GCA":"A",
					"GCC":"A",
					"GCG":"A",
					"GCT":"A",
					"TGC":"C",
					"TGT":"C",
					"GAC":"D",
					"GAT":"D",
					"GAA":"E",
					"GAG":"E",
					"TTC":"F",
					"TTT":"F",
					"GGA":"G",
					"GGC":"G",
					"GGG":"G",
					"GGT":"G",
					"CAC":"H",
					"CAT":"H",
					"ATA":"I",
					"ATC":"I",
					"ATT":"I",
					"AAA":"K",
					"AAG":"K",
					"CTA":"L",
					"CTC":"L",
					"CTT":"L",
					"CTG":"L",
					"TTG":"L",
					"TTA":"L",
					"ATG":"M",
					"AAC":"N",
					"AAT":"N",
					"CCA":"P",
					"CCC":"P",
					"CCG":"P",
					"CCT":"P",
					"CAA":"Q",
					"CAG":"Q",
					"AGA":"R",
					"AGG":"R",
					"CGA":"R",
					"CGC":"R",
					"CGG":"R",
					"CGT":"R",
					"AGC":"S",
					"AGT":"S",
					"TCA":"S",
					"TCC":"S",
					"TCG":"S",
					"TCT":"S",
					"TAG":"*",
					"TGA":"*",
					"TAA":"*",
					"ACA":"T",
					"ACC":"T",
					"ACG":"T",
					"ACT":"T",
					"GTA":"V",
					"GTC":"V",
					"GTG":"V",
					"GTT":"V",
					"TGG":"W",
					"TAC":"Y",
					"TAT":"Y",

                    "GCA":"A",
                    "GCC":"A",
                    "GCG":"A",
                    "GCU":"A",
                    "UGC":"C",
                    "UGU":"C",
                    "GAC":"D",
                    "GAU":"D",
                    "GAA":"E",
                    "GAG":"E",
                    "UUC":"F",
                    "UUU":"F",
                    "GGA":"G",
                    "GGC":"G",
                    "GGG":"G",
                    "GGU":"G",
                    "CAC":"H",
                    "CAU":"H",
                    "AUA":"I",
                    "AUC":"I",
                    "AUU":"I",
                    "AAA":"K",
                    "AAG":"K",
                    "CUA":"L",
                    "CUC":"L",
                    "CUU":"L",
                    "CUG":"L",
                    "UUG":"L",
                    "UUA":"L",
                    "AUG":"M",
                    "AAC":"N",
                    "AAU":"N",
                    "CCA":"P",
                    "CCC":"P",
                    "CCG":"P",
                    "CCU":"P",
                    "CAA":"Q",
                    "CAG":"Q",
                    "AGA":"R",
                    "AGG":"R",
                    "CGA":"R",
                    "CGC":"R",
                    "CGG":"R",
                    "CGU":"R",
                    "AGC":"S",
                    "AGU":"S",
                    "UCA":"S",
                    "UCC":"S",
                    "UCG":"S",
                    "UCU":"S",
                    "UAG":"*",
                    "UGA":"*",
                    "UAA":"*",
                    "ACA":"T",
                    "ACC":"T",
                    "ACG":"T",
                    "ACU":"T",
                    "GUA":"V",
                    "GUC":"V",
                    "GUG":"V",
                    "GUU":"V",
                    "UGG":"W",
                    "UAC":"Y",
                    "UAU":"Y",
					}

def get_amino_acid_for_codon(codon):
	return translate_table[codon]

def translate_with_maxvalue(ntseq):
	'''
	translates nucleotide to protein sequence
	'''
	dnaseq, protseq = ("", "")
	for baseposition in range(0, len(ntseq), 3):
		codon = ntseq[baseposition:baseposition + 3]
		print("codon ", codon)
		aa_codon, nucl_codon, max_value = get_max_codon_value_for_aminoacid(codon)
		dnaseq+=nucl_codon
		protseq+=aa_codon

	return dnaseq, protseq



def get_max_codon_value_for_aminoacid(codon):

	aminoacid = get_amino_acid_for_codon(codon)

	# get the codon probability values

	aa_values = aa_codon_values[aminoacid]
	max_value=0
	for key in aa_values.keys():
		if aa_values[key] > max_value:
			max_value = aa_values[key]
			nucl_codon=key

	return aminoacid, nucl_codon, max_value

codon_probability_table=sys.argv[1]
dna_sequence = sys.argv[2]

codon_table_fh=open(codon_probability_table)
dna_fh = open(dna_sequence)





# create codon value data structure

nucleotide_codon_values={}
for line in codon_table_fh:
	line=line.rstrip()
	codon=line.split()[0]
	value=float(line.split()[1])

	nucleotide_codon_values[codon] = value
codon_table_fh.close()

print("nucleotide codon values ", nucleotide_codon_values)
# get nucleotide_codon_values
aa_codon_values={}
for key in nucleotide_codon_values.keys():
	aa = get_amino_acid_for_codon(key)
	if aa in aa_codon_values.keys():
		aa_codon_values[aa].update({key:nucleotide_codon_values[key]})
	else:
		aa_codon_values[aa]={key:nucleotide_codon_values[key]}

print("amino acid values ", aa_codon_values)



for record in SeqIO.parse(dna_sequence, 'fasta'):
	seqid=record.description
	ntseq=str(record.seq)

	dnaseq, protseq = translate_with_maxvalue(ntseq)

	print(dnaseq, protseq)
