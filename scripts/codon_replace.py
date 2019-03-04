#!/usr/bin/env python3

# script to replace codon with another codon with higher probability value (as provided in the input table) without changing the amino acid sequence
# usage: python3 scriptname codon_probability_table.txt input.fasta intron_table.txt
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

def get_nucleotide_seq_with_maxvalue(ntseq):

	dnaseq, protseq = ("", "")
	for baseposition in range(0, len(ntseq), 3):
		codon = ntseq[baseposition:baseposition + 3]
		print("codon ", codon)
		aa_codon, nucl_codon, max_value = get_max_codon_value_for_aminoacid(codon)
		dnaseq+=nucl_codon
		protseq+=aa_codon

	return dnaseq

def get_nucleotide_seq_with_introns_with_maxvalue(ntseq, intron_positions):

	dnaseq= ""
	intron_started = False
	for baseposition in range(0, len(ntseq), 3):
		print("basepositions ", baseposition)
		codon = ntseq[baseposition:baseposition + 3]
		print("codon " , codon)
		if len(intron_positions) > 0 and baseposition+1 >= intron_positions[0][0] and baseposition+1 < intron_positions[0][1]:
			intron_started = True
			print("intron started")
		else:
			intron_started = False


		if intron_started == False:
			aa_codon, nucl_codon, max_value = get_max_codon_value_for_aminoacid(codon)
			dnaseq+=nucl_codon
		else:
			dnaseq+=codon

			intron_positions.pop(0)
			print("Intron ended ", intron_positions)
			intron_started = False

	return dnaseq


def get_max_codon_value_for_aminoacid(codon):

	aminoacid = get_amino_acid_for_codon(codon)

	# get the codon probability values

	aa_values = aa_codon_values[aminoacid]
	print("aa values ", aa_values)
	max_value=0
	for key in aa_values.keys():
		if aa_values[key] > max_value:
			max_value = aa_values[key]
			nucl_codon=key

	return aminoacid, nucl_codon, max_value

def get_intron_positions(intron_data):

	intron_fh = open(intron_data)
	introns = {}
	for line in intron_fh:
		line=line.rstrip()
		linearray = line.split()
		seqid = linearray[0]
		if len(linearray) == 2:
			intron_start_end = line.split()[1]
			intron_start = int(intron_start_end.split("-")[0])
			intron_end = int(intron_start_end.split("-")[1])
		elif len(linearray) == 3:
			intron_start = int(linearray[1])
			intron_end = int(linearray[2])
		else:
			print("Intron format incorrect. The intron position input file should be 2 column data or 3 column data. If 2 column data, the format is geneid start-end; if 3 column data the format is geneid intron_start intron_end")

		if seqid in introns.keys():
			introns[seqid].append((intron_start, intron_end) )
		else:
			introns[seqid] = [(intron_start, intron_end)]

	return introns

codon_probability_table=sys.argv[1]
dna_sequence = sys.argv[2]
if len(sys.argv) == 4:
	intron_positions=sys.argv[3]
	introns_provided=True
	introns = get_intron_positions(intron_positions)
	print("introns positions ", introns)

else:
	intron_positions = False

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
	if "dna" in record.description.lower() or "scaffold" in record.description.lower() or "contig" in record.description.lower():
		dnaseq = get_nucleotide_seq_with_introns_with_maxvalue(ntseq, introns[record.id])
	else:
		dnaseq = get_nucleotide_seq_with_maxvalue(ntseq)

	print(dnaseq)
