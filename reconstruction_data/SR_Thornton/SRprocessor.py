# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 21:43:37 2017

@author: xjw1001001
"""

from Bio import SeqIO
seq_dictprotein = SeqIO.to_dict(SeqIO.parse( './ER/processed alignmentER.fasta', "fasta" ))
seq_dictnucleotide = SeqIO.to_dict(SeqIO.parse( './ER/ERaERb.fasta', "fasta" ))
seq_dictprotein.keys()
seq_dictprotein.values()

bases = 'tcag'.upper()
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

codon_to_protein    = dict(zip(codons, amino_acids))
codon_nonstop  = [a for a in codon_to_protein.keys() if not codon_to_protein[a]=='*']
codon_to_state = {a.upper() : i for (i, a) in enumerate(codon_nonstop)}
state_to_codon = {i : a.upper() for (i, a) in enumerate(codon_nonstop)}

for keys in seq_dictprotein.keys():
    for i in range((len(seq_dictnucleotide[keys].seq)+1)/3):
        seq_dictnucleotide[keys].seq[3*i:3*i+3]
    seq_dictprotein[keys].seq[i]
    