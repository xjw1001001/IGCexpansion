# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 21:43:37 2017

@author: xjw1001001
"""

from Bio import SeqIO
seq_dictThoronton = SeqIO.to_dict(SeqIO.parse( 'SR213 Final_ER_AR.fasta', "fasta" ))
seq_dictnucleotide = SeqIO.to_dict(SeqIO.parse( 'ERaERb.fasta', "fasta" ))
seq_dictThoronton.keys()
seq_dictThoronton.values()

bases = 'tcag'.upper()
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

codon_table    = dict(zip(codons, amino_acids))
codon_nonstop  = [a for a in self.codon_table.keys() if not codon_table[a]=='*']
codon_to_state = {a.upper() : i for (i, a) in enumerate(codon_nonstop)}
state_to_codon = {i : a.upper() for (i, a) in enumerate(codon_nonstop)}

for keys in seq_dictThoronton.keys():
    for i in range(len(seq_dictThoronton[keys].seq)):
        
    seq_dictThoronton[keys].seq
    seq_dictnucleotide[keys].seq