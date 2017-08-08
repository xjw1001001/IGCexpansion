# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 21:43:37 2017

@author: xjw1001001
"""

from Bio import SeqIO
seq_dictprotein = SeqIO.to_dict(SeqIO.parse('./ER/ERpreprocess/processed alignmentER.fasta', "fasta" ))
seq_dictnucleotide = SeqIO.to_dict(SeqIO.parse( './ER/ERpreprocess/ERaERb.fasta', "fasta" ))
seq_dictprotein.keys()
seq_dictprotein.values()

bases = 'tcag'.upper()
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

codon_to_protein    = dict(zip(codons, amino_acids))
codon_nonstop  = [a for a in codon_to_protein.keys() if not codon_to_protein[a]=='*']
codon_to_state = {a.upper() : i for (i, a) in enumerate(codon_nonstop)}
state_to_codon = {i : a.upper() for (i, a) in enumerate(codon_nonstop)}

outgroup = 'Branchiostoma_floridaeERa'
species_keys=[]
Thoronton_keys = []
for keys in seq_dictprotein.keys():
    if keys[-1] != 'T' and keys != outgroup:
        species_keys.append(keys)
        Thoronton_keys.append(keys + 'T')

Thoronton_dict={}
for keys in Thoronton_keys:
    Thoronton_dict[keys] = str(seq_dictprotein[keys].seq)

nucleotide_dict={}
for keys in species_keys:
    nucleotide_dict[keys] = str(seq_dictnucleotide[keys].seq)
nucleotide_dict[outgroup] = str(seq_dictnucleotide[outgroup].seq)

protein_dict={}
alligned_nucleotide_dict={}
for keys in species_keys:
    protein_dict[keys] = str(seq_dictprotein[keys].seq)
    alligned_nucleotide_dict[keys] = ''
    alligned_nucleotide_dict[keys + 'T'] = ''
    j=0
    for i in range(len(protein_dict[keys])):
        if protein_dict[keys][i] == '-':
            alligned_nucleotide_dict[keys]+='---'
            alligned_nucleotide_dict[keys + 'T']+='---'
        else:
            alligned_nucleotide_dict[keys]+=nucleotide_dict[keys][3*j:3*j+3]
            if protein_dict[keys][i] == Thoronton_dict[keys+'T'][i]:
                alligned_nucleotide_dict[keys + 'T']+=nucleotide_dict[keys][3*j:3*j+3]
            else:
                alligned_nucleotide_dict[keys + 'T']+='---'
            j=j+1

protein_dict[outgroup] = str(seq_dictprotein[outgroup].seq)
alligned_nucleotide_dict[outgroup] = ''
j=0
for i in range(len(protein_dict[outgroup])):
    if protein_dict[outgroup][i] == '-':
        alligned_nucleotide_dict[outgroup]+='---'
    else:
        alligned_nucleotide_dict[outgroup]+=nucleotide_dict[outgroup][3*j:3*j+3]
        j=j+1



f1=open('./ER/ERpreprocess/ERnucleotide.fasta','w+')
f2=open('./ER/ERpreprocess/ERnucleotide_Thoronton.fasta','w+')
for keys in alligned_nucleotide_dict.keys():
    f1.write('>'+keys+'\n')
    f1.write(alligned_nucleotide_dict[keys])
    f1.write('\n')

for keys in Thoronton_keys:
    f2.write('>'+keys+'\n')
    f2.write(alligned_nucleotide_dict[keys])
    f2.write('\n')

f1.close()
f2.close()
            
        
            
