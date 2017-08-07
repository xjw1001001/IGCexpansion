# -*- coding: utf-8 -*-
"""
Created on Mon Aug 07 21:43:37 2017

@author: xjw1001001
"""

from Bio import SeqIO
seq_dictThoronton = SeqIO.to_dict(SeqIO.parse( 'SR213 Final_ER_AR.fasta', "fasta" ))
seq_dictnucleotide = SeqIO.to_dict(SeqIO.parse( 'ERaERb.fasta', "fasta" ))
