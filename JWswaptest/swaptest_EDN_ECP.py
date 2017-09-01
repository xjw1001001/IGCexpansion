# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 22:22:17 2017

@author: xjw1001001
"""

from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from IGCexpansion.CodonGeneconFunc import *

if __name__ == '__main__':
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = './data/EDN_ECP/primateoutcome processed.fasta'
    newicktree = './data/EDN_ECP/primate_EDN_ECP_modified.newick'
    seq_dict = SeqIO.to_dict(SeqIO.parse( alignment_file, "fasta" ))
    name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
    tree, edge_list, node_to_num = read_newick(newicktree, 'N2')
    
    species1_list = ['Rhinopithecus_bieti','Colobus_angolensis']
    species2_list = ['Cercocebus_atys','Macaca_fascicularis','Macaca_mulatta',
                     'Macaca_nemestrina','Mandrillus_leucophaeus','Papio_anubis']
    
    path = './data/EDN_ECP/swaptest/'
    outgroup = 'Saguinus_oedipus'
    swap_dict = {}
    for species1 in species1_list:
        for species2 in species2_list:
            swap_dict[species1 + '_' + species2]=[]
            filename = open(path + species1 + '_' + species2 + '.fasta' ,'w')
            filename.write('>'+ outgroup + paralog[0] +'\n')
            filename.write(name_to_seq[outgroup + paralog[0]]+'\n')
            filename.write('>'+ 'species1' + paralog[0] +'\n')
            filename.write(name_to_seq[species1 + paralog[0]]+'\n')
            filename.write('>'+ 'species1' + paralog[1] +'\n')
            filename.write(name_to_seq[species1 + paralog[1]]+'\n')
            filename.write('>'+ 'species2' + paralog[0] +'\n')
            filename.write(name_to_seq[species2 + paralog[0]]+'\n')
            filename.write('>'+ 'species2' + paralog[1] +'\n')
            filename.write(name_to_seq[species2 + paralog[1]]+'\n')
            filename.close()
            

            
            filename = open(path + species1 + '_' + species2 + '_swap.fasta' ,'w')
            filename.write('>'+ outgroup + paralog[0] +'\n')
            filename.write(name_to_seq[outgroup + paralog[0]]+'\n')
            filename.write('>'+ 'species2' + paralog[0] +'\n')
            filename.write(name_to_seq[species1 + paralog[0]]+'\n')
            filename.write('>'+ 'species1' + paralog[1] +'\n')
            filename.write(name_to_seq[species1 + paralog[1]]+'\n')
            filename.write('>'+ 'species1' + paralog[0] +'\n')
            filename.write(name_to_seq[species2 + paralog[0]]+'\n')
            filename.write('>'+ 'species2' + paralog[1] +'\n')
            filename.write(name_to_seq[species2 + paralog[1]]+'\n')
            filename.close()
    
    
            alignment_file = path + species1 + '_' + species2 + '.fasta'
            MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/', post_dup = 'N2')
            MG94_tau.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_tau.x_process))
            
            alignment_file = path + species1 + '_' + species2 + '_swap.fasta'
            MG94_taus = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/', post_dup = 'N2')
            MG94_taus.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_taus.x_process))

    
            