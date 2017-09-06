# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 22:22:17 2017

@author: xjw1001001
"""
import cPickle
import os
import numpy as np
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from IGCexpansion.CodonGeneconFunc import *
import argparse

def main(args):
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
            swap_dict[species1 + '_' + species2]={}
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
            
            swap_dict[species1 + '_' + species2]['length'] = len(name_to_seq[outgroup + paralog[0]])
            swap_dict[species1 + '_' + species2]['same'] = {}
            swap_dict[species1 + '_' + species2]['same']['1A1B'] = 0
            for i in range(len(name_to_seq[outgroup + paralog[0]])):
                if name_to_seq[species1 + paralog[0]][i] == name_to_seq[species1 + paralog[1]][i]:
                    swap_dict[species1 + '_' + species2]['same']['1A1B']+=1
            
            swap_dict[species1 + '_' + species2]['same']['2A2B'] = 0
            for i in range(len(name_to_seq[outgroup + paralog[0]])):
                if name_to_seq[species2 + paralog[0]][i] == name_to_seq[species2 + paralog[1]][i]:
                    swap_dict[species1 + '_' + species2]['same']['2A2B']+=1
            
            swap_dict[species1 + '_' + species2]['same']['1A2B'] = 0
            for i in range(len(name_to_seq[outgroup + paralog[0]])):
                if name_to_seq[species1 + paralog[0]][i] == name_to_seq[species2 + paralog[1]][i]:
                    swap_dict[species1 + '_' + species2]['same']['1A2B']+=1
            
            swap_dict[species1 + '_' + species2]['same']['2A1B'] = 0
            for i in range(len(name_to_seq[outgroup + paralog[0]])):
                if name_to_seq[species2 + paralog[0]][i] == name_to_seq[species1 + paralog[1]][i]:
                    swap_dict[species1 + '_' + species2]['same']['2A1B']+=1

            
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
            if not os.path.isdir('./save/' + species1 + '_' + species2+ str(args.tau) + '/'):
                os.mkdir('./save/' + species1 + '_' + species2+ str(args.tau) + '/')

            MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/' + species1 + '_' + species2 + str(args.tau) + '/', post_dup = 'N2')
            MG94_tau.tau = args.tau
            MG94_tau.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2][str(args.tau)] = np.exp(MG94_tau.x_process)
            
            alignment_file = path + species1 + '_' + species2 + '_swap.fasta'
            if not os.path.isdir('./save/' + species1 + '_' + species2+ str(args.tau) + '_swap/'):
                os.mkdir('./save/' + species1 + '_' + species2+ str(args.tau) + '_swap/')
            MG94_taus = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/'+ species1 + '_' + species2 + str(args.tau) + '_swap/', post_dup = 'N2')
            MG94_taus.tau = args.tau
            MG94_taus.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2][str(args.tau) + 'swap'] = np.exp(MG94_taus.x_process)
    
    cPickle.dump(swap_dict,open("./result/EDNECPswap.pkl","wb")) 
    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tau', required = True, help = 'tau')
    main(parser.parse_args())
    
    #data = cPickle.load(open("./result/EDNECPswap.pkl","r"))