# -*- coding: utf-8 -*-
"""
Created on Fri Sep 01 11:05:00 2017

@author: xjw1001001
"""

import cPickle
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from IGCexpansion.CodonGeneconFunc import *

if __name__ == '__main__':
    paralog = ['ERa', 'ERb']
    Force = None
    alignment_file = './data/ERa_ERb/input_ERaERb.fasta'
    newicktree = './data/ERa_ERb/(N8,N14).newick'
    seq_dict = SeqIO.to_dict(SeqIO.parse( alignment_file, "fasta" ))
    name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
    tree, edge_list, node_to_num = read_newick(newicktree, 'N14')
    
    species1_list = ['Acanthopagrus_schlegelii','Sparus_aurata']
    species2_list = ['Oryzias_latipes','Micropterus_salmoides','Halichoeres_tenuispinis',
                     'Ictalurus_punctatus','Danio_rerio','Carassius_auratus',
                     'Xenopus_tropicalis','Coturnix_japonica','Gallus_gallus',
                     'Bos_taurus','Sus_scrofa','Mus_musculus','Homo_sapiens']
    
    path = './data/ERa_ERb/swaptest/'
    outgroup = 'Branchiostoma_floridae'
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
            if not os.path.isdir('./save/' + species1 + '_' + species2 + '/'):
                os.mkdir('./save/' + species1 + '_' + species2 + '/')

            MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/' + species1 + '_' + species2 + '/', post_dup = 'N2')
            MG94_tau.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_tau.x_process))
            
            alignment_file = path + species1 + '_' + species2 + '_swap.fasta'
            if not os.path.isdir('./save/' + species1 + '_' + species2 + '_swap/'):
                os.mkdir('./save/' + species1 + '_' + species2 + '_swap/')
            MG94_taus = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/'+ species1 + '_' + species2 + '_swap/', post_dup = 'N2')
            MG94_taus.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_taus.x_process))
    
    cPickle.dump(swap_dict,open("./result/EDNECP_swap(N8,N14).pkl","wb")) 
    
    
    paralog = ['ERa', 'ERb']
    Force = None
    alignment_file = './data/ERa_ERb/input_ERaERb.fasta'
    newicktree = './data/ERa_ERb/(N11,Mus_musculus).newick'
    seq_dict = SeqIO.to_dict(SeqIO.parse( alignment_file, "fasta" ))
    name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
    tree, edge_list, node_to_num = read_newick(newicktree, 'N10')
    
    species1_list = ['Mus_musculus']
    species2_list = ['Bos_taurus','Sus_scrofa','Homo_sapiens']
    
    path = './data/ERa_ERb/swaptest/'
    outgroup = 'Branchiostoma_floridae'
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
            if not os.path.isdir('./save/' + species1 + '_' + species2 + '/'):
                os.mkdir('./save/' + species1 + '_' + species2 + '/')

            MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/' + species1 + '_' + species2 + '/', post_dup = 'N2')
            MG94_tau.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_tau.x_process))
            
            alignment_file = path + species1 + '_' + species2 + '_swap.fasta'
            if not os.path.isdir('./save/' + species1 + '_' + species2 + '_swap/'):
                os.mkdir('./save/' + species1 + '_' + species2 + '_swap/')
            MG94_taus = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/'+ species1 + '_' + species2 + '_swap/', post_dup = 'N2')
            MG94_taus.get_mle(True, True, 0, 'BFGS')
            swap_dict[species1 + '_' + species2].append(np.exp(MG94_taus.x_process))
    
    cPickle.dump(swap_dict,open("./result/EDNECPswap.pkl","wb")) 
    
    #data = cPickle.load(open("./result/EDNECPswap.pkl","r"))