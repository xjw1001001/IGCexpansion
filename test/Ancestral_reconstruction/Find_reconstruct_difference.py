# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 22:19:16 2017

@author: xjw1001001
"""

from IGCexpansion.CodonGeneconFunc import *
import matplotlib.pyplot as plt

EDNECP_newicktree ='/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
Yeast_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/YeastTree.newick'
series_path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/'

paralog_list = [['YLR406C', 'YDL075W'],#TODO: other data
 ['YER131W', 'YGL189C'],
 ['YML026C', 'YDR450W'],
 ['YNL301C', 'YOL120C'],
 ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'],
 ['YJL177W', 'YKL180W'],
 ['YBR191W', 'YPL079W'],
 ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'],
 ['YBL087C', 'YER117W'],
 ['YLR333C', 'YGR027C'],
 ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'],
 ]

bases = 'tcag'.upper()
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table    = dict(zip(codons, amino_acids))
codon_nonstop  = [a for a in codon_table.keys() if not codon_table[a]=='*']
codon_to_state = {a.upper() : i for (i, a) in enumerate(codon_nonstop)}
state_to_codon = {i : a.upper() for (i, a) in enumerate(codon_nonstop)}
pair_to_state  = {pair:i for i, pair in enumerate(product(codon_nonstop, repeat = 2))}

##IGC and PAML and IGC_0 yeast
tree = Yeast_newicktree
outgroup = 'kluyveri'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
for pair in paralog_list:
    seq_dict={}
    Model = ['IGC','IGC_0','PAML']
    
    seq_dict['IGC'] = SeqIO.to_dict(SeqIO.parse( series_path + 'series/ancestral_reconstruction_'+'_'.join(pair)+'_MG94+IGC.fasta', "fasta" ))
    seq_dict['IGC_0'] = SeqIO.to_dict(SeqIO.parse( series_path + 'series/ancestral_reconstruction_'+'_'.join(pair)+'_MG94.fasta', "fasta" ))
    seq_dict['PAML'] = SeqIO.to_dict(SeqIO.parse( series_path + 'PAMLfasta/PAML_'+'_'.join(pair)+'.fasta', "fasta" ))
    
    name_to_seq = {model:{name:str(seq_dict[model][name].seq) for name in seq_dict[model].keys()}for model in Model}
    
    nsites = len(name_to_seq['IGC'][name_to_seq['IGC'].keys()[0]])/3 #amio acids
    
    #对比不同,返回有不同位点记录
    difference_dic = []
    for site in range(nsites):
        for node in node_to_num.keys():
            if node == 'N0' or node == outgroup:
                bool_1=name_to_seq['IGC'][node+pair[0]][3*site:3*site+3]!=name_to_seq['IGC_0'][node+pair[0]][3*site:3*site+3]
                bool_2=name_to_seq['IGC'][node+pair[0]][3*site:3*site+3]!=name_to_seq['PAML'][node+pair[0]][3*site:3*site+3]
                bool_3=name_to_seq['PAML'][node+pair[0]][3*site:3*site+3]!=name_to_seq['IGC_0'][node+pair[0]][3*site:3*site+3]
                if bool_1 or bool_2 or bool_3:
                    difference_dic.append(site)
                    break
            else:
                paralog =pair[0]
                bool_1=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
                bool_2=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['PAML'][node+paralog][3*site:3*site+3]
                bool_3=name_to_seq['PAML'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
                if bool_1 or bool_2 or bool_3:
                    difference_dic.append(site)
                    break
                paralog =pair[1]
                bool_1=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
                bool_2=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['PAML'][node+paralog][3*site:3*site+3]
                bool_3=name_to_seq['PAML'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
                if bool_1 or bool_2 or bool_3:
                    difference_dic.append(site)
                    break
            
    #每个位点重建树 
    tree_dict={}
    for site in difference_dic:
        tree_dict[site] = {}
        for model in Model:
            tree_dict[site][model] = Phylo.read(tree, 'newick')
            for clade in tree_dict[site][model].get_terminals():
                if clade.name == 'N0' or clade.name == outgroup:
                    line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3]
                    clade.name = line
                else:
                    line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3] +'\n '+ name_to_seq[model][clade.name+pair[1]][3*site:3*site+3]
                    clade.name = line
            for clade in tree_dict[site][model].get_nonterminals():
                if clade.name == 'N0' or clade.name == outgroup:
                    line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3]
                    clade.name = line
                else:
                    line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3] +'\n '+ name_to_seq[model][clade.name+pair[1]][3*site:3*site+3]
                    clade.name = line
    
    
    # object-oriented plot
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    
    fig_list = []
    canvas_list = []
    for n in range(len(difference_dic)):
        fig_list.append(Figure())
        canvas_list.append(FigureCanvas(fig_list[n]))
        ax={}
        ax['tree'] = fig_list[n].add_axes([0.1,0.6,0.3,0.3])# [left, bottom, width, height] 
        ax['IGC'] = fig_list[n].add_axes([0.6,0.6,0.3,0.3])
        ax['IGC_0'] = fig_list[n].add_axes([0.1,0.1,0.3,0.3])
        ax['PAML'] = fig_list[n].add_axes([0.6,0.1,0.3,0.3])
        for key in Model:
            ax[key].set_title(key+' site:'+ str(difference_dic[n]+1),fontsize='4')
            Phylo.draw(tree_dict[difference_dic[n]][key],axes = ax[key])
        ax['tree'].set_title('tree',fontsize='4')
        Phylo.draw(Phylo.read(tree,'newick'),axes = ax['tree'])
        if not os.path.isdir('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair)):
            os.mkdir('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair))
        canvas_list[n].print_figure('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair) + '/' + str(difference_dic[n]+1) + '.jpg',dpi=800)
    


##IGC and PAML and IGC_0 EDNECP
tree = EDNECP_newicktree
outgroup = 'Saguinus_oedipus'
pair = ['EDN', 'ECP']
ktree, edge_list, node_to_num = read_newick(tree, 'N1')



seq_dict={}
Model = ['IGC','IGC_0','PAML']

seq_dict['IGC'] = SeqIO.to_dict(SeqIO.parse( series_path + 'series/ancestral_reconstruction_'+'_'.join(pair)+'_MG94+IGC.fasta', "fasta" ))
seq_dict['IGC_0'] = SeqIO.to_dict(SeqIO.parse( series_path + 'series/ancestral_reconstruction_'+'_'.join(pair)+'_MG94.fasta', "fasta" ))
seq_dict['PAML'] = SeqIO.to_dict(SeqIO.parse( series_path + 'PAMLfasta/PAML_'+'_'.join(pair)+'.fasta', "fasta" ))

name_to_seq = {model:{name:str(seq_dict[model][name].seq) for name in seq_dict[model].keys()}for model in Model}

nsites = len(name_to_seq['IGC'][name_to_seq['IGC'].keys()[0]])/3 #amio acids

#对比不同,返回有不同位点记录
difference_dic = []
for site in range(nsites):
    for node in node_to_num.keys():
        if node == 'N0' or node == outgroup:
            bool_1=name_to_seq['IGC'][node+pair[0]][3*site:3*site+3]!=name_to_seq['IGC_0'][node+pair[0]][3*site:3*site+3]
            bool_2=name_to_seq['IGC'][node+pair[0]][3*site:3*site+3]!=name_to_seq['PAML'][node+pair[0]][3*site:3*site+3]
            bool_3=name_to_seq['PAML'][node+pair[0]][3*site:3*site+3]!=name_to_seq['IGC_0'][node+pair[0]][3*site:3*site+3]
            if bool_1 or bool_2 or bool_3:
                difference_dic.append(site)
                break
        else:
            paralog =pair[0]
            bool_1=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
            bool_2=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['PAML'][node+paralog][3*site:3*site+3]
            bool_3=name_to_seq['PAML'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
            if bool_1 or bool_2 or bool_3:
                difference_dic.append(site)
                break
            paralog =pair[1]
            bool_1=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
            bool_2=name_to_seq['IGC'][node+paralog][3*site:3*site+3]!=name_to_seq['PAML'][node+paralog][3*site:3*site+3]
            bool_3=name_to_seq['PAML'][node+paralog][3*site:3*site+3]!=name_to_seq['IGC_0'][node+paralog][3*site:3*site+3]
            if bool_1 or bool_2 or bool_3:
                difference_dic.append(site)
                break
        
#每个位点重建树 EDN_ECP
tree_dict={}
for site in difference_dic:
    tree_dict[site] = {}
    for model in Model:
        tree_dict[site][model] = Phylo.read(tree, 'newick')
        for clade in tree_dict[site][model].get_terminals():
            if clade.name == 'N0' or clade.name == outgroup:
                line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3]
                clade.name = line
            else:
                line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3] +'\n '+ name_to_seq[model][clade.name+pair[1]][3*site:3*site+3]
                clade.name = line
        for clade in tree_dict[site][model].get_nonterminals():
            if clade.name == 'N0' or clade.name == outgroup:
                line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3]
                clade.name = line
            else:
                line=name_to_seq[model][clade.name+pair[0]][3*site:3*site+3] +'\n '+ name_to_seq[model][clade.name+pair[1]][3*site:3*site+3]
                clade.name = line
'''            
Phylo.draw(tree_dict[103]['IGC'])     
Phylo.draw(tree_dict[103]['IGC_0'])
Phylo.draw(tree_dict[103]['PAML'])
'''

# object-oriented plot
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

fig_list = []
canvas_list = []
for n in range(len(difference_dic)):
    fig_list.append(Figure())
    canvas_list.append(FigureCanvas(fig_list[n]))
    ax={}
    ax['tree'] = fig_list[n].add_axes([0.1,0.6,0.3,0.3])# [left, bottom, width, height] 
    ax['IGC'] = fig_list[n].add_axes([0.6,0.6,0.3,0.3])
    ax['IGC_0'] = fig_list[n].add_axes([0.1,0.1,0.3,0.3])
    ax['PAML'] = fig_list[n].add_axes([0.6,0.1,0.3,0.3])
    for key in Model:
        ax[key].set_title(key+' site:'+ str(difference_dic[n]+1),fontsize='4')
        Phylo.draw(tree_dict[difference_dic[n]][key],axes = ax[key])
    ax['tree'].set_title('tree',fontsize='4')
    Phylo.draw(Phylo.read(tree,'newick'),axes = ax['tree'])
    if not os.path.isdir('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair)):
        os.mkdir('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair))
    canvas_list[n].print_figure('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/pictures_different/' + '_'.join(pair) + '/' + str(difference_dic[n]+1) + '.jpg',dpi=800)













