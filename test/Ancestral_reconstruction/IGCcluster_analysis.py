# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 22:00:31 2017

@author: xjw1001001
"""

import numpy as np
from IGCexpansion.CodonGeneconFunc import *

llpath = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/'
path = llpath + 'matrix/sitewise_IGC_statmatrix/'
paralog_list = [['YLR406C', 'YDL075W'],#pair#TODO: other data
 ['YER131W', 'YGL189C'], ['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'], ['YLR333C', 'YGR027C'], ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'], ['EDN', 'ECP'],['ERa', 'ERb'],['ARa', 'ERa'],['AR', 'MR'],['AR', 'GR'],['AR', 'PR'],
 ['MR', 'GR'],['MR', 'PR'],['PR', 'GR'] ] 
ARMRGRPR_list = [['AR', 'MR'],['AR', 'GR'],['AR', 'PR'],['MR', 'GR'],['MR', 'PR'],['PR', 'GR']]
Yeast_list = [['YLR406C', 'YDL075W'], ['YER131W', 'YGL189C'],['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'], ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'],
 ['YLR333C', 'YGR027C'],['YMR142C', 'YDL082W'], ['YER102W', 'YBL072C']]

EDNECP_newicktree ='/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
Yeast_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/YeastTree.newick'
ERa_ERb_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/ER/species.newick'
ARa_ERa_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/ARa_ERa/ERa_ARa_species.newick'
ARMRGRPR_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/AR_MR_GR_PR/species_common/species_common.newick'

bases = 'tcag'.upper()
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table    = dict(zip(codons, amino_acids))
codon_nonstop  = [a for a in codon_table.keys() if not codon_table[a]=='*']
codon_to_state = {a.upper() : i for (i, a) in enumerate(codon_nonstop)}
state_to_codon = {i : a.upper() for (i, a) in enumerate(codon_nonstop)}
state_to_codon[61] = 'xxx'
pair_to_state  = {pair:i for i, pair in enumerate(product(codon_nonstop, repeat = 2))}
Yeast_node = 13
EDNECP_node = 25
Model_list = ['IGC','tau=0']#model
#'arg_' +'_'.join(pair) + '_MG94_' + model + '.npy'
def state_to_compositecodons(state):
    state_1, state_2 = divmod(state, 62)
    state_1 = int(state_1)
    state_2 = int(state_2)
    return (state_to_codon[state_1],state_to_codon[state_2])


#read data of posterior etc
Expected_tau = {} #paralog array
ExpectedIGC = {}
ExpectedIGC['num'] = {}
ExpectedIGC['1to2'] = {}
ExpectedIGC['2to1'] = {}

model = {}

posterior = {}
posterior['1to2'] = {}
posterior['2to1'] = {}
posterior['IGC'] = {}
ExpectedIGC['point'] = {}
ExpectedIGC['proportion'] = {}
#'_'.join(pair)
for pair in paralog_list:
    model['_'.join(pair)] = {}
    model['_'.join(pair)]['IGC'] = np.loadtxt(open(llpath + 'model_likelihood/ancestral_reconstruction_' + '_'.join(pair) + '_MG94_IGCBFGS.txt','r'))
    model['_'.join(pair)]['tau=0'] = np.loadtxt(open(llpath + 'model_likelihood/ancestral_reconstruction_' + '_'.join(pair) + '_MG94_tau=0BFGS.txt','r'))
    model['_'.join(pair)]['PAML'] = np.loadtxt(open(llpath + 'PAML/output/summary/' + '_'.join(pair) + '.txt','r'))
    Expected_tau['_'.join(pair)] = np.loadtxt(open(path + 'Expected_tau/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['num']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['1to2']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum1_2/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['2to1']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum2_1/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['point']['_'.join(pair)] = np.loadtxt(open(path + 'SitewiseExpectedpointMutation/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['proportion']['_'.join(pair)] = np.loadtxt(open(path + 'Sitewiseporpotion/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    posterior['1to2']['_'.join(pair)] = np.loadtxt(open(path + 'posterior/' + '_'.join(pair) + '_MG94_IGC_1to2.txt','r'))
    posterior['2to1']['_'.join(pair)] = np.loadtxt(open(path + 'posterior/' + '_'.join(pair) + '_MG94_IGC_2to1.txt','r'))
    posterior['IGC']['_'.join(pair)] = np.loadtxt(open(path + 'posterior/' + '_'.join(pair) + '_MG94_IGC_IGC.txt','r'))

#generate 0 1 x sequence for each data
reconstruct_path = 'matrix/reconstruction_likelihood/npy/'
dict_all = {}
difference_threshold_begin = 0.5 #threshold for difference
difference_threshold_end   = 0.5
point_mutation_threshold   = 0.2
IGC_high_threshold = 0.5
IGC_low_threshold  = 0.1

for pair in paralog_list:
    # read data
    dict_all['_'.join(pair)]={}
    for model in Model_list:
        dict_all['_'.join(pair)][model]={}
        dict_all['_'.join(pair)][model]['arg'] = np.load(llpath+reconstruct_path+'arg_' +'_'.join(pair) + '_MG94_' + model + '.npy')
        dict_all['_'.join(pair)][model]['likelihood'] = np.load(llpath+reconstruct_path+'likelihood_' +'_'.join(pair) + '_MG94_' + model + '.npy')

branchwise_information = {}#1.how about begin difference near  0.5
branchwise_assign_1to2 = {}
branchwise_assign_2to1 = {}
branchwise_assign_IGC = {}
branchwise_display = {}
##Yeast
plist = Yeast_list
tree = Yeast_newicktree
outgroup = 'kluyveri'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
edge_to_num = {edge_list[i]:i for i in range(len(edge_list))}


for pair in plist:
    branchwise_information['_'.join(pair)] = {}
    branchwise_assign_1to2['_'.join(pair)] = {}
    branchwise_assign_2to1['_'.join(pair)] = {}
    branchwise_assign_IGC['_'.join(pair)] = {}
    branchwise_display['_'.join(pair)] = {}
    filename = open(llpath+ 'cluster_result/' + '_'.join(pair)  + '.txt' ,'w')
    for branch in edge_list:
        if branch[1] == outgroup:
            continue
        printflag = 0
        branchwise_display['_'.join(pair)][branch] = [0 for site in range(len(posterior['1to2']['_'.join(pair)]))]
        branchwise_information['_'.join(pair)][branch] = []
        branchwise_assign_1to2['_'.join(pair)][branch] = ''
        branchwise_assign_2to1['_'.join(pair)][branch] = ''
        branchwise_assign_IGC['_'.join(pair)][branch] = ''
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            begin_difference = 0
            end_difference = 0
            for i in range(10):#probability of first state difference
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][i])
                if state1 != state2:
                    begin_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[0]]][i]
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][i])
                if state1 != state2:
                    end_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[1]]][i]
            branchwise_information['_'.join(pair)][branch].append({})
            branchwise_information['_'.join(pair)][branch][site]['begin_difference'] = begin_difference
            branchwise_information['_'.join(pair)][branch][site]['end_difference'] = end_difference
            branchwise_information['_'.join(pair)][branch][site]['point_mutation'] = ExpectedIGC['point']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] = posterior['1to2']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] = posterior['2to1']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC']     = posterior['IGC']['_'.join(pair)][site][edge_to_num[branch]]
            
            if branchwise_information['_'.join(pair)][branch][site]['begin_difference'] < difference_threshold_begin:
                if branchwise_information['_'.join(pair)][branch][site]['end_difference'] < difference_threshold_end and branchwise_information['_'.join(pair)][branch][site]['point_mutation'] < point_mutation_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='x'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='x'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='x'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
            else:
                if branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_high_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_low_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_high_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_low_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_high_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_low_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
        
        for site in range(len(posterior['1to2']['_'.join(pair)])-5):
            flag = 0
            for i in range(5):
                if branchwise_assign_IGC['_'.join(pair)][branch][site+i] == '1':
                    flag += 1
            if flag >= 2:
                for i in range(5):
                    branchwise_display['_'.join(pair)][branch][site+i] = 1
                    printflag = 1
        
        if printflag == 0:
            continue
        
        filename.write(str(branch)+ '\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(str(site) + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1]] + '\t')
        filename.write('\n')
    filename.close()
            


##EDN
plist = [['EDN','ECP']]
tree = EDNECP_newicktree
outgroup = 'Saguinus_oedipus'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
edge_to_num = {edge_list[i]:i for i in range(len(edge_list))}


for pair in plist:
    branchwise_information['_'.join(pair)] = {}
    branchwise_assign_1to2['_'.join(pair)] = {}
    branchwise_assign_2to1['_'.join(pair)] = {}
    branchwise_assign_IGC['_'.join(pair)] = {}
    branchwise_display['_'.join(pair)] = {}
    filename = open(llpath+ 'cluster_result/' + '_'.join(pair)  + '.txt' ,'w')
    for branch in edge_list:
        if branch[1] == outgroup:
            continue
        printflag = 0
        branchwise_display['_'.join(pair)][branch] = [0 for site in range(len(posterior['1to2']['_'.join(pair)]))]
        branchwise_information['_'.join(pair)][branch] = []
        branchwise_assign_1to2['_'.join(pair)][branch] = ''
        branchwise_assign_2to1['_'.join(pair)][branch] = ''
        branchwise_assign_IGC['_'.join(pair)][branch] = ''
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            begin_difference = 0
            end_difference = 0
            for i in range(10):#probability of first state difference
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][i])
                if state1 != state2:
                    begin_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[0]]][i]
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][i])
                if state1 != state2:
                    end_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[1]]][i]
            branchwise_information['_'.join(pair)][branch].append({})
            branchwise_information['_'.join(pair)][branch][site]['begin_difference'] = begin_difference
            branchwise_information['_'.join(pair)][branch][site]['end_difference'] = end_difference
            branchwise_information['_'.join(pair)][branch][site]['point_mutation'] = ExpectedIGC['point']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] = posterior['1to2']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] = posterior['2to1']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC']     = posterior['IGC']['_'.join(pair)][site][edge_to_num[branch]]
            
            if branchwise_information['_'.join(pair)][branch][site]['begin_difference'] < difference_threshold_begin:
                if branchwise_information['_'.join(pair)][branch][site]['end_difference'] < difference_threshold_end and branchwise_information['_'.join(pair)][branch][site]['point_mutation'] < point_mutation_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='x'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='x'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='x'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
            else:
                if branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_high_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_low_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_high_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_low_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_high_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_low_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
        
        for site in range(len(posterior['1to2']['_'.join(pair)])-5):
            flag = 0
            for i in range(5):
                if branchwise_assign_IGC['_'.join(pair)][branch][site+i] == '1':
                    flag += 1
            if flag >= 2:
                for i in range(5):
                    branchwise_display['_'.join(pair)][branch][site+i] = 1
                    printflag = 1
        
        if printflag == 0:
            continue
        
        filename.write(str(branch)+ '\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(str(site) + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1]] + '\t')
        filename.write('\n')
    filename.close()
            

##ERaERb
plist = [['ERa', 'ERb']]
tree = ERa_ERb_newicktree
outgroup = 'Branchiostoma_floridae'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
edge_to_num = {edge_list[i]:i for i in range(len(edge_list))}


for pair in plist:
    branchwise_information['_'.join(pair)] = {}
    branchwise_assign_1to2['_'.join(pair)] = {}
    branchwise_assign_2to1['_'.join(pair)] = {}
    branchwise_assign_IGC['_'.join(pair)] = {}
    branchwise_display['_'.join(pair)] = {}
    filename = open(llpath+ 'cluster_result/' + '_'.join(pair)  + '.txt' ,'w')
    for branch in edge_list:
        if branch[1] == outgroup:
            continue
        printflag = 0
        branchwise_display['_'.join(pair)][branch] = [0 for site in range(len(posterior['1to2']['_'.join(pair)]))]
        branchwise_information['_'.join(pair)][branch] = []
        branchwise_assign_1to2['_'.join(pair)][branch] = ''
        branchwise_assign_2to1['_'.join(pair)][branch] = ''
        branchwise_assign_IGC['_'.join(pair)][branch] = ''
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            begin_difference = 0
            end_difference = 0
            for i in range(10):#probability of first state difference
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][i])
                if state1 != state2:
                    begin_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[0]]][i]
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][i])
                if state1 != state2:
                    end_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[1]]][i]
            branchwise_information['_'.join(pair)][branch].append({})
            branchwise_information['_'.join(pair)][branch][site]['begin_difference'] = begin_difference
            branchwise_information['_'.join(pair)][branch][site]['end_difference'] = end_difference
            branchwise_information['_'.join(pair)][branch][site]['point_mutation'] = ExpectedIGC['point']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] = posterior['1to2']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] = posterior['2to1']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC']     = posterior['IGC']['_'.join(pair)][site][edge_to_num[branch]]
            
            if branchwise_information['_'.join(pair)][branch][site]['begin_difference'] < difference_threshold_begin:
                if branchwise_information['_'.join(pair)][branch][site]['end_difference'] < difference_threshold_end and branchwise_information['_'.join(pair)][branch][site]['point_mutation'] < point_mutation_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='x'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='x'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='x'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
            else:
                if branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_high_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_low_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_high_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_low_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_high_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_low_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
        
        for site in range(len(posterior['1to2']['_'.join(pair)])-5):
            flag = 0
            for i in range(5):
                if branchwise_assign_IGC['_'.join(pair)][branch][site+i] == '1':
                    flag += 1
            if flag >= 2:
                for i in range(5):
                    branchwise_display['_'.join(pair)][branch][site+i] = 1
                    printflag = 1
        
        if printflag == 0:
            continue
        
        filename.write(str(branch)+ '\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(str(site) + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1]] + '\t')
        filename.write('\n')
    filename.close()
    
##ARaERa
plist = [['ARa', 'ERa']]
tree = ARa_ERa_newicktree
outgroup = 'Mus_musculus'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
edge_to_num = {edge_list[i]:i for i in range(len(edge_list))}


for pair in plist:
    branchwise_information['_'.join(pair)] = {}
    branchwise_assign_1to2['_'.join(pair)] = {}
    branchwise_assign_2to1['_'.join(pair)] = {}
    branchwise_assign_IGC['_'.join(pair)] = {}
    branchwise_display['_'.join(pair)] = {}
    filename = open(llpath+ 'cluster_result/' + '_'.join(pair)  + '.txt' ,'w')
    for branch in edge_list:
        if branch[1] == outgroup:
            continue
        printflag = 0
        branchwise_display['_'.join(pair)][branch] = [0 for site in range(len(posterior['1to2']['_'.join(pair)]))]
        branchwise_information['_'.join(pair)][branch] = []
        branchwise_assign_1to2['_'.join(pair)][branch] = ''
        branchwise_assign_2to1['_'.join(pair)][branch] = ''
        branchwise_assign_IGC['_'.join(pair)][branch] = ''
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            begin_difference = 0
            end_difference = 0
            for i in range(10):#probability of first state difference
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][i])
                if state1 != state2:
                    begin_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[0]]][i]
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][i])
                if state1 != state2:
                    end_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[1]]][i]
            branchwise_information['_'.join(pair)][branch].append({})
            branchwise_information['_'.join(pair)][branch][site]['begin_difference'] = begin_difference
            branchwise_information['_'.join(pair)][branch][site]['end_difference'] = end_difference
            branchwise_information['_'.join(pair)][branch][site]['point_mutation'] = ExpectedIGC['point']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] = posterior['1to2']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] = posterior['2to1']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC']     = posterior['IGC']['_'.join(pair)][site][edge_to_num[branch]]
            
            if branchwise_information['_'.join(pair)][branch][site]['begin_difference'] < difference_threshold_begin:
                if branchwise_information['_'.join(pair)][branch][site]['end_difference'] < difference_threshold_end and branchwise_information['_'.join(pair)][branch][site]['point_mutation'] < point_mutation_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='x'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='x'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='x'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
            else:
                if branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_high_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_low_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_high_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_low_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_high_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_low_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
        
        for site in range(len(posterior['1to2']['_'.join(pair)])-5):
            flag = 0
            for i in range(5):
                if branchwise_assign_IGC['_'.join(pair)][branch][site+i] == '1':
                    flag += 1
            if flag >= 2:
                for i in range(5):
                    branchwise_display['_'.join(pair)][branch][site+i] = 1
                    printflag = 1
        
        if printflag == 0:
            continue
        
        filename.write(str(branch)+ '\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(str(site) + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1]] + '\t')
        filename.write('\n')
    filename.close()
    
##ARMRGRPR
plist =ARMRGRPR_list
tree = ARMRGRPR_newicktree
outgroup = 'Aplysia_californica'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
edge_to_num = {edge_list[i]:i for i in range(len(edge_list))}


for pair in plist:
    branchwise_information['_'.join(pair)] = {}
    branchwise_assign_1to2['_'.join(pair)] = {}
    branchwise_assign_2to1['_'.join(pair)] = {}
    branchwise_assign_IGC['_'.join(pair)] = {}
    branchwise_display['_'.join(pair)] = {}
    filename = open(llpath+ 'cluster_result/' + '_'.join(pair)  + '.txt' ,'w')
    for branch in edge_list:
        if branch[1] == outgroup:
            continue
        printflag = 0
        branchwise_display['_'.join(pair)][branch] = [0 for site in range(len(posterior['1to2']['_'.join(pair)]))]
        branchwise_information['_'.join(pair)][branch] = []
        branchwise_assign_1to2['_'.join(pair)][branch] = ''
        branchwise_assign_2to1['_'.join(pair)][branch] = ''
        branchwise_assign_IGC['_'.join(pair)][branch] = ''
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            begin_difference = 0
            end_difference = 0
            for i in range(10):#probability of first state difference
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][i])
                if state1 != state2:
                    begin_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[0]]][i]
                state1, state2 = state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][i])
                if state1 != state2:
                    end_difference += dict_all['_'.join(pair)]['IGC']['likelihood'][site][node_to_num[branch[1]]][i]
            branchwise_information['_'.join(pair)][branch].append({})
            branchwise_information['_'.join(pair)][branch][site]['begin_difference'] = begin_difference
            branchwise_information['_'.join(pair)][branch][site]['end_difference'] = end_difference
            branchwise_information['_'.join(pair)][branch][site]['point_mutation'] = ExpectedIGC['point']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] = posterior['1to2']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] = posterior['2to1']['_'.join(pair)][site][edge_to_num[branch]]
            branchwise_information['_'.join(pair)][branch][site]['IGC']     = posterior['IGC']['_'.join(pair)][site][edge_to_num[branch]]
            
            if branchwise_information['_'.join(pair)][branch][site]['begin_difference'] < difference_threshold_begin:
                if branchwise_information['_'.join(pair)][branch][site]['end_difference'] < difference_threshold_end and branchwise_information['_'.join(pair)][branch][site]['point_mutation'] < point_mutation_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='x'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='x'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='x'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
            else:
                if branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_high_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC1to2'] > IGC_low_threshold:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_1to2['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_high_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC2to1'] > IGC_low_threshold:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_2to1['_'.join(pair)][branch]+='0'
                if branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_high_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='1'
                elif branchwise_information['_'.join(pair)][branch][site]['IGC'] > IGC_low_threshold:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='X'
                else:
                    branchwise_assign_IGC['_'.join(pair)][branch]+='0'
        
        for site in range(len(posterior['1to2']['_'.join(pair)])-5):
            flag = 0
            for i in range(5):
                if branchwise_assign_IGC['_'.join(pair)][branch][site+i] == '1':
                    flag += 1
            if flag >= 2:
                for i in range(5):
                    branchwise_display['_'.join(pair)][branch][site+i] = 1
                    printflag = 1
        
        if printflag == 0:
            continue
        
        filename.write(str(branch)+ '\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(str(site) + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[0]]][0])[1]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[0]] + '\t')
        filename.write('\n')
        
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1] + '\t')
        filename.write('\n')
        for site in range(len(posterior['1to2']['_'.join(pair)])):
            if branchwise_display['_'.join(pair)][branch][site] == 1:
                filename.write(codon_table[state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site][node_to_num[branch[1]]][0])[1]] + '\t')
        filename.write('\n')
    filename.close()