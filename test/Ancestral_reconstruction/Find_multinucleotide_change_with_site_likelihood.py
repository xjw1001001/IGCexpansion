# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 21:25:36 2017

@author: xjw1001001
"""
from IGCexpansion.CodonGeneconFunc import *
import numpy as np

path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/'
reconstruct_path = 'matrix/reconstruction_likelihood/npy/'
paralog_list = [['YLR406C', 'YDL075W'],#pair#TODO: other data
 ['YER131W', 'YGL189C'], ['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'], ['YLR333C', 'YGR027C'], ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'], ['EDN', 'ECP'],['AR', 'MR'],['AR', 'GR'],['AR', 'PR'],['MR', 'GR'],
 ['MR', 'PR'],['PR', 'GR'] ,['ERa', 'ERb'],['ARa', 'ERa']] #


ARMRGRPR_list = [['AR', 'MR'],['AR', 'GR'],['AR', 'PR'],['MR', 'GR'],['MR', 'PR'],['PR', 'GR']]

Yeast_list = [['YLR406C', 'YDL075W'], ['YER131W', 'YGL189C'],['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'], ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'],
 ['YLR333C', 'YGR027C'],['YMR142C', 'YDL082W'], ['YER102W', 'YBL072C'],]

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

#read PAML reconstruct resule #TODO: other data
PAMLseq_dict = {}
name_to_seq = {}
for pair in paralog_list:
    PAMLseq_dict['_'.join(pair)] = SeqIO.to_dict(SeqIO.parse( '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/' + 'PAMLfasta/PAML_'+'_'.join(pair)+'.fasta', "fasta" ))
    name_to_seq['_'.join(pair)] = {name:str(PAMLseq_dict['_'.join(pair)][name].seq) for name in PAMLseq_dict['_'.join(pair)].keys()}


dict_all = {}
for pair in paralog_list:
    # read data
    dict_all['_'.join(pair)]={}
    for model in Model_list:
        dict_all['_'.join(pair)][model]={}
        dict_all['_'.join(pair)][model]['arg'] = np.load(path+reconstruct_path+'arg_' +'_'.join(pair) + '_MG94_' + model + '.npy')
        dict_all['_'.join(pair)][model]['likelihood'] = np.load(path+reconstruct_path+'likelihood_' +'_'.join(pair) + '_MG94_' + model + '.npy')

#从likelihood中找到有点arg变化,或似然变化大于0.05的位点，返回site node n， 然后建立一个关于这些位点的list
site_list_dict = {}#[pair][site][node]
for pair in paralog_list:
    site_list_dict['_'.join(pair)]=[]
    for site in range(len(dict_all['_'.join(pair)]['IGC']['arg'][:,0,0])):
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            #首先比较arg
            if dict_all['_'.join(pair)]['IGC']['arg'][site,node,0] != dict_all['_'.join(pair)]['tau=0']['arg'][site,node,0]:
                site_list_dict['_'.join(pair)].append(site)
                break
            #arg1相等,比较似然
            elif abs(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,0] - dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,0])>0.05:
                site_list_dict['_'.join(pair)].append(site)
                break
            #在上面都相等，似然相差小于0.05，比较第二大的arg，且在似然大于0.05时发现不同报告
            elif dict_all['_'.join(pair)]['IGC']['arg'][site,node,1] != dict_all['_'.join(pair)]['tau=0']['arg'][site,node,1]:
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,1]>0.05 or dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,1]>0.05:
                    site_list_dict['_'.join(pair)].append(site)
                    break

#记录这些位点所有不同的地方，存为site_information_dict，记录两种模型下所有大于0.05似然的arg和似然
site_information_dict = {}#[pair][site]{[node]:{['IGC']:[(arg1,lnl),(arg2,lnl)],['tau=0']:[]}
 #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，但tau=0中未出现多核酸替换的位点，组成 multi_site_dict 
site_multianalysis_dict = {}
multi_site_dict = {}#[pair][site]{[branch]:['ATG->CGG']}
##Yeast
plist = Yeast_list
tree = Yeast_newicktree
outgroup = 'kluyveri'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
for pair in plist:
    site_information_dict['_'.join(pair)] = {}
    site_multianalysis_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        site_information_dict['_'.join(pair)][site] = {}
        site_multianalysis_dict['_'.join(pair)][site] = {}        
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_information_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = {}
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'] = {}
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'][state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i])]=dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]
            
            for i in range(len(dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'][state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node,i])]=dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,i]
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = []
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]>0.05:
                    site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'].append((state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i]),dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]))
    
    #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，且疑似IGC，组成 multi_site_dict 
    multi_site_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        flag1 = 0
        for branch in edge_list:
            i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
            i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
            for i in range(i0):
                for j in range(i1):
                    for k in range(2):
                        count = 0
                        for l in range(3):
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                count = count + 1
                        if count >= 2:
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    flag1 = 1
            if flag1 == 1:
                break   
        if flag1 == 1:
            multi_site_dict['_'.join(pair)][site]={}
            for branch in edge_list:
                flag = 0
                i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
                i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
                for i in range(i0):
                    for j in range(i1):
                        for k in range(2):
                            count = 0
                            for l in range(3):
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                    count = count + 1
                            if count >= 2:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        flag = 1
                if flag == 1:
                    multi_site_dict['_'.join(pair)][site][branch]={}
                    for i in range(i0):
                        for j in range(i1):
                            for k in range(2):#paralog
                                count = 0
                                for l in range(3):
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                        count = count + 1
                                if count >= 2: 
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                            event = site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] + ' ' + pair[0]+' '+site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1] + ' ' + pair[1]
                                            multi_site_dict['_'.join(pair)][site][branch][event] = {}
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node2 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[1]]['IGC'][(event[5:8],event[8 + len(pair[0]) + 7:8 + len(pair[0]) + 10])]
                                            if not (event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5]) in site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'].keys():
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            else:
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]-site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherwise arg'] = state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node_to_num[branch[0]],0])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherarg prob'] = dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node_to_num[branch[0]],0]



          
##EDNECP
plist = [['EDN', 'ECP']]
tree = EDNECP_newicktree
outgroup = 'Saguinus_oedipus'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
for pair in plist:
    site_information_dict['_'.join(pair)] = {}
    site_multianalysis_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        site_information_dict['_'.join(pair)][site] = {}
        site_multianalysis_dict['_'.join(pair)][site] = {}        
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_information_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = {}
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'] = {}
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'][state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i])]=dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]
            
            for i in range(len(dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'][state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node,i])]=dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,i]
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = []
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]>0.05:
                    site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'].append((state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i]),dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]))
    
    #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，且疑似IGC，组成 multi_site_dict 
    multi_site_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        flag1 = 0
        for branch in edge_list:
            i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
            i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
            for i in range(i0):
                for j in range(i1):
                    for k in range(2):
                        count = 0
                        for l in range(3):
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                count = count + 1
                        if count >= 2:
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    flag1 = 1
            if flag1 == 1:
                break   
        if flag1 == 1:
            multi_site_dict['_'.join(pair)][site]={}
            for branch in edge_list:
                flag = 0
                i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
                i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
                for i in range(i0):
                    for j in range(i1):
                        for k in range(2):
                            count = 0
                            for l in range(3):
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                    count = count + 1
                            if count >= 2:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        flag = 1
                if flag == 1:
                    multi_site_dict['_'.join(pair)][site][branch]={}
                    for i in range(i0):
                        for j in range(i1):
                            for k in range(2):#paralog
                                count = 0
                                for l in range(3):
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                        count = count + 1
                                if count >= 2: 
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                            event = site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] + ' ' + pair[0]+' '+site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1] + ' ' + pair[1]
                                            multi_site_dict['_'.join(pair)][site][branch][event] = {}
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node2 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[1]]['IGC'][(event[5:8],event[8 + len(pair[0]) + 7:8 + len(pair[0]) + 10])]
                                            if not (event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5]) in site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'].keys():
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            else:
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]-site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherwise arg'] = state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node_to_num[branch[0]],0])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherarg prob'] = dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node_to_num[branch[0]],0]


##ERaERb
plist = [['ERa', 'ERb']]
tree = ERa_ERb_newicktree
outgroup = 'Branchiostoma_floridae'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
for pair in plist:
    site_information_dict['_'.join(pair)] = {}
    site_multianalysis_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        site_information_dict['_'.join(pair)][site] = {}
        site_multianalysis_dict['_'.join(pair)][site] = {}        
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_information_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = {}
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'] = {}
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'][state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i])]=dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]
            
            for i in range(len(dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'][state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node,i])]=dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,i]
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = []
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]>0.05:
                    site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'].append((state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i]),dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]))
    
    #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，且疑似IGC，组成 multi_site_dict 
    multi_site_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        flag1 = 0
        for branch in edge_list:
            i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
            i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
            for i in range(i0):
                for j in range(i1):
                    for k in range(2):
                        count = 0
                        for l in range(3):
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                count = count + 1
                        if count >= 2:
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    flag1 = 1
            if flag1 == 1:
                break   
        if flag1 == 1:
            multi_site_dict['_'.join(pair)][site]={}
            for branch in edge_list:
                flag = 0
                i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
                i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
                for i in range(i0):
                    for j in range(i1):
                        for k in range(2):
                            count = 0
                            for l in range(3):
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                    count = count + 1
                            if count >= 2:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        flag = 1
                if flag == 1:
                    multi_site_dict['_'.join(pair)][site][branch]={}
                    for i in range(i0):
                        for j in range(i1):
                            for k in range(2):#paralog
                                count = 0
                                for l in range(3):
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                        count = count + 1
                                if count >= 2: 
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                            event = site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] + ' ' + pair[0]+' '+site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1] + ' ' + pair[1]
                                            multi_site_dict['_'.join(pair)][site][branch][event] = {}
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node2 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[1]]['IGC'][(event[5:8],event[8 + len(pair[0]) + 7:8 + len(pair[0]) + 10])]
                                            if not (event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5]) in site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'].keys():
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            else:
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]-site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherwise arg'] = state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node_to_num[branch[0]],0])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherarg prob'] = dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node_to_num[branch[0]],0]

##ARaERa
plist = [['ARa', 'ERa']]
tree = ARa_ERa_newicktree
outgroup = 'Mus_musculus'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
for pair in plist:
    site_information_dict['_'.join(pair)] = {}
    site_multianalysis_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        site_information_dict['_'.join(pair)][site] = {}
        site_multianalysis_dict['_'.join(pair)][site] = {}        
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_information_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = {}
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'] = {}
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'][state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i])]=dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]
            
            for i in range(len(dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'][state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node,i])]=dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,i]
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = []
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]>0.05:
                    site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'].append((state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i]),dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]))
    
    #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，且疑似IGC，组成 multi_site_dict 
    multi_site_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        flag1 = 0
        for branch in edge_list:
            i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
            i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
            for i in range(i0):
                for j in range(i1):
                    for k in range(2):
                        count = 0
                        for l in range(3):
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                count = count + 1
                        if count >= 2:
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    flag1 = 1
            if flag1 == 1:
                break   
        if flag1 == 1:
            multi_site_dict['_'.join(pair)][site]={}
            for branch in edge_list:
                flag = 0
                i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
                i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
                for i in range(i0):
                    for j in range(i1):
                        for k in range(2):
                            count = 0
                            for l in range(3):
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                    count = count + 1
                            if count >= 2:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        flag = 1
                if flag == 1:
                    multi_site_dict['_'.join(pair)][site][branch]={}
                    for i in range(i0):
                        for j in range(i1):
                            for k in range(2):#paralog
                                count = 0
                                for l in range(3):
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                        count = count + 1
                                if count >= 2: 
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                            event = site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] + ' ' + pair[0]+' '+site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1] + ' ' + pair[1]
                                            multi_site_dict['_'.join(pair)][site][branch][event] = {}
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node2 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[1]]['IGC'][(event[5:8],event[8 + len(pair[0]) + 7:8 + len(pair[0]) + 10])]
                                            if not (event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5]) in site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'].keys():
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            else:
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]-site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherwise arg'] = state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node_to_num[branch[0]],0])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherarg prob'] = dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node_to_num[branch[0]],0]


##ARMRGRPR
plist =ARMRGRPR_list
tree = ARMRGRPR_newicktree
outgroup = 'Aplysia_californica'

ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
for pair in plist:
    site_information_dict['_'.join(pair)] = {}
    site_multianalysis_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        site_information_dict['_'.join(pair)][site] = {}
        site_multianalysis_dict['_'.join(pair)][site] = {}        
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_information_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = {}
            site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'] = {}
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['IGC'][state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i])]=dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]
            
            for i in range(len(dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,:])):
                site_information_dict['_'.join(pair)][site][num_to_node[node]]['tau=0'][state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node,i])]=dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node,i]
        for node in range(len(dict_all['_'.join(pair)]['IGC']['arg'][0,:,0])):
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]] = {}
            #记录似然
            site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'] = []
            for i in range(len(dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,:])):
                if dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]>0.05:
                    site_multianalysis_dict['_'.join(pair)][site][num_to_node[node]]['IGC'].append((state_to_compositecodons(dict_all['_'.join(pair)]['IGC']['arg'][site,node,i]),dict_all['_'.join(pair)]['IGC']['likelihood'][site,node,i]))
    
    #在   site_list_dict 中指示的位点寻找IGC中出现多核酸替换，且疑似IGC，组成 multi_site_dict 
    multi_site_dict['_'.join(pair)] = {}
    for site in site_list_dict['_'.join(pair)]:
        flag1 = 0
        for branch in edge_list:
            i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
            i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
            for i in range(i0):
                for j in range(i1):
                    for k in range(2):
                        count = 0
                        for l in range(3):
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                count = count + 1
                        if count >= 2:
                            if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    flag1 = 1
            if flag1 == 1:
                break   
        if flag1 == 1:
            multi_site_dict['_'.join(pair)][site]={}
            for branch in edge_list:
                flag = 0
                i0 = len(site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'])
                i1 = len(site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'])#几个备选
                for i in range(i0):
                    for j in range(i1):
                        for k in range(2):
                            count = 0
                            for l in range(3):
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                    count = count + 1
                            if count >= 2:
                                if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        flag = 1
                if flag == 1:
                    multi_site_dict['_'.join(pair)][site][branch]={}
                    for i in range(i0):
                        for j in range(i1):
                            for k in range(2):#paralog
                                count = 0
                                for l in range(3):
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][k][l]!=site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][k][l]:
                                        count = count + 1
                                if count >= 2: 
                                    if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] != site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] and site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                        if site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0]==site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] or site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] == site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1]:
                                            event = site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][0] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][0] + ' ' + pair[0]+' '+site_multianalysis_dict['_'.join(pair)][site][branch[0]]['IGC'][i][0][1] + '->' + site_multianalysis_dict['_'.join(pair)][site][branch[1]]['IGC'][j][0][1] + ' ' + pair[1]
                                            multi_site_dict['_'.join(pair)][site][branch][event] = {}
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node2 prob in IGC'] = site_information_dict['_'.join(pair)][site][branch[1]]['IGC'][(event[5:8],event[8 + len(pair[0]) + 7:8 + len(pair[0]) + 10])]
                                            if not (event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5]) in site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'].keys():
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            else:
                                                multi_site_dict['_'.join(pair)][site][branch][event]['node1 prob dif'] = abs(site_information_dict['_'.join(pair)][site][branch[0]]['IGC'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])]-site_information_dict['_'.join(pair)][site][branch[0]]['tau=0'][(event[0:3],event[8 + len(pair[0]) + 2:8 + len(pair[0]) + 5])])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherwise arg'] = state_to_compositecodons(dict_all['_'.join(pair)]['tau=0']['arg'][site,node_to_num[branch[0]],0])
                                            multi_site_dict['_'.join(pair)][site][branch][event]['node1 otherarg prob'] = dict_all['_'.join(pair)]['tau=0']['likelihood'][site,node_to_num[branch[0]],0]

#dict_all['_'.join(pair)][model]['arg']
#EXCEL 表: dataset, site, branch, event, node1 prob in IGC, node1 prob dif, node 2 prob in IGC, node1 otherwise arg, node1 otherarg prob 
#EDN_ECP    10  'N3','N7'   'TCC->TGT EDN TGT->TGT ECP'
f=open('IGCevent_result.txt','w+')
f.write('dataset\tsite\tbranch\tevent\tnode1 prob in IGC\tnode1 prob dif\tnode2 prob in IGC\tnode1 tau=0 arg\tnode1 otherarg prob\tnode1 PAML arg\n')
for paralog in multi_site_dict.keys():
    for site in multi_site_dict[paralog].keys():
        for branch in multi_site_dict[paralog][site].keys():
            f.write(str(paralog) + '\t' + str(site) + '\t' + str(branch) + '\t' + str(multi_site_dict[paralog][site][branch].keys()[0]) + '\t')
            f.write(str(multi_site_dict[paralog][site][branch][multi_site_dict[paralog][site][branch].keys()[0]]['node1 prob in IGC']) + '\t')
            f.write(str(multi_site_dict[paralog][site][branch][multi_site_dict[paralog][site][branch].keys()[0]]['node1 prob dif']) + '\t')
            f.write(str(multi_site_dict[paralog][site][branch][multi_site_dict[paralog][site][branch].keys()[0]]['node2 prob in IGC']) + '\t')
            f.write(str(multi_site_dict[paralog][site][branch][multi_site_dict[paralog][site][branch].keys()[0]]['node1 otherwise arg']) + '\t')
            f.write(str(multi_site_dict[paralog][site][branch][multi_site_dict[paralog][site][branch].keys()[0]]['node1 otherarg prob'])+ '\t')
            if branch[0] == 'N0':
                f.write(name_to_seq[paralog][branch[0]+paralog[0:(len(paralog)-1)/2]][3*site:3*site+3]+'\t')
            else:
                f.write(name_to_seq[paralog][branch[0]+paralog[0:(len(paralog)-1)/2]][3*site:3*site+3]+', '+name_to_seq[paralog][branch[0]+paralog[(len(paralog)+1)/2:len(paralog)]][3*site:3*site+3]+'\t')
            f.write('\n')
f.close()
    
#in IGC_event result.txt