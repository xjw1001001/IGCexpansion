# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 21:47:06 2017

@author: xjw1001001
"""
import numpy as np
from IGCexpansion.CodonGeneconFunc import *
import xlwt
import datetime
ezxf = xlwt.easyxf
def write_xls(file_name, sheet_name_list, headings, data_dict, heading_xf, data_xfs):
    book = xlwt.Workbook()
    for sheet_name in sheet_name_list:
        sheet = book.add_sheet(sheet_name)
        rowx = 0
        for colx, value in enumerate(headings):
            sheet.write(rowx, colx, value, heading_xf)
        sheet.set_panes_frozen(True) # frozen headings instead of split panes
        sheet.set_horz_split_pos(rowx+1) # in general, freeze after last heading row
        sheet.set_remove_splits(True) # if user does unfreeze, don't leave a split there
        for row in data[sheet_name]:
            rowx += 1
            for colx, value in enumerate(row):
                sheet.write(rowx, colx, value, data_xfs[colx])
    book.save(file_name)


path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/matrix/sitewise_IGC_statmatrix/'
paralog_list = [['YLR406C', 'YDL075W'],#pair#TODO: other data
 ['YER131W', 'YGL189C'], ['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'], ['YLR333C', 'YGR027C'], ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'], ['EDN', 'ECP'] ] 
EDNECP_newicktree ='/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
Yeast_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/YeastTree.newick'

#read data
Expected_tau = {} #paralog array
ExpectedIGC = {}
ExpectedIGC['num'] = {}
ExpectedIGC['1to2'] = {}
ExpectedIGC['2to1'] = {}

posterior = {}
posterior['1to2'] = {}
posterior['2to1'] = {}
#'_'.join(pair)
for pair in paralog_list:
    Expected_tau['_'.join(pair)] = np.loadtxt(open(path + 'Expected_tau/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['num']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['1to2']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum1_2/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    ExpectedIGC['2to1']['_'.join(pair)] = np.loadtxt(open(path + 'ExpectedIGCnum2_1/' + '_'.join(pair) + '_MG94_IGC.txt','r'))
    posterior['1to2']['_'.join(pair)] = np.loadtxt(open(path + 'posterior/' + '_'.join(pair) + '_MG94_IGC_no1to2.txt','r'))
    posterior['2to1']['_'.join(pair)] = np.loadtxt(open(path + 'posterior/' + '_'.join(pair) + '_MG94_IGC_no2to1.txt','r'))
    #posterior['1to2']['_'.join(pair)] = np.ones([len(posterior['1to2']['_'.join(pair)][:,0]),len(posterior['1to2']['_'.join(pair)][0,:])]) - np.exp(posterior['1to2']['_'.join(pair)])
    #posterior['2to1']['_'.join(pair)] = np.ones([len(posterior['2to1']['_'.join(pair)][:,0]),len(posterior['2to1']['_'.join(pair)][0,:])]) - np.exp(posterior['2to1']['_'.join(pair)])
    
#write in excels  m = ExpectedIGC['num']['_'.join(pair)].tolist()
#EDNECP
plist = [['EDN', 'ECP']]
tree = EDNECP_newicktree
path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/matrix/excel_results/'
ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
control_list = ['IGCnum','1to2num','2to1num','posterior1to2','posterior2to1']
data = {}
for pair in plist:
    for control in control_list:    
        if control == 'IGCnum':
            data[control] = ExpectedIGC['num']['_'.join(pair)].tolist()
        elif control == '1to2num':
            data[control] = ExpectedIGC['1to2']['_'.join(pair)].tolist()
        elif control == '2to1num':
            data[control] = ExpectedIGC['2to1']['_'.join(pair)].tolist()
        elif control == 'posterior1to2':
            data[control] = posterior['1to2']['_'.join(pair)].tolist()
        elif control == 'posterior2to1':
            data[control] = posterior['2to1']['_'.join(pair)].tolist()
    hdngs = [branch[0]+ ', ' + branch[1] for branch in edge_list]
    heading_xf = ezxf('font: bold on; align: wrap on, vert centre, horiz center')
    data_xfs = [ezxf(num_format_str='#0.000000') for k in range(len(hdngs))]
        
    write_xls(path + '_'.join(pair) + '_IGCanalysis.xls', control_list, hdngs, data, heading_xf, data_xfs)
    
#yeast
plist = [['YLR406C', 'YDL075W'],
 ['YER131W', 'YGL189C'], ['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'], ['YLR333C', 'YGR027C'], ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'] ] 
tree = Yeast_newicktree
path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/matrix/excel_results/'
ktree, edge_list, node_to_num = read_newick(tree, 'N1')
num_to_node = {node_to_num[i]:i for i in node_to_num}
control_list = ['IGCnum','1to2num','2to1num','posterior1to2','posterior2to1']
data = {}
for pair in plist:
    for control in control_list:    
        if control == 'IGCnum':
            data[control] = ExpectedIGC['num']['_'.join(pair)].tolist()
        elif control == '1to2num':
            data[control] = ExpectedIGC['1to2']['_'.join(pair)].tolist()
        elif control == '2to1num':
            data[control] = ExpectedIGC['2to1']['_'.join(pair)].tolist()
        elif control == 'posterior1to2':
            data[control] = posterior['1to2']['_'.join(pair)].tolist()
        elif control == 'posterior2to1':
            data[control] = posterior['2to1']['_'.join(pair)].tolist()
    hdngs = [branch[0]+ ', ' + branch[1] for branch in edge_list]
    heading_xf = ezxf('font: bold on; align: wrap on, vert centre, horiz center')
    data_xfs = [ezxf(num_format_str='#0.000000') for k in range(len(hdngs))]
        
    write_xls(path + '_'.join(pair) + '_IGCanalysis.xls', control_list, hdngs, data, heading_xf, data_xfs)

