# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 21:25:36 2017

@author: xjw1001001
"""

import numpy as np

path = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/'
reconstruct_path = 'matrix/reconstruction_likelihood/ancestral_reconstruction_'
paralog_list = [['YLR406C', 'YDL075W'],#pair
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
EDN_pair = ['EDN', 'ECP']
Yeast_node = 13
EDNECP_node = 25
Model_list = ['IGC','tau=0']#model
#'_'.join(pair) + '_MG94_' + model + '_node_' + node +'.txt'

dict_all = {}
for pair in paralog_list:
    # yeast read data
    dict_all[pair]={}
    for model in Model_list:
        dict_all[pair][model]=[]
        
    
    
    
    a=np.loadtxt(path + reconstruct_path + 'EDN_ECP_MG94_IGC_node_0.txt')