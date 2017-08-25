# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 08:23:33 2017

@author: xjw1001001
"""
#only when PAML in desktop is available,the yeast version only
from Bio import Seq, SeqIO, AlignIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np
paralog_list = [['YLR406C', 'YDL075W'],
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

for pair in paralog_list:
    primalline=[]
    fastaline=[]
    with open('/Users/xjw1001001/Desktop/PAML/output/' + '_'.join(pair) +'/out/construct.fasta','r') as f:
        for line in f.readlines():
            primalline.append(line)
            sline = '>' + line
            sline=sline.replace('node #14','Root'+pair[0])
            sline=sline.replace(' ','')
            sline=sline.replace('\n','')
            sline=sline.replace('node#15','N0'+pair[0])
            for i in range(5):
                sline=sline.replace('node#' + str(15+1+i),'N'+str(1+i)+pair[1])
                sline=sline.replace('node#' + str(20+1+i),'N'+str(1+i)+pair[0])
            sline=sline.replace(pair[0],pair[0] + '\n')
            sline=sline.replace(pair[1],pair[1] + '\n')
            fastaline.append(sline)      
    f1 = open('/Users/xjw1001001/Desktop/PAML/PAMLfasta/PAML_' + '_'.join(pair) +'.fasta','w+')
    for line in fastaline:
        f1.write(line)
        f1.write('\n')
    f1.close()

#ERa_ERb
pair = ['ERa','ERb']
primalline=[]
fastaline=[]
substitution_dict = {'node#39':'N14ERa','node#38':'N8ERa','node#37':'N7ERa','node#36':'N6ERa','node#41':'N9ERa','node#40':'N5ERa'
                     ,'node#35':'N4ERa','node#44':'N13ERa','node#46':'N12ERa','node#47':'N11ERa','node#45':'N10ERa'
                     ,'node#43':'N3ERa','node#42':'N2ERa','node#34':'N1ERa'
                     ,'node#53':'N14ERb','node#52':'N8ERb','node#51':'N7ERb','node#50':'N6ERb','node#55':'N9ERb','node#54':'N5ERb'
                     ,'node#49':'N4ERb','node#58':'N13ERb','node#60':'N12ERb','node#61':'N11ERb','node#59':'N10ERb'
                     ,'node#57':'N3ERb','node#56':'N2ERb','node#48':'N1ERb'}
with open('/Users/xjw1001001/Desktop/PAML/output/' + '_'.join(pair) +'/out/construct.fasta','r') as f:
    for line in f.readlines():
        primalline.append(line)
        sline = '>' + line
        sline=sline.replace('node #32','Root'+pair[0])
        sline=sline.replace(' ','')
        sline=sline.replace('\n','')
        sline=sline.replace('node#33','N0'+pair[0])
        for i in substitution_dict.keys():
            sline=sline.replace(i,substitution_dict[i])
        sline=sline.replace(pair[0],pair[0] + '\n')
        sline=sline.replace(pair[1],pair[1] + '\n')
        fastaline.append(sline)      
f1 = open('/Users/xjw1001001/Desktop/PAML/PAMLfasta/PAML_' + '_'.join(pair) +'.fasta','w+')
for line in fastaline:
    f1.write(line)
    f1.write('\n')
f1.close()

#ARa_ERa
pair = ['ARa','ERa']
primalline=[]
fastaline=[]
substitution_dict = {'node#36':'N12ERa','node#35':'N11ERa','node#34':'N7ERa','node#33':'N6ERa','node#32':'N5ERa','node#37':'N8ERa'
                     ,'node#31':'N4ERa','node#41':'N10ERa','node#40':'N9ERa','node#39':'N3ERa','node#38':'N2ERa'
                     ,'node#30':'N1ERa'
                     ,'node#48':'N12ARa','node#47':'N11ARa','node#46':'N7ARa','node#45':'N6ARa','node#44':'N5ARa','node#49':'N8ARa'
                     ,'node#43':'N4ARa','node#53':'N10ARa','node#52':'N9ARa','node#51':'N3ARa','node#50':'N2ARa'
                     ,'node#42':'N1ARa','node#29':'N0ERa','node#28':'RootERa'}
with open('/Users/xjw1001001/Desktop/PAML/output/' + '_'.join(pair) +'/out/construct.fasta','r') as f:
    for line in f.readlines():
        primalline.append(line)
        sline = '>' + line
        sline=sline.replace(' ','')
        sline=sline.replace('\n','')
        for i in substitution_dict.keys():
            sline=sline.replace(i,substitution_dict[i])
        sline=sline.replace(pair[0],pair[0] + '\n')
        sline=sline.replace(pair[1],pair[1] + '\n')
        fastaline.append(sline)      
f1 = open('/Users/xjw1001001/Desktop/PAML/PAMLfasta/PAML_' + '_'.join(pair) +'.fasta','w+')
for line in fastaline:
    f1.write(line)
    f1.write('\n')
f1.close()

#ARGRMRPR
pairlist = [['AR', 'MR'],
                     ['AR', 'GR'],
                     ['AR', 'PR'],
                     ['MR', 'GR'],
                     ['MR', 'PR'],
                     ['PR', 'GR']]
for pair in pairlist:
    primalline=[]
    fastaline=[]
    substitution_dict = {'node#25':'N4'+pair[0],'node#31':'N9'+pair[0],'node#30':'N7'+pair[0]
                         ,'node#32':'N8'+pair[0],'node#29':'N6'+pair[0],'node#28':'N5'+pair[0]
                         ,'node#27':'N3'+pair[0],'node#26':'N2'+pair[0],'node#24':'N1'+pair[0]
                         ,'node#34':'N4'+pair[1],'node#40':'N9'+pair[1],'node#39':'N7'+pair[1]
                         ,'node#41':'N8'+pair[1],'node#38':'N6'+pair[1],'node#37':'N5'+pair[1]
                         ,'node#36':'N3'+pair[1],'node#35':'N2'+pair[1],'node#33':'N1'+pair[1]
                         ,'node#23':'N0'+pair[0],'node#22':'ROOT'+pair[0]
                         }
    with open('/Users/xjw1001001/Desktop/PAML/output/' + '_'.join(pair) +'/out/construct.fasta','r') as f:
        for line in f.readlines():
            primalline.append(line)
            sline = '>' + line
            sline=sline.replace(' ','')
            sline=sline.replace('\n','')
            for i in substitution_dict.keys():
                sline=sline.replace(i,substitution_dict[i])
            sline=sline.replace(pair[0],pair[0] + '\n')
            sline=sline.replace(pair[1],pair[1] + '\n')
            fastaline.append(sline)      
    f1 = open('/Users/xjw1001001/Desktop/PAML/PAMLfasta/PAML_' + '_'.join(pair) +'.fasta','w+')
    for line in fastaline:
        f1.write(line)
        f1.write('\n')
    f1.close()

PAML_parameter_dict = {}
path = '/Users/xjw1001001/Desktop/PAML/'
paralog_list = [['YLR406C', 'YDL075W'],#pair#TODO: other data
 ['YER131W', 'YGL189C'], ['YML026C', 'YDR450W'], ['YNL301C', 'YOL120C'], ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'], ['YJL177W', 'YKL180W'], ['YBR191W', 'YPL079W'], ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'], ['YBL087C', 'YER117W'], ['YLR333C', 'YGR027C'], ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'], ['EDN', 'ECP'],['ERa', 'ERb'],['AR', 'MR'],['AR', 'GR'],['AR', 'PR'],
 ['MR', 'GR'],['MR', 'PR'],['PR', 'GR'] ] 
for pair in paralog_list:#parameters: kappa(-5), omega(-1), tau,branches 
    PAML_parameter_dict['_'.join(pair)] = {}
    codeml_result = codeml.read(path+'output/' + '_'.join(pair) + '/out/' + '_'.join(pair) + '_codeml')
    #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
    parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
    PAML_parameter_dict['_'.join(pair)]['kappa'] =  parameter_list[-5]
    PAML_parameter_dict['_'.join(pair)]['omega'] =  parameter_list[-1]
    
        