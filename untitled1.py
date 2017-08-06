# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 03:02:47 2017

@author: xjw1001001
"""
A={}
A[18]=states_matrix_IGC[18,:,:]
A[26]=states_matrix_IGC[26,:,:]
A[29]=states_matrix_IGC[29,:,:]
A[47]=states_matrix_IGC[47,:,:]
A[55]=states_matrix_IGC[55,:,:]
A[61]=states_matrix_IGC[61,:,:]
A[65]=states_matrix_IGC[65,:,:]
A[76]=states_matrix_IGC[76,:,:]
A[84]=states_matrix_IGC[84,:,:]
A[85]=states_matrix_IGC[85,:,:]
A[87]=states_matrix_IGC[87,:,:]
A[90]=states_matrix_IGC[90,:,:]
A[92]=states_matrix_IGC[92,:,:]
A[93]=states_matrix_IGC[93,:,:]
A[95]=states_matrix_IGC[95,:,:]
A[99]=states_matrix_IGC[99,:,:]
A[102]=states_matrix_IGC[102,:,:]
A[115]=states_matrix_IGC[115,:,:]
A[120]=states_matrix_IGC[120,:,:]
A[125]=states_matrix_IGC[125,:,:]
A[128]=states_matrix_IGC[128,:,:]

A[122]=states_matrix_IGC[122,:,:]
A[155]=states_matrix_IGC[155,:,:]

B={}
B[18]=states_matrix[18,:,:]
B[26]=states_matrix[26,:,:]
B[29]=states_matrix[29,:,:]
B[47]=states_matrix[47,:,:]
B[55]=states_matrix[55,:,:]
B[61]=states_matrix[61,:,:]
B[65]=states_matrix[65,:,:]
B[76]=states_matrix[76,:,:]
B[84]=states_matrix[84,:,:]
B[85]=states_matrix[85,:,:]
B[87]=states_matrix[87,:,:]
B[90]=states_matrix[90,:,:]
B[92]=states_matrix[92,:,:]
B[93]=states_matrix[93,:,:]
B[95]=states_matrix[95,:,:]
B[99]=states_matrix[99,:,:]
B[102]=states_matrix[102,:,:]
B[115]=states_matrix[115,:,:]
B[120]=states_matrix[120,:,:]
B[125]=states_matrix[125,:,:]
B[128]=states_matrix[128,:,:]

B[122]=states_matrix[122,:,:]
B[155]=states_matrix[155,:,:]







IGClikelihood={}
sitelist = [18,26,29,47,55,61,65,76,84,85,87,90,92,93,95,99,102,115,120,125,128,122,155]
for site in sitelist:  
    IGClikelihood[site]={}  
    for node in range(len(self.node_to_num)):   
        IGClikelihood[site][node]={}      
        for i in bn.argpartition(-A[site][:,node], 7)[:7]:  
            IGClikelihood[site][node][i] = A[site][i,node]
f = open('./test/Ancestral_reconstruction/matrix/IGClikelihooddict.txt', 'wb')
pickle.dump(IGClikelihood, f)
f.close()

NoIGClikelihood={}
sitelist = [18,26,29,47,55,61,65,76,84,85,87,90,92,93,95,99,102,115,120,125,128,122,155]
for site in sitelist:  
    NoIGClikelihood[site]={}  
    for node in range(len(self.node_to_num)):     
        NoIGClikelihood[site][node]={}      
        for i in bn.argpartition(-B[site][:,node], 7)[:7]:      
            NoIGClikelihood[site][node][i] = B[site][i,node]
            
f = open('./test/Ancestral_reconstruction/matrix/NoIGClikelihooddict.txt', 'wb')
pickle.dump(NoIGClikelihood, f)
f.close()

sites=930
state_1, state_2 = divmod(sites, 61)
state_1 = self.state_to_codon[int(state_1)]
state_2 = self.state_to_codon[int(state_2)]
[state_1,state_2]

import numpy as np
AIGCmax=np.load('./test/Ancestral_reconstruction/matrix/EDNECPwithIGC.npy')
AnoIGCmax=np.load('./test/Ancestral_reconstruction/matrix/EDNECPwithoutIGC.npy')

import cPickle as pickle
f1 = file('./test/Ancestral_reconstruction/matrix/IGClikelihooddict.txt', 'rb') 
f2 = file('./test/Ancestral_reconstruction/matrix/NoIGClikelihooddict.txt', 'rb') 
AIGClikelihood = pickle.load(f1) 
AnoIGClikelihood = pickle.load(f2) 

from Bio import SeqIO
seq_dictMAPL = SeqIO.to_dict(SeqIO.parse( './reconstruction_data/mafft-win/MAPLreconst EDNECP.fasta', "fasta" ))
seq_dictIGC = SeqIO.to_dict(SeqIO.parse( './reconstruction_data/mafft-win/ancestral_reconstruction_EDN_ECP_MG94+IGC.fasta', "fasta" ))
seq_dictIGC_tau_0 = SeqIO.to_dict(SeqIO.parse( './reconstruction_data/mafft-win/ancestral_reconstruction_EDN_ECP_MG94.fasta', "fasta" ))

seq_dictIGC['N1ECP'].seq#like string

flag=0#paml reubild
result = {}
for nodes in seq_dictIGC_tau_0.keys():
    result[nodes] = {}
    for site in range(len(seq_dictMAPL[nodes].seq)):
        if seq_dictMAPL[nodes].seq[site] == seq_dictIGC_tau_0[nodes].seq[site]:
            continue
        else:
            
            result[nodes][site] = {}
            result[nodes][site]['IGC'] = seq_dictMAPL[nodes].seq[site]
            result[nodes][site]['PAML'] = seq_dictIGC_tau_0[nodes].seq[site]
            flag=flag+1
print (flag)    

#find the max
tt=7
likelihood_temp=[]
argmatrix=np.zeros((156,25,tt))  
likelihood_matrix=np.zeros((156,25,tt))  

for site in range(self.nsites):
    likelihood_temp.append([])
    for node in range(len(self.node_to_num)):
        likelihood_temp[site].append([])
        likelihood_temp[site][node]={}
        for i in range(tt):      
            likelihood_temp[site][node][np.argpartition(-states_matrix[site,0:3721,node], tt)[i]]=states_matrix[site,0:3721,node][np.argpartition(-states_matrix[site,0:3721,node], tt)[i]]
#sort
likelihood_dict=[]
for site in range(self.nsites):
    likelihood_dict.append([])
    for node in range(len(self.node_to_num)):
        likelihood_dict[site].append([])
        likelihood_dict[site][node]=sorted(likelihood_temp[site][node].items(), lambda x, y: cmp(x[1], y[1]), reverse=True)
        for i in range(tt):  
            (argmatrix[site,node,i],likelihood_matrix[site,node,i])=likelihood_dict[site][node][i]
self.likelihood_dict=likelihood_dict
#save as numpy
if self.tau == 0:
    model = self.Model + '_tau=0'
else:
    model = self.Model + '_IGC'

for node in range(len(self.node_to_num)):    
    np.savetxt(open('./test/Ancestral_reconstruction/matrix/' + 'ancestral_reconstruction_' + self.paralog[0] + '_' + self.paralog[1] + '_' +model + '_node_' + str(node) +'.txt', 'w+'), argmatrix[:,node,:])
    np.savetxt(open('./test/Ancestral_reconstruction/matrix/' + 'ancestral_reconstruction_' + self.paralog[0] + '_' + self.paralog[1] + '_' +model + '_node_' + str(node) +'.txt', 'w+'), likelihood_matrix[:,node,:])
        