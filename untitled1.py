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
B[95]=states_matrix[95,:,:]
B[99]=states_matrix[99,:,:]
B[102]=states_matrix[102,:,:]
B[115]=states_matrix[115,:,:]
B[120]=states_matrix[120,:,:]
B[125]=states_matrix[125,:,:]
B[128]=states_matrix[128,:,:]

B[122]=states_matrix[122,:,:]
B[155]=states_matrix[155,:,:]

for site in A.keys():
     np.save('./test/Ancestral_reconstruction/matrix/likelihood_IGC_' + str(site+1) + 'th.npy', A[site])
for site in B.keys():
     np.save('./test/Ancestral_reconstruction/matrix/likelihood_no_IGC_' + str(site+1) + 'th.npy', B[site])

