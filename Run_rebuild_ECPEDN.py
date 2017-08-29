# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 03:13:20 2017

@author: xjw1001001
"""
#EDNECP树中有两个枝长很小
# -*- coding: utf-8 -*-
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
# AIC=2k-2ln(L)
    
if __name__ == '__main__':
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primateoutcome processed.fasta'
    newicktree = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
    
    #MG94+tau
    MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './test/save/')
    
    #MG94_tau.get_mle(True, True, 0, 'BFGS')
    MG94_tau.site_reconstruction()
    
    #MG94_tau.Expected_tau_for_sitewise_and_branchwise()
    MG94_tau_series = MG94_tau.reconstruction_series
    MG94_tau_likelihooddict = MG94_tau.likelihood_dict
    
    #MG94_tau.get_individual_summary(summary_path = './test/Summary/')
    #MG94_tau.get_SitewisePosteriorSummary(summary_path = './test/Summary/')
    
    
    
    #MG94
    MG94 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = None, save_path = './test/save/')
    #MG94.get_mle(True, True, 0, 'BFGS')
    MG94.site_reconstruction()#存了fasta,Pamlrebuild要自己处理
    #MG94.Expected_tau_for_sitewise_and_branchwise()
    MG94_series = MG94.reconstruction_series
    MG94_likelihooddict = MG94.likelihood_dict
    result = MG94_tau.find_differences_between(MG94_tau_series, MG94_series)
    