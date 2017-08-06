# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 03:13:20 2017

@author: xjw1001001
"""

# -*- coding: utf-8 -*-
from IGCexpansion.CodonGeneconv import ReCodonGeneconv

    
if __name__ == '__main__':
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primateoutcome processed.fasta'
    newicktree = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
    
    #MG94+tau
    MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './test/save/')
    #MG94_tau.get_mle(True, True, 0, 'BFGS')
    MG94_tau.site_reconstruction()
    MG94_tau_series = MG94_tau.reconstruction_series
    MG94_tau.get_individual_summary(summary_path = './test/Summary/')
    MG94_tau.get_SitewisePosteriorSummary(summary_path = './test/Summary/')
    
    
    '''
    #MG94
    MG94 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = None, save_path = './test/save/')
    #MG94.get_mle(True, True, 0, 'BFGS')
    MG94.site_reconstruction()#存了fasta,Pamlrebuild要自己处理
    MG94_series = MG94.reconstruction_series
    result = MG94_tau.find_differences_between(MG94_tau_series, MG94_series)
    ''''''
    HKY_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = None, save_path = './test/save/')
    HKY_tau.get_mle(True, True, 0, 'BFGS')
    HKY_tau.site_reconstruction()
    HKY_tau_series = HKY_tau.reconstruction_series
    
    #MG94
    HKY = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = {4:0.0}, clock = None, save_path = './test/save/')
    HKY.get_mle(True, True, 0, 'BFGS')
    HKY.site_reconstruction()
    HKY_series = HKY.reconstruction_series
    result = HKY_tau.find_differences_between(HKY_tau_series, HKY_series)
    '''
    