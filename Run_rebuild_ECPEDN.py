# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 03:13:20 2017

@author: xjw1001001
"""

# -*- coding: utf-8 -*-
from IGCexpansion.CodonGeneconv import ReCodonGeneconv

def main():
    paralog = ['EDN', 'ECP']
    Force = None
    alignment_file = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primateoutcome processed.fasta'
    newicktree = './reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
    
    #MG94+tau
    MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './test/save/')
    MG94_tau.get_mle(True, True, 0, 'BFGS')
    MG94_tau.site_reconstruction()
    MG94_tau_series = MG94_tau.reconstruction_series
    #MG94_tau.get_individual_summary(summary_path = './test/Summary/')
    #MG94_tau.get_SitewisePosteriorSummary(summary_path = './test/Summary/')

    #MG94
    MG94 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = None, save_path = './test/save/')
    MG94.get_mle(True, True, 0, 'BFGS')
    MG94.site_reconstruction()
    MG94_series = MG94.reconstruction_series
    result = MG94_tau.find_differences_between(MG94_tau_series, MG94_series)
    
if __name__ == '__main__':
    main()
    