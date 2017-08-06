# -*- coding: utf-8 -*-
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    Force = None
    alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './YeastTree.newick'
    
    #MG94+tau
    MG94_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './test/save/')
    MG94_tau.get_mle(True, True, 0, 'BFGS')
    MG94_tau.site_reconstruction()
    MG94_tau_series = MG94_tau.reconstruction_series
    MG94_tau_likelihooddict = MG94_tau.likelihood_dict
    MG94_tau.get_individual_summary(summary_path = './test/Summary/')
    MG94_tau.get_SitewisePosteriorSummary(summary_path = './test/Summary/')

    #MG94
    MG94 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = None, save_path = './test/save/')
    MG94.get_mle(True, True, 0, 'BFGS')
    MG94.site_reconstruction()
    MG94_series = MG94.reconstruction_series
    MG94_likelihooddict = MG94.likelihood_dict
    result = MG94_tau.find_differences_between(MG94_tau_series, MG94_series)
    
    '''
    #HKY+tau
    HKY_tau = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = None, save_path = './test/save/')
    HKY_tau.get_mle(True, True, 0, 'BFGS')
    HKY_tau.site_reconstruction()
    HKY_tau_series = HKY_tau.reconstruction_series

    #HKY
    HKY = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = {4:0.0}, clock = None, save_path = './test/save/')
    HKY.get_mle(True, True, 0, 'BFGS')
    HKY.site_reconstruction()
    HKY_series = HKY.reconstruction_series
    result = HKY_tau.find_differences_between(HKY_tau_series, HKY_series)
    '''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
    
    main(parser.parse_args())

#python Runrebuild.py --paralog1 YBL087C --paralog2 YER117W --no-clock
#chmod +x rebuild_IGC.sh
#./rebuild_IGC.sh