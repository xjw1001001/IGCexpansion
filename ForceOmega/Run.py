## This script is to run model under omega being fixed as 1.0
## This is to test how much omega value would affect the inferrence
## 05/10/2016
## Xiang Ji

from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './YeastTree.newick'
    if args.force:
        if args.model == 'MG94':
            Force = {4:1.0}
        elif args.model == 'HKY':
            print ('HKY model does not have omega parameter!')            
            Force = {3:1.0}
    else:
        Force = None
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = args.model, Force = Force, clock = args.clock)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './Summary/')
    test.get_SitewisePosteriorSummary(summary_path = './Summary/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required = True, help = 'Substitution Model')
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
    
    main(parser.parse_args())
