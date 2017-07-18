# -*- coding: utf-8 -*-
import os

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    sh_line = 'sbatch -o ./cluster_outs/IGC-%j.out --mail-type=Fail --mail-user=1176434052@qq.com ./ShFiles/'

    model = 'rebuild'
    IGC_bash_file = './' + model + '_IGC.sh'
    
    with open(IGC_bash_file, 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for paralog in pairs:
            f.write(sh_line + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh \n')
            with open('./ShFiles/' + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run.py ' + ' --paralog1 ' + paralog[0]
                        + ' --paralog2 ' + paralog[1] + ' --no-clock' + '\n')


##    sh_line = 'sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles_multiply3/'
##
##    for IGC_geo in IGC_geo_list:
##        IGC_geo_sh_file = './' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_multiply3.sh'
##        with open(IGC_geo_sh_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for num_sim in sim_num_list:
##                f.write(sh_line + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh \n')
##                with open('./ShFiles_multiply3/' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python RunSimulation_multiply3.py --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')
