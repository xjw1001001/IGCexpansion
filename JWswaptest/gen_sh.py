# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 01:12:02 2017

@author: xjw1001001
"""

import os

if __name__ == '__main__':
    sh_line = 'sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/'

    model = 'swapEDNECP'
    IGC_bash_file = './' + model + '_IGC.sh'
    
    with open(IGC_bash_file, 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for tau in [0.2*i for i in range(11)]:
            f.write(sh_line + 'tau=' + str(tau) + '.sh \n')
            with open('./ShFiles/' + 'tau=' + str(tau) + '.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python swaptest_EDN_ECP.py ' + ' --tau ' + str(tau) + '\n')