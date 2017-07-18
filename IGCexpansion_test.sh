#!/bin/bash
#SBATCH -J IGCtest1
sbatch -o IGCtest-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/test_MG94_HKY_YL.sh 
