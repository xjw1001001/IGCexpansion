#!/bin/bash
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=0.0.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=0.2.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=0.4.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=0.6.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=0.8.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=1.0.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=1.2.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=1.4.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=1.6.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=1.8.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/tau=2.0.sh 
