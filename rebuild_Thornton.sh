#!/bin/bash
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonARPR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonARMR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonARGR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonMRGR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonMRPR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonPRGR.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonERaERb.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/ThorntonARaERa.sh 