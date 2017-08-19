#!/bin/bash
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/Thornton1.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/Thornton2.sh 
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/Thornton3.sh 