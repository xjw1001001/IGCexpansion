#!/bin/bash
sbatch -o ./cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com ./ShFiles/1.sh
