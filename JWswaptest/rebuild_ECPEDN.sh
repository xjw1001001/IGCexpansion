#!/bin/bash
sbatch -c 5 -o ../cluster_outs/IGC-%j.out --mail-type=ALL --mail-user=1176434052@qq.com 1.sh 