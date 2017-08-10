# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 08:23:33 2017

@author: xjw1001001
"""
#only when PAML in desktop is available,the yeast version only

paralog_list = [['YLR406C', 'YDL075W'],
 ['YER131W', 'YGL189C'],
 ['YML026C', 'YDR450W'],
 ['YNL301C', 'YOL120C'],
 ['YNL069C', 'YIL133C'],
 ['YMR143W', 'YDL083C'],
 ['YJL177W', 'YKL180W'],
 ['YBR191W', 'YPL079W'],
 ['YER074W', 'YIL069C'],
 ['YDR418W', 'YEL054C'],
 ['YBL087C', 'YER117W'],
 ['YLR333C', 'YGR027C'],
 ['YMR142C', 'YDL082W'],
 ['YER102W', 'YBL072C'],
 ]

for pair in paralog_list:
    primalline=[]
    fastaline=[]
    with open('/Users/xjw1001001/Desktop/PAML/output/' + '_'.join(pair) +'/out/construct.fasta','r') as f:
        for line in f.readlines():
            primalline.append(line)
            if 'node #14' in line:
                continue
            sline = '>' + line
            sline=sline.replace(' ','')
            sline=sline.replace('\n','')
            sline=sline.replace('node#15','N0'+pair[0])
            for i in range(5):
                sline=sline.replace('node#' + str(15+1+i),'N'+str(1+i)+pair[1])
                sline=sline.replace('node#' + str(20+1+i),'N'+str(1+i)+pair[0])
            sline=sline.replace(pair[0],pair[0] + '\n')
            sline=sline.replace(pair[1],pair[1] + '\n')
            fastaline.append(sline)      
    f1 = open('/Users/xjw1001001/Desktop/PAML/PAMLfasta/PAML_' + '_'.join(pair) +'.fasta','w+')
    for line in fastaline:
        f1.write(line)
        f1.write('\n')
    f1.close()

