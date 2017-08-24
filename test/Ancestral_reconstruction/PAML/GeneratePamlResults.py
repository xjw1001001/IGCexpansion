import os
import subprocess
from Bio import Seq, SeqIO, AlignIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np

def partition_seqfile(seqfile, partitioned_seqfile):
    assert(os.path.isfile(seqfile))
    align = AlignIO.read(seqfile, 'fasta')
    align_length = align.get_alignment_length()
    with open(partitioned_seqfile, 'w+') as f:
        f.write('    '.join([' ', '13', str(align_length), 'GC \n']))#13 for Yeast, 25 for EDN
        with open(seqfile, 'r') as g:
            for line in g:
                if line[0] == '>':
                    name = line[1:-1]
                else:
                    f.write(name + '    ' + line)
                
    

if __name__ == '__main__':
    #yeast datas
    path = '/Users/xjw1001001/Desktop/PAML/'
    
    pair_file = '/Users/xjw1001001/Desktop/PAML/Filtered_pairs.txt'
    pairs = []
    SRlist = [['AR', 'MR'],
                     ['AR', 'GR'],
                     ['AR', 'PR'],
                     ['MR', 'GR'],
                     ['MR', 'PR'],
                     ['PR', 'GR']]
    with open(pair_file, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    tree_pair = ['YML026C', 'YDR450W']
    with open(path + 'YeastTree.newick', 'r') as f:
        all_tree_lines = f.readlines()

    with open(path + 'codeml_tail.ctl', 'r') as f:
        all_codeml_ctl_lines = f.readlines()

    with open(path + 'baseml_tail.ctl', 'r') as f:
        all_baseml_ctl_lines = f.readlines()
    



    #pairs = [pairs[0]]
    for pair in pairs:
        codeml_dir = path + 'output/' + '_'.join(pair) + '/codeml'
        baseml_dir = path + 'output/' + '_'.join(pair) + '/baseml'

        print 'Now run paml on pair ' + ' '.join(pair)
        seqfile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta'
        partitioned_seqfile = seqfile.replace('input.fasta', 'partitioned.fasta')
        if not os.path.isdir(path + 'output/' + '_'.join(pair)):
            os.mkdir(path + 'output/' + '_'.join(pair))
            
        if not os.path.isfile(seqfile):
            mafft_aligned_file = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/MafftAlignment/' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta'
            os.system(' '.join(['cp', mafft_aligned_file, seqfile]))

        #if not os.path.isfile(partitioned_seqfile):
        partition_seqfile(seqfile, partitioned_seqfile)
        
        treefile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_tree.newick'        
        with open(treefile, 'w+') as f:
            for line in all_tree_lines:
                new_line = line.replace(tree_pair[0], ''+pair[0])
                new_line = new_line.replace(tree_pair[1], ''+pair[1])
                f.write(new_line)


        outfile_baseml = path + 'output/' + '_'.join(pair) + '/out/' + '_'.join(pair) + '_baseml'
        outfile_codeml = path + 'output/' + '_'.join(pair) + '/out/' + '_'.join(pair) + '_codeml'
        baseml_ctlfile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml_control.ctl'
        codeml_ctlfile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_control.ctl'
        with open(baseml_ctlfile, 'w+') as f:
            f.writelines(['seqfile = ' + partitioned_seqfile + '\n', 'treefile = ' + treefile + '\n', 'outfile = ' + outfile_baseml + '\n'])
            f.writelines(all_baseml_ctl_lines)

        with open(codeml_ctlfile, 'w+') as f:
            f.writelines(['seqfile = ' + partitioned_seqfile + '\n', 'treefile = ' + treefile + '\n', 'outfile = ' + outfile_codeml + '\n'])
            f.writelines(all_codeml_ctl_lines)
        
        baseml_cmd = [baseml_dir, path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml_control.ctl']
        codeml_cmd = [codeml_dir, path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_control.ctl']
        #subprocess.check_output(baseml_cmd)
        #workpath =path + 'output/' + '_'.join(pair) + '/out/'
        #os.chdir(workpath)
        #os.system(' '.join(codeml_cmd))

    summary_mat = []
    finished_list = []
    label = ['MG94_baseml_tree_length', 'MG94_baseml_lnL', 'MG94_baseml_kappa', 'MG94_r2', 'MG94_r3']
    footer = ' '.join(label)


    #pairs = pairs[0:2]
    for pair in pairs:#parameters:%AG, %A, %C, kappa(-5), omega(-1), tau,branches 
        codeml_result = codeml.read(path+'output/' + '_'.join(pair) + '/out/' + '_'.join(pair) + '_codeml')
        #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
        parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                            codeml_result['NSsites'][0]['lnL'],
                            float(parameter_list[-1]),
                            float(parameter_list[-3]),
                            float(parameter_list[-2])]
                           )
        finished_list.append(pair)
        np.savetxt(open(path+'output/summary/' + '_'.join(pair) + '.txt', 'w+'), np.array([codeml_result['NSsites'][0]['lnL'],2*(len(codeml_result['NSsites'][0]['parameters']['branches'])+5)-2*codeml_result['NSsites'][0]['lnL']]))
 
    codeml_result = codeml.read('/Users/xjw1001001/Desktop/PAML/output/EDN_ECP/out/EDN_ECP_codeml')
    #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
    parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
    summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                        codeml_result['NSsites'][0]['lnL'],
                        float(parameter_list[-1]),
                        float(parameter_list[-3]),
                        float(parameter_list[-2])]
                       )
    pair = ['EDN','ECP']
    np.savetxt(open(path+'output/summary/' + '_'.join(pair) + '.txt', 'w+'), np.array([codeml_result['NSsites'][0]['lnL'],2*(len(codeml_result['NSsites'][0]['parameters']['branches'])+5)-2*codeml_result['NSsites'][0]['lnL']]))
    finished_list.append(['EDN','ECP'])
    
    codeml_result = codeml.read('/Users/xjw1001001/Desktop/PAML/output/ERa_ERb/out/ERa_ERb_codeml')
    #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
    parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
    summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                        codeml_result['NSsites'][0]['lnL'],
                        float(parameter_list[-1]),
                        float(parameter_list[-3]),
                        float(parameter_list[-2])]
                       )
    pair = ['ERa','ERb']
    np.savetxt(open(path+'output/summary/' + '_'.join(pair) + '.txt', 'w+'), np.array([codeml_result['NSsites'][0]['lnL'],2*(len(codeml_result['NSsites'][0]['parameters']['branches'])+5)-2*codeml_result['NSsites'][0]['lnL']]))
    finished_list.append(['ERa','ERb'])
    
    codeml_result = codeml.read('/Users/xjw1001001/Desktop/PAML/output/ARa_ERa/out/ARa_ERa_codeml')
    #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
    parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
    summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                        codeml_result['NSsites'][0]['lnL'],
                        float(parameter_list[-1]),
                        float(parameter_list[-3]),
                        float(parameter_list[-2])]
                       )
    pair = ['ARa','ERa']
    np.savetxt(open(path+'output/summary/' + '_'.join(pair) + '.txt', 'w+'), np.array([codeml_result['NSsites'][0]['lnL'],2*(len(codeml_result['NSsites'][0]['parameters']['branches'])+5)-2*codeml_result['NSsites'][0]['lnL']]))
    finished_list.append(['ARa','ERa'])


    for pair in SRlist:
        codeml_cmd = [codeml_dir, path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_control.ctl']
        #workpath =path + 'output/' + '_'.join(pair) + '/out/'
        #os.chdir(workpath)
        #os.system(' '.join(codeml_cmd))
        codeml_result = codeml.read(path+'output/' + '_'.join(pair) + '/out/' + '_'.join(pair) + '_codeml')
        #baseml_result = baseml.read('/Users/xjw1001001/Documents/GitHub/IGCexpansion2/test/Ancestral_reconstruction/PAML/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
        parameter_list = codeml_result['NSsites'][0]['parameters']['parameter list'].split(' ')
        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                            codeml_result['NSsites'][0]['lnL'],
                            float(parameter_list[-1]),
                            float(parameter_list[-3]),
                            float(parameter_list[-2])]
                           )
        finished_list.append(pair)
        np.savetxt(open(path+'output/summary/' + '_'.join(pair) + '.txt', 'w+'), np.array([codeml_result['NSsites'][0]['lnL'],2*(len(codeml_result['NSsites'][0]['parameters']['branches'])+5)-2*codeml_result['NSsites'][0]['lnL']]))

    
    
    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    np.savetxt(open(path+'output/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
