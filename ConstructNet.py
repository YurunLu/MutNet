import pandas as pd
import numpy as np
import getembed

def construct_prop_net(pathwaypath, ppi_path):
    geneset = list(pd.read_csv(result_path + notes + '/gene_info.csv', index_col=0).index)
    geneset = [g[5:] for g in geneset]
    ppinet = pd.DataFrame(0, index=geneset, columns=geneset)
    pathwaynet = pd.DataFrame(0, index=geneset, columns=geneset)
    pathwaynetn = pd.DataFrame(0, index=geneset, columns=geneset)
    print(ppinet.shape)

    with open(ppi_path, 'r') as f:
        for line in f:
            l = line[:-1].split('\t')
            if l[1] in ppinet and l[0] in geneset:
                ppinet.loc[l[0], l[1]] = 1
                ppinet.loc[l[1], l[0]] = 1
    print('network ppi done')
    ppinet.to_csv(sourcepath + '/Network_PPI.csv')

    with open(pathwaypath + '/ReactomePathways.gmt', 'r') as f:
        for line in f:
            l = line[:-1].split('\t')
            genes = set(l[2:]) & set(geneset)
            if len(genes) > 1:
                w = 2 / (len(genes) * (len(genes) - 1))
                w1 = 1 / (len(genes) - 1)
                pathwaynet.loc[list(genes), list(genes)] += w
                pathwaynet.loc[list(genes), list(genes)] -= np.diag([w] * len(genes))
                pathwaynetn.loc[list(genes), list(genes)] += w1
                pathwaynetn.loc[list(genes), list(genes)] -= np.diag([w1] * len(genes))
    print('network pathway done')
    pathwaynet.to_csv(sourcepath + '/Network_pathway.csv')
    pathwaynetn.to_csv(sourcepath + '/Network_pathwayn.csv')

if __name__ == '__main__':
    datapath, pathwaypath, ppi_path, result_path, sourcepath, \
    window, dim, WL, WN, notes, UseTraits, Traits_path, Drug, DrugPath, withV, Kernel = getembed.getargs(
        argfile='./arg_pan cancer.txt')
    construct_prop_net(pathwaypath, ppi_path)
