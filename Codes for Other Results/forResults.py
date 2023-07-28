import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle as plk
from adjustText import adjust_text

c1 = pd.read_csv('./CancerGeneCosmic.csv', index_col=0).index
c2 = pd.read_csv('./cancerGeneList.csv', index_col=0).index
c3 = pd.read_csv('./MutPanningGeneTumorPairs_0.05.csv', index_col=0)
c3 = c3[c3['Q-value (false-discovery rate)']<0.25].index
c4 = [s.split('\t')[1].split(';')[0] for s in open('./KEGG.txt','r')]
Cancergene = set(c1)|set(c2)|set(c3)|set(c4)
print(len(Cancergene))

def getPPI():
    pathway = pd.read_csv('Results/PanCancer/example_pathway.csv', index_col=0)
    ppi = pd.read_csv('Results/PanCancer/example_ppi.csv', index_col=0)
    genes = ['FRMPD2','AGRN','EGF','HELZ2','NEB']
    for g in genes:
        print([ppig for ppig in ppi.index if ppi.loc[g,ppig]>1])
        print([pathwayg for pathwayg in pathway.index if pathway.loc[g, pathwayg] > 0.25])
    len([ppig for ppig in ppi.index if ppi.loc[ppig, genes[1]] >= 1 and ppi.loc[ppig, genes[2]] >= 1])

def drawdrivers():
    MutNetresults = pd.read_csv('./Results/PanCancer/dim128_window4_L10N10_V_T/gene_info_proped.csv', index_col=0)
    thres = np.percentile(MutNetresults['RawScore'], 90)
    drivers = list(MutNetresults[(MutNetresults['ScoreAfterProp'] > thres) & (MutNetresults['RawScore'] > thres)].index)
    print(len(set(drivers)&Cancergene))
    MutNetresults = MutNetresults.loc[drivers]
    plt.figure(figsize=(4,4))
    KCGs = list(Cancergene&set(MutNetresults.index))
    NGs = list(set(MutNetresults.index)-Cancergene)
    plt.scatter(MutNetresults.loc[NGs]['RawScore'], MutNetresults.loc[NGs]['ScoreAfterProp'], alpha=0.5, c='gray')
    plt.scatter(MutNetresults.loc[KCGs]['RawScore'], MutNetresults.loc[KCGs]['ScoreAfterProp'], alpha=0.5, c='red',
                label='Known Cancer Gene')
    texts = [plt.text(MutNetresults.loc[g,'RawScore'],MutNetresults.loc[g,'ScoreAfterProp'],g)
             for g in MutNetresults.index if
             MutNetresults.loc[g,'ScoreAfterProp']>0.105]
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$MutScore_{Raw}$')
    plt.ylabel(r'$MutScore_{AfterNP}$')
    plt.xticks([0.01,0.1,1])
    plt.yticks([0.01, 0.1, 1])
    adjust_text(texts,arrowprops=dict(arrowstyle='-',lw=0.5,color='grey'))
    plt.tight_layout()
    plt.savefig('./alldrivers.pdf')

def drawdriverwithppifrac():
    f = open('Results/PanCancer/Network_ppi.pkl', 'rb')
    ppi = plk.load(f)

if __name__=='__main__':
    #getPPI()
    drawdrivers()
    #drawdriverwithppifrac()