import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib import cm

if __name__=='__main__':
    c1 = pd.read_csv('./CancerGeneCosmic.csv', index_col=0).index
    c2 = pd.read_csv('./cancerGeneList.csv', index_col=0).index
    c3 = pd.read_csv('./MutPanningGeneTumorPairs_0.05.csv', index_col=0)
    c3 = c3[c3['Q-value (false-discovery rate)'] < 0.25].index
    c4 = [s.split('\t')[1].split(';')[0] for s in open('./KEGG.txt', 'r')]
    Cancergene = set(c1) | set(c2) | set(c3) | set(c4)

    labelspecific = [g
                     for l in open('./Results/PanCancer/dim128_window4_L10N10_V_T/label_marker_gene.txt')
                     for g in l.split('\t')[1].split('|')]
    labelspecific = set(labelspecific)

    N = 8
    colors = cm.get_cmap('Set3', N)
    genes = [set([g.strip('\n') for g in open('./Results/PanCancer/dim128_window4_L10N10_V_T/drivers_community%sof8.txt'%i,'r')])
             for i in range(1,N+1)]

    drivers = pd.read_csv('./Results/PanCancer/dim128_window4_L10N10_V_T/gene_info_proped.csv',index_col=0)
    drivers['gene'] = drivers.index
    drivers['community'] = 0
    drivers['NP to Raw'] = drivers['ScoreAfterProp']/drivers['RawScore']
    for i in range(N):
        for g in genes[i]:
            drivers.loc[g,'community'] = i+1
    drivers=drivers[drivers['community']!=0]

    plt.figure(figsize=(4.5,8))

    plt.subplot(4, 1, 1)
    plt.bar(range(1,N+1),[len(s) for s in genes],color=[colors((i-0.5)/N) for i in range(1,N+1)])
    plt.xlim(xmin=0.5, xmax=8.5)
    plt.ylabel('Number of Genes')

    plt.subplot(4, 1, 2)
    plt.bar(range(1, N + 1), [len(s & Cancergene)/len(s) for s in genes], color=[colors((i - 0.5) / N) for i in range(1, N + 1)])
    plt.xlim(xmin=0.5, xmax=8.5)
    plt.ylabel('Percentage of \nKnown Cancer Genes')

    plt.subplot(4, 1, 3)
    sns.boxplot(x=drivers['community'], y=drivers['NP to Raw'], palette=[colors((i - 0.5) / N) for i in range(1, N + 1)])
    plt.yscale('log')
    plt.ylim(ymax=10)
    plt.xlim(xmin=-0.5,xmax=7.5)
    plt.xlabel('')
    plt.ylabel(r'$MutScore_{AfterNP}/MutScore_{Raw}$')

    plt.subplot(4, 1, 4)
    plt.bar(range(1, N + 1), [len(s & labelspecific) / len(s) for s in genes],
            color=[colors((i - 0.5) / N) for i in range(1, N + 1)])
    plt.xlim(xmin=0.5, xmax=8.5)
    plt.ylabel('Percentage of \nLabel Specific Genes')
    plt.xlabel('community')

    plt.tight_layout()
    plt.savefig('./Statistics for Communities.pdf')