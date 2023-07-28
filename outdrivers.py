import pandas as pd
import numpy as np
from gensim.models import Word2Vec
import gensim
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.cluster import AgglomerativeClustering
import pickle
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import os, sys
from multiprocessing import Process, Queue, Manager
from scipy import stats
from random import choice
from scipy.cluster.hierarchy import linkage
import importdatas
import getembed
import networkx

q = 0.2
p = 0.33

def calsimi_new(model, gene_wgs, sample_wgs, wgs_sample):
    samplelist = list(sample_wgs.keys())
    genelist = list(gene_wgs.keys())

    gene_sample_simi = pd.DataFrame(0, index=genelist, columns=samplelist)
    gene_sample_simi_count = pd.DataFrame(0, index=genelist, columns=samplelist)
    for gene in gene_wgs.keys():
        for v in gene_wgs[gene]:
            vgs = model.wv.similarity(gene,v) + 1
            for sample in wgs_sample[v]:
                simi = model.wv.similarity(gene,sample) + 1
                gene_sample_simi.loc[gene,sample] += simi * vgs
                gene_sample_simi_count.loc[gene, sample] += 1
    gene_sample_simi.to_csv(result_path+notes+'/gene_sample_simi_new.csv')
    gene_sample_simi_count.to_csv(result_path + notes + '/gene_sample_num.csv')

    gene_sample_simi = pd.read_csv(result_path + notes + '/gene_sample_simi_new.csv', index_col=0)
    gene_sample_simi_count = pd.read_csv(result_path + notes + '/gene_sample_num.csv', index_col=0)
    gene_info = pd.DataFrame(0, index=gene_sample_simi.index, columns=['RawScore','MutRate'])

    gene_info['RawScore'] = gene_sample_simi.sum(axis=1) / len(samplelist)
    gene_info['MutRate'] = gene_sample_simi_count.sum(axis=1) / len(samplelist)
    gene_info.to_csv(result_path + notes + '/gene_info.csv')

def gene_netprop(result_path,notes):
    print('getinfo')

    gene_info = pd.read_csv(result_path + notes + '/gene_info.csv', index_col=0)
    gene_info.index = [g[5:] for g in gene_info.index]

    pathwaynet = pickle.load(open(sourcepath + '/Network_pathwayn.pkl','rb'))
    ppinet = pickle.load(open(sourcepath + '/Network_ppi.pkl', 'rb'))

    networkad = q * ppinet + (1 - q) * pathwaynet
    genelist = [g for g in gene_info.index if g in networkad.index]
    networkad = networkad.loc[genelist][genelist]
    networkad = networkad.loc[~(networkad==0).all(axis=1)]
    networkad = networkad[list(networkad.index)]
    networkad = (networkad / networkad.sum()).T
    print(networkad.shape)
    genelist = list(networkad.index)
    gene_info = gene_info.loc[genelist]

    seedvec = list()
    for gene in genelist:
        seedvec.append(gene_info.loc[gene, 'MutRate']) #'RawScore''MutRate'

    seedvec = np.array(seedvec)
    vec = seedvec
    ad = np.array(networkad)

    step = 0
    vec = np.array(seedvec)
    vec = (1 - p) * np.dot(vec, ad) + p * seedvec
    while step < 5000:
        vec = (1 - p) * np.dot(vec, ad) + p * seedvec
        step += 1
    gene_info['ScoreAfterProp']=vec
    gene_info.to_csv(result_path + notes + '/gene_info_proped.csv')

def outdrivergenes():
    drivergenes = pd.read_csv(result_path + notes + '/gene_info_proped.csv', index_col=0)
    thres = np.percentile(drivergenes['RawScore'], 99)
    drivergenes = drivergenes[(drivergenes['RawScore'] > thres) & (drivergenes['ScoreAfterProp'] > thres)]
    drivergenes.to_csv(result_path + notes + '/driver_genes.csv')

if __name__ == '__main__':
    datapath, pathwaypath, ppi_path, result_path, sourcepath, \
    window, dim, WL, WN, notes, UseTraits, Traits_path, Drug, DrugPath, withV, Kernel = getembed.getargs(argfile = './arg_pan cancer.txt')

    model = gensim.models.Word2Vec.load(result_path+notes+'/trained_model.model')

    f = open(sourcepath + 'sample_wgs.pkl', 'rb')
    sample_wgs = pickle.load(f)
    f = open(sourcepath + 'var_gene.pkl', 'rb')
    var_gene = pickle.load(f)
    gene_wgs, wgs_sample = importdatas.inverse(sample_wgs, var_gene)
    calsimi_new(model, gene_wgs, sample_wgs, wgs_sample)
    gene_netprop(result_path,notes)
    outdrivergenes()