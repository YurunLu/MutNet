import pandas as pd
import gc
from random import choice, shuffle
import random
from gensim.models import Word2Vec
from sklearn.metrics.pairwise import cosine_similarity
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
from multiprocessing import Process, Queue, Manager
from scipy import stats
from sklearn.metrics import roc_curve, auc
import importdatas

def getppi(ppi_path, geneset):
    gene_gene = dict()
    with open(ppi_path, 'r') as f:
        for line in f:
            g1 = line[:-1].split('\t')[0]
            g2 = line[:-1].split('\t')[1]
            if 'gene_' + g1 not in geneset or 'gene_' + g2 not in geneset:
                continue
            if 'gene_' + g1 not in gene_gene:
                gene_gene['gene_' + g1] = set()
            if 'gene_' + g2 not in gene_gene:
                gene_gene['gene_' + g2] = set()
            gene_gene['gene_' + g1].add('gene_' + g2)
            gene_gene['gene_' + g2].add('gene_' + g1)
    return gene_gene


def getnodes(sample_wgs, var_gene):
    node = dict()
    node['sample'] = set(sample_wgs.keys())
    node['variant'] = set(var_gene.keys())
    print(len(node['variant']))
    return node


def getpathmeans():
    means = dict()
    with open('./pathways/ReactomePathways.gmt', 'r') as f:
        for line in f:
            l = line.split('\t')
            means[l[1]] = l[0]
    return means


# return gene_wgs, wgs_sample, pathway_gene (in dict)
def getVGPPP(v, var_gene, gene_pathway, pathways_relation):
    if withV:
        walktemp = [v]
    else:
        walktemp = []
    if len(var_gene[v]) != 0:
        gene = choice(list(var_gene[v]))
    else:
        return walktemp

    walktemp.append(gene)
    if gene in gene_pathway.keys():
        p = choice(list(gene_pathway[gene]))
        walktemp.append(p)
        while p in pathways_relation.keys():
            p = choice(list(pathways_relation[p]))
            walktemp.append(p)
    return walktemp


def walkGVSVG(walksgene, gene_wgs, wgs_sample, sample_wgs, var_gene, geneset, label):
    c = 0
    for gene in geneset:  # GVSVG
        for v1 in list(gene_wgs[gene]):
            for tempi in range(len(wgs_sample[v1])):
                tempwalk = [gene, v1]
                v = v1
                while len(tempwalk) < len(wgs_sample[v1]) * WL:  # len(wgs_sample[v])*len(gene_wgs[gene])*L:
                    s = choice(list(wgs_sample[v]))
                    tempwalk.append(s)
                    if s in sample_trait:
                        tr = choice(list(sample_trait[s]))
                        tempwalk.append(tr)
                        s = choice(list(trait_sample[tr]))
                        tempwalk.append(s)
                    v = choice(list(sample_wgs[s]))
                    tempwalk.append(v)
                    g = choice(list(var_gene[v]))
                    tempwalk.append(g)
                    v = choice(list(gene_wgs[g]))
                    tempwalk.append(v)
                v = choice(list(gene_wgs[gene]))
                s = choice(list(wgs_sample[v]))
                tempwalk = [s, v] + tempwalk
                if not withV:
                    tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
                walksgene.append(tempwalk)
        c += 1
        if c % 10 == 0:
            print('\r%s  %d / %d in SVGVS' % (label, c, len(geneset)), end="")
    print(label+' SVGVS done')
    return walksgene


def walkPGVSVGP(walkspath, pathway_gene, gene_wgs, wgs_sample, sample_wgs, var_gene, gene_pathway, pathwayset,
                label):
    temp = 0
    for pathway in pathwayset:  # PGVSVGP
        for g1 in list(pathway_gene[pathway]):
            for tempi in range(len(gene_wgs[g1])):
                if g1 in gene_wgs.keys():
                    tempwalk = [pathway, g1]
                    templen = 0
                    g = g1
                    while templen + len(tempwalk) < len(gene_wgs[g1]) * WL:
                        v = choice(list(gene_wgs[g]))
                        tempwalk.append(v)
                        s = choice(list(wgs_sample[v]))
                        tempwalk.append(s)
                        if s in sample_trait:
                            tr = choice(list(sample_trait[s]))
                            tempwalk.append(tr)
                            s = choice(list(trait_sample[tr]))
                            tempwalk.append(s)
                        v = choice(list(sample_wgs[s]))
                        tempwalk.append(v)
                        if len(var_gene[v] & set(gene_pathway.keys())) == 0:
                            g = choice(list(var_gene[v]))
                            tempwalk.append(g)
                            if not withV:
                                tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
                            walkspath.append(tempwalk)
                            templen += len(tempwalk)
                            tempwalk = [pathway, g1]
                            g = g1
                        else:
                            g = choice(list(var_gene[v] & set(gene_pathway.keys())))
                            tempwalk.append(g)
                            p = choice(list(gene_pathway[g]))
                            tempwalk.append(p)
                            g = choice(list(pathway_gene[p] & set(gene_wgs.keys())))
                            tempwalk.append(g)
                    if not withV:
                        tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
                    walkspath.append(tempwalk)
        temp += 1
        if temp % 20 == 0:
            print('\r%s  %d / %d in SVGPGVS' % (label, temp, len(pathwayset)), end="")
    print(label + ' SVGPGVS done')
    return walkspath


def walkGGG(walksgene, gene_gene, geneset, label):
    c = 0
    for gene in geneset:
        for g1 in gene_gene[gene]:
            temp = 0
            while temp < WL:
                tempwalk = [gene,g1]
                temp += 1
                ng = choice(list(gene_gene[g1]))
                while ng not in tempwalk and temp < WL:
                    tempwalk.append(ng)
                    temp += 1
                    ng = choice(list(gene_gene[ng]))
                walksgene.append(tempwalk)
        c += 1
        if c % 100 == 0:
            print('\r%s  %d / %d in GGG' % (label, c, len(geneset)), end="")
    print(label + ' GGG done')
    return walksgene


def walkSVGGVS(walksgene, gene_wgs, wgs_sample, sample_wgs, var_gene, gene_gene, geneset, label):
    # twalks = list()
    c = 0
    for gene in geneset:
        if len(gene_gene[gene] & set(gene_wgs.keys())) == 0:
            continue
        temp = 0
        tempwalk = [gene]
        g = gene
        while temp + len(tempwalk) < len(gene_gene[gene]) * WL:
            g = choice(list(gene_gene[g]))
            tempwalk.append(g)
            v = choice(list(gene_wgs[g]))
            tempwalk.append(v)
            s = choice(list(wgs_sample[v]))
            tempwalk.append(s)
            if s in sample_trait:
                tr = choice(list(sample_trait[s]))
                tempwalk.append(tr)
                s = choice(list(trait_sample[tr]))
                tempwalk.append(s)
            v = choice(list(sample_wgs[s]))
            tempwalk.append(v)
            if len(var_gene[v] & set(gene_gene.keys())) == 0:
                g = choice(list(var_gene[v]))
                tempwalk.append(g)
                temp += len(tempwalk)
                if not withV:
                    tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
                walksgene.append(tempwalk)
                tempwalk = [gene]
                g = gene
            else:
                g = choice(list(var_gene[v] & set(gene_gene.keys())))
                tempwalk.append(g)
                '''
                if len(gene_gene[g] & set(gene_wgs.keys())) == 0:
                    temp += len(tempwalk)
                    if not withV:
                        tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
                    walksgene.append(tempwalk)
                    tempwalk = [gene]
                    g = gene'''
        if not withV:
            tempwalk = list(filter(lambda x: x.split('_')[0] != 'variant', tempwalk))
        walksgene.append(tempwalk)
        c += 1
        if c % 100 == 0:
            print('\r%s  %d / %d in SVGGVS' % (label, c, len(geneset)), end="")
    print(label + ' SVGGVS done')
    return walksgene


def walkGGPPGG(walksgene, gene_gene, gene_pathway, pathway_gene, geneset, label):
    # twalks = []
    c = 0
    for gene in geneset:
        temp = 0
        if len(gene_gene[gene] & set(gene_pathway.keys())) == 0:
            continue
        while temp < len(gene_gene[gene]) * WL // 3:
            g = choice(list(gene_gene[gene] & set(gene_pathway.keys())))
            tempwalk = [gene]
            while (g not in tempwalk) and (temp + len(tempwalk) < len(gene_gene[gene]) * WL//3):
                tempwalk.append(g)
                p = choice(list(gene_pathway[g]))
                tempwalk.append(p)
                if p in pathways_relation.keys():
                    p = choice(list(pathways_relation[p]))
                    tempwalk.append(p)
                if len(pathway_gene[p] & set(gene_gene.keys())) == 0:
                    break
                else:
                    g = choice(list(pathway_gene[p] & set(gene_gene.keys())))
                    if (len(gene_gene[g] & set(gene_pathway.keys())) == 0) or (g in tempwalk):
                        break
                    tempwalk.append(g)
                    g = choice(list(gene_gene[gene] & set(gene_pathway.keys())))
            temp += len(tempwalk)
            walksgene.append(tempwalk)
        c += 1
        if c % 100 == 0:
            print('\r%s  %d / %d in GGPGG' % (label, c, len(geneset)), end="")
    print(label + ' GGPGG done')
    return walksgene


def walk(sample_wgs, var_gene, gene_wgs, wgs_sample, gene_pathway, pathways_relation, pathway_gene,traits):
    global walks, sample_trait, trait_sample
    walks = list()
    sample_trait = {'sample_'+s : {'trait_' + traits[s]} for s in traits.keys()}
    trait_sample = {'trait_'+t : set({'sample_'+s for s in traits.keys() if traits[s] == t})
                    for t in set(traits.values())}
    t = 0
    if Drug:
        walks += walkDurg(gene_wgs, wgs_sample, sample_wgs, var_gene)
        print('SVGD', len(walks))
    for sample in sample_wgs.keys():  # GVSVGPPPP
        for v1 in sample_wgs[sample]:
            temp = 0
            vtgenenum = len(var_gene[v1])
            while temp < vtgenenum * WL:
                walkv1 = getVGPPP(v1, var_gene, gene_pathway, pathways_relation)
                s1 = sample
                center = [sample]
                if sample in sample_trait:
                    tr = choice(list(sample_trait[sample]))
                    s1 = choice(list(trait_sample[tr]))
                    center = [s1, tr, sample]
                #if len(sample_wgs[s1]) != 0:
                v2 = choice(list(sample_wgs[s1]))
                walkv2 = getVGPPP(v2, var_gene, gene_pathway, pathways_relation)
                walkv2.reverse()
                #else:
                    #walkv2 = []
                walktemp = walkv2 + center + walkv1
                temp += 1
                walks.append(walktemp)
        t += 1
        if t % 20 == 0:
            print('\r%d / %d in SVGPPP'%(t, len(sample_wgs.keys())),end="")
    print('SVGPP end with', len(walks))

    walksgene = Manager().list()
    N = Kernel
    geneset = []
    for i in range(N):
        geneset.append(set())
    for i in gene_wgs.keys():
        geneset[random.randint(0, N - 1)].add(i)
    processes = []
    for i in range(N):
        processes.append(
            Process(target=walkGVSVG,
                    args=(walksgene, gene_wgs, wgs_sample, sample_wgs, var_gene, geneset[i], 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walksgene)

    print('GVSVG', len(walks))

    walkspath = Manager().list()
    N = Kernel
    pathwayset = []
    for i in range(N):
        pathwayset.append(set())
    for i in pathway_gene.keys():
        pathwayset[random.randint(0, N - 1)].add(i)
    processes = []
    for i in range(N):
        processes.append(
            Process(target=walkPGVSVGP,
                    args=(walkspath, pathway_gene, gene_wgs, wgs_sample, sample_wgs, var_gene, gene_pathway,
                          pathwayset[i], 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walkspath)
    print('PGVSVG', len(walks))

    N = Kernel
    geneset = []
    for i in range(N):
        geneset.append(set())
    for i in gene_gene.keys():
        geneset[random.randint(0, N - 1)].add(i)

    walksgene = Manager().list()
    processes = []
    for i in range(N):
        processes.append(
            Process(target=walkGGG,
                    args=(walksgene, gene_gene, geneset[i], 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walksgene)
    print(len(walks))

    walksgene = Manager().list()
    processes = []
    for i in range(N):
        processes.append(
            Process(target=walkSVGGVS,
                    args=(walksgene, gene_wgs, wgs_sample, sample_wgs, var_gene, gene_gene, geneset[i], 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walksgene)
    print(len(walks))

    walksgene = Manager().list()
    processes = []
    for i in range(N):
        processes.append(
            Process(target=walkGGPPGG,
                    args=(
                        walksgene, gene_gene, gene_pathway, pathway_gene, geneset[i], 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walksgene)

    print('PPI', len(walks))
    return walks


def walkDurg(gene_wgs, wgs_sample, sample_wgs, var_gene):
    walk = []
    gdsc = pd.read_csv(DrugPath, index_col=0)
    durg_gene = dict()
    gene_durg = dict()
    durgable = set()
    variants = set()
    for d in gdsc.index:
        # print(gdsc.loc[d,'targets'])
        genes = set(gdsc.loc[d, 'targets'].split('|'))
        genes = {'gene_' + x for x in genes}
        genes = genes & set(gene_wgs.keys())
        durgable = durgable | genes
        if len(genes) != 0:
            durg_gene['drug_' + str(d)] = genes
            for g in genes:
                variants = variants | set(gene_wgs[g])
                if g not in gene_durg:
                    gene_durg[g] = set()
                gene_durg[g].add('drug_' + str(d))
    print('Drug num', len(durg_gene))
    for s in samples:
        if len(sample_wgs[s] & variants) != 0:
            for v in sample_wgs[s] & variants:
                for g in var_gene[v] & durgable:
                    for d in gene_durg[g]:
                        for i in range(WL):
                            g1 = choice(list(durg_gene[d]))
                            v1 = choice(list(gene_wgs[g1]))
                            s1 = choice(list(wgs_sample[v1]))
                            walk.append([s, v, g, d, g1, v1, s1])

    return walk


def train(sample_wgs, var_gene, gene_wgs, wgs_sample, gene_pathway, pathways_relation, pathway_gene, dim,traits):
    walks = walk(sample_wgs, var_gene, gene_wgs, wgs_sample, gene_pathway, pathways_relation, pathway_gene,traits)
    gc.collect()
    f = open(result_path + notes + '/walks.pkl', 'wb')
    pickle.dump(walks, f)
    f.close()
    print('walk done')
    # global model
    model = Word2Vec(walks, size=dim, window=window, min_count=0, sg=1, negative=5, workers=Kernel, iter=15)
    model.save(result_path + notes + '/trained_model.model')
    model.wv.save_word2vec_format(result_path + notes + '/vecs.txt')
    del walks


def importvecs(notes, dim):
    vecs = dict()
    with open(result_path + notes + '/vecs.txt', 'r') as f:
        next(f)
        for line in f:
            line = line.split()
            vecs[' '.join(line[:-dim])] = [float(x) for x in line[-dim:]]
    return vecs


def samplecluster(samples, vecs):
    sample_info = pd.read_csv('./patients_info.csv')
    sample_info.index = sample_info['ID']
    types = {'PDAC', 'ASC', 'ACC', 'IPMN', 'SPT'}
    samplesvecs = {sample: vecs[sample] for sample in samples if
                   sample_info.loc[sample.split('_')[1], 'Patho_type'] in types}
    samples = list(samplesvecs.keys())
    print(len(samples))
    vec = pd.DataFrame(samplesvecs)
    vcorr = pd.DataFrame(0, columns=samples, index=samples)
    for i in samples:
        for j in samples:
            vcorr.loc[i, j] = cosine_similarity([vecs[i], vecs[j]])[0][1]
    print(vcorr)
    # vcorr = vec.corr()
    print('calculate corr done')

    vcorr.to_csv(result_path + notes + '/sample_simi.csv', columns=samples, index=True)

    plt.figure()
    sns.clustermap(vcorr,
                   vmax=0.8,
                   method='ward',
                   # col_colors=[p_color[sample_info.loc[x,'Patho_type']] for x in samples],
                   cmap="hot_r")

    plt.savefig(result_path + notes + '/samplecluster.jpg')


def getargs(argfile = './args_profile.txt'):
    with open(argfile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if l[0] == 'Mutation_file':
                datapath = l[1].strip('\n')
            elif l[0] == 'Pathway_Dir':
                pathwaypath = l[1].strip('\n')
            elif l[0] == 'PPI_file':
                ppi_path = l[1].strip('\n')
            elif l[0] == 'Result_Path':
                result_path = l[1].strip('\n')
            elif l[0] == 'Source_Path':
                sourcepath = l[1].strip('\n')
            elif l[0] == 'Windowize':
                window = int(l[1].strip('\n'))
            elif l[0] == 'EmbedDim':
                dim = int(l[1].strip('\n'))
            elif l[0] == 'WalkLen':
                WL = int(l[1].strip('\n'))
            elif l[0] == 'WalkNum':
                WN = int(l[1].strip('\n'))
            elif l[0] == 'Use_Triats':
                if l[1].strip('\n') in {'False', 'F', 'FALSE', 'false'}:
                    UseTraits = False
                    Traits_path = ''
                else:
                    UseTraits = True
                    Traits_path = l[1].strip('\n')
            elif l[0] == 'Use_Drugs':
                if l[1].strip('\n') in {'False', 'F', 'FALSE', 'false'}:
                    Drug = False
                    DrugPath = ''
                else:
                    Drug = True
                    DrugPath = l[1].strip('\n')
            elif l[0] == 'With_Variants':
                withV = l[1].strip('\n') in {'True', 'T', 'TRUE', 'true'}
            elif l[0] == 'Kernel':
                Kernel = int(l[1].strip('\n'))
            elif l[0] == 'Result_note':
                if l[1].strip('\n') == 'Default':
                    notes = 'dim%d_window%d_L%dN%d' % (dim, window, WL, WN)
                    if withV:
                        notes = notes + '_V'
                    if UseTraits:
                        notes = notes + '_T'
                    if Drug:
                        notes = notes + '_D'
                else:
                    notes = l[1].strip('\n')
                if not notes in os.listdir(result_path):
                    os.mkdir(result_path + notes)

    return datapath, pathwaypath, ppi_path, result_path, sourcepath, window, dim, WL, WN, notes, \
           UseTraits, Traits_path, Drug, DrugPath, withV, Kernel


if __name__ == '__main__':
    datapath, pathwaypath, ppi_path, result_path, sourcepath, window, dim, WL, WN, notes, \
    UseTraits, Traits_path, Drug, DrugPath, withV, Kernel = getargs()
    sample_wgs, var_gene, gene_wgs, wgs_sample = importdatas.inputwgs(datapath)
    wgsnodes = importdatas.getnodes(sample_wgs, var_gene)

    f = open(sourcepath + '/sample_wgs.pkl', 'wb')
    pickle.dump(sample_wgs, f)
    f.close()
    f = open(sourcepath + '/var_gene.pkl', 'wb')
    pickle.dump(var_gene, f)
    f.close()
    print('import wgs done')
    print('sample num:', len(sample_wgs))
    print('variants num:', len(var_gene))
    print('gene num:', len(gene_wgs))

    gene_pathway, pathways_relation, pathway_gene, pathwaynodes = importdatas.importpathways(pathwaypath,
                                                                                             set(gene_wgs.keys()))
    print('import pathways done')
    gene_gene = getppi(ppi_path, set(gene_wgs.keys()))
    print('import ppi done')

    f = open(sourcepath + '/nodes.pkl', 'wb')
    pickle.dump({'wgsnodes': wgsnodes, 'pathwaynodes': pathwaynodes}, f)
    f.close()
    samples = list(sample_wgs.keys())
    if UseTraits:
        traits = pd.read_csv(Traits_path, index_col=0)
        traits = {s: traits.loc[s, 'Type'] for s in traits.index
                  if traits.loc[s, 'Type'] != 'other' and 'sample_'+s in sample_wgs.keys()}
    else:
        traits = dict()

    train(sample_wgs, var_gene, gene_wgs, wgs_sample, gene_pathway, pathways_relation, pathway_gene, dim,traits)

    print(notes)
