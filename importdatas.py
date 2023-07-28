import pandas as pd


def inputwgs(datapath):
    sample_wgs = dict()
    var_gene = dict()

    muta_coding = pd.read_csv(datapath)
    for j in muta_coding.index:
        sample = muta_coding.loc[j, 'Tumor_Sample']
        if not 'sample_' + sample in sample_wgs:
            sample_wgs['sample_' + sample] = set()
        mutaID = muta_coding.loc[j, 'Chromosome'] + '_' + str(muta_coding.loc[j, 'Start_Position']) + \
                 str(muta_coding.loc[j, 'Reference_Allele'])
        sample_wgs['sample_' + sample].add('variant_MutaC_' + mutaID)
        var_gene['variant_MutaC_' + mutaID] = {'gene_' + str(muta_coding.loc[j, 'Gene'])}
    gene_wgs, wgs_sample = inverse(sample_wgs, var_gene)
    return sample_wgs, var_gene, gene_wgs, wgs_sample


def inverse(sample_wgs, var_gene):
    gene_wgs = dict()
    wgs_sample = dict()
    for sample in sample_wgs.keys():
        for variant in sample_wgs[sample]:
            if not variant in wgs_sample:
                wgs_sample[variant] = set()
            wgs_sample[variant].add(sample)
            for gene in var_gene[variant]:
                if not gene in gene_wgs:
                    gene_wgs[gene] = set()
                gene_wgs[gene].add(variant)
    return gene_wgs, wgs_sample


'''
sample_wgs: {sample:{'cnv':{'cnv_ID'}, 
                     'sv':{'sv_svID'},
                     'MutaC':{'MutaC_ID'},
                     'ncv':{'ncv_snvID','ncv_indelID'}}
                     }

'''


def importpathways(pathwaypath, geneset):
    gene_pathway = dict()
    pathways_relation = dict()
    pathway_gene = dict()
    node = dict()

    node['pathway'] = set()
    with open(pathwaypath + '/ReactomePathways.gmt', 'r') as f:
        for line in f:
            l = line[:-1].split('\t')
            pathway_gene['pathway_' + l[1]] = {'gene_' + x for x in l[2:]} & geneset
            for gene in pathway_gene['pathway_' + l[1]]:
                if not gene in gene_pathway:
                    gene_pathway[gene] = set()
                gene_pathway[gene].add('pathway_' + l[1])
                node['pathway'].add('pathway_' + l[1])

    with open(pathwaypath + '/ReactomePathwaysRelation.txt', 'r') as f:
        for line in f:
            l = line[:-1].split('\t')
            if 'pathway_' + l[1] in pathway_gene and 'pathway_' + l[0] in pathway_gene:
                if not 'pathway_' + l[1] in pathways_relation:
                    pathways_relation['pathway_' + l[1]] = set()
                pathways_relation['pathway_' + l[1]].add('pathway_' + l[0])
    node['gene'] = set(gene_pathway.keys())
    return gene_pathway, pathways_relation, pathway_gene, node

'''
gene_pathway: {gene:{pathways}}
pathways_relation: {pathway:{pathways}}
'''


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


def getpathmeans(pathwaypath):
    means = dict()
    with open(pathwaypath + '/ReactomePathways.gmt', 'r') as f:
        for line in f:
            l = line.split('\t')
            means[l[1]] = l[0]
    return means
