

def getinfo():
    geneinfo = dict()
    with open('./9606.protein.info.v11.5.txt','r') as f:
        next(f)
        for line in f:
            l = line[:-1].split('\t')
            if l[-1]=='annotation not available' or (l[1][:4]=='ENSP' and len(l[1])==15):
                continue
            geneinfo[l[0]]=l[1]
    return geneinfo

def outlink(geneinfo):
    g = open('./ppi.txt','w')
    KRAS = set()
    with open('./9606.protein.links.v11.5.txt','r') as f:
        next(f)
        for line in f:
            l = line[:-1].split()
            if l[0] in geneinfo and l[1] in geneinfo:
                g.write(geneinfo[l[0]]+'\t'+geneinfo[l[1]]+'\t'+l[2]+'\n')
                if geneinfo[l[0]]=='KRAS' or geneinfo[l[1]]=='KRAS':
                    KRAS = KRAS | {geneinfo[l[0]],geneinfo[l[1]]}
    print(len(KRAS))
    g.close()

if __name__=='__main__':
    outlink(getinfo())