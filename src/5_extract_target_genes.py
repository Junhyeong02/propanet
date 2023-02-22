import os
import argparse

import networkx as nx
import numpy as np

from utils.target_gene import Target_genes

def makeGeneSet(DEGFile,geneSetF):
    ##-control
    DEGSet = set(np.genfromtxt(DEGFile, dtype=np.str))
    if geneSetF==None: 
        geneSet = DEGSet

    else: 
        geneSet = set(np.genfromtxt(geneSetF, dtype=np.str))
    
    return set(DEGSet)&geneSet

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nwk',help='network file')
    parser.add_argument('--deg')
    parser.add_argument('--tf')
    parser.add_argument('--tfrank')
    parser.add_argument('--geneSetFile')
    parser.add_argument('--out')
    args=parser.parse_args()

    if not os.path.exists(os.path.dirname(args.out)): 
        os.mkdir(os.path.dirname(args.out))

    geneSet = makeGeneSet(args.deg, args.geneSetFile)
    network = nx.read_edgelist(args.nwk,data=(('weight',float),),create_using=nx.DiGraph())
    TFSet = set(np.genfromtxt(args.tf, dtype=np.str, delimiter='\t'))

    with open(args.deg) as f1, open(args.tf) as f2:
        DEGli=f1.read().strip().split()
        TFli=f2.read().strip().split()
        edgeSet = set()

    for idx,TF in enumerate(TFli):
        edges, TGenes, TGs = Target_genes(TF,network,DEGli,TFSet,geneSet)
        edgeSet |= edges

        with open(args.out,'w') as f4:
            for (a, b, c) in edgeSet: 
                f4.write(a+'\t'+b+'\t'+str(c)+'\n')
                    
