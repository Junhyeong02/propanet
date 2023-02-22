import numpy as np
import networkx as nx

def Target_genes_noTF(TF, g, TFSet) :
    visited = set()
    res=set()
    stack = []
    stack.append(TF)
    while len(stack)>0:
        v = stack.pop()
        if v not in visited :
            visited.add(v)
            if v not in TFSet:
                res.add(v)
            for s in g.successors(v) :
                if s not in visited:
                    stack.append(s)
    return res

def Target_genes(TF, g, DEGli, TFSet, geneSet) :
    visited = set()
    res=set()
    DEGset=set(DEGli)
    coverTGs = Target_genes_noTF(TF, g, TFSet)
    allEdges = False if len(coverTGs&geneSet)/float(len(geneSet))>0.5 else True
    stack = []
    stack.append(TF)
    subnet = nx.DiGraph()
    
    while len(stack)>0:
        v = stack.pop()
        if v not in visited :
            visited.add(v)
            if v in DEGset and v not in TFSet:
                res.add(v)
            for s in g.successors(v) :
                if s not in visited and (allEdges or g[v][s]['weight']>0.8):
                    subnet.add_edge(v, s, weight=g[v][s]['weight'])
                    stack.append(s)

    # include non-DE TFs in a shortest path w/ the highest weight
    edgeSet = set()
    midTFs = set()
    for v in res :
        currPath = []
        currPathScore = -np.inf
        for p in nx.all_shortest_paths(subnet, source=TF, target=v) :
            pscore = 0
            for i in range(len(p)-1) :
                pscore += g[p[i]][p[i+1]]['weight']
            pscore = pscore / float(len(p)-1)
            if pscore > currPathScore : 
                currPathScore = pscore
                currPath = p
        midTFs = midTFs | set(currPath)
        for i in range(len(currPath)-1):
            edgeSet.add((currPath[i], currPath[i+1], subnet[currPath[i]][currPath[i+1]]['weight']))
    return edgeSet, (res|midTFs)-set([TF]), res