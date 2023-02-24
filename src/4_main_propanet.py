import pandas as pd
import numpy as np
import os
import argparse
import networkx as nx
import time
from tqdm import tqdm
import utils.propagation as NP
import utils.target_gene as TG 

def infByNode(TF, g, TFset) :
    visited = set()
    profit = 0
    stack = []
    stack.append(TF)
    while len(stack)>0:
        v = stack.pop()
        if v not in visited :
            visited.add(v)
            if v not in TFset: profit += abs(g.nodes[v]['weight'])
            for s in g.successors(v) :
                if s not in visited:
                    stack.append(s)
    visitedTGs = visited - TFset
    return profit, len(visited)-1, len(visitedTGs) 

def IM(nwk,TFset,repeat) :
    '''
    Influence Maximization
    ==============================
    Input 
        nwk (nx.object)
        TFset (set)
        repeat (int) 
    Output
        infNo (dict): reachability for each TF
        TFRank (list): sorted TF list
    '''

    infNo = {}
    for n in TFset :
        infNo[n]=0.0
    for _ in tqdm(range(repeat)):
        # Produce G'
        g = nwk.copy()
        #for (a,b) in network.edges() :
        #    if np.random.randint(1,1000)>abs(g[a][b]['weight']*1000) :
        #        g.remove_edge(a, b)
                #Calculate influcence (profit)
        for TF in TFset :
            profit, lenInf, _ = infByNode(TF, g, TFset)
            #print TF,profit,lenInf,lenInfTG
            if lenInf>0 and not np.isnan(profit/float(lenInf)) : infNo[TF] += profit/float(lenInf)
    for n in TFset:
        infNo[n]=infNo[n]/float(repeat)
    TFRank = sorted(infNo.keys(), key=lambda x: infNo[x], reverse=True)
    for key in infNo.keys() :
        if (infNo[key]==0.0) : TFRank.remove(key)
    return TFRank, infNo

def TF_adding_NP(DEGli, geneSet, TFli, TFrankFile, DEGnetFile, seed, coverNo=200, coverage=None):
    '''
    Trim TF list with improving correlation using Network propagation    
    =================================================================
    Input
        DEGli (li)
        geneSet (set)
        TFli (li) 
        TFrankFile (str)
        DEGnetFile (str)
        seed (pd.DataFrame)
        coverNo (int)
        coverage (float)           
    Output
        spearmanLi (li)
        cover (li)
        TF_trimmed (li)
        lst_node_weight (li)
    '''
    DEGnet = nx.read_edgelist(DEGnetFile,data=(('weight',float),),create_using=nx.DiGraph())
    TFset=set(TFli)
    nodeCnt=len(set(DEGnet.nodes())-TFset)
    TFli_rank=pd.read_csv(TFrankFile,sep='\t',header=None)[0].tolist()
    corr=[]
    cover=[]
    TF_trimmed=[]
    TG_set=set()
    prevCor = 0
    currTol = 0
    spearmanLi = []
    for iTF,TF in tqdm(enumerate(TFli_rank)):
        #seedFile: only selected TF( TFli[:iTF] ) is marked, otherwise 0
        TF_trimmed.append(TFli_rank[iTF])
        seed2weight = seed.loc[TF_trimmed,:].T.to_dict('records')[0]
        wk = NP.Walker(DEGnetFile, absWeight=True)
        spearman, lst_node_weight = wk.run_exp(seed2weight, TFset, 0.1, normalize=False)
        spearmanLi.append(spearman)
        corTmp=0
        corTmp=spearman[0]
        corr.append(corTmp)
        #TG_set |= TG.Target_genes(TF,DEGnet,DEGli,TFset)
        edges, TGalltmp, TGtmp = TG.Target_genes(TF,DEGnet,DEGli,TFset,geneSet)
        if TG_set == (TG_set | TGtmp): 
            TF_trimmed.pop()
            continue
        TG_set |= TGtmp

        c=float(len(TG_set))/nodeCnt
        cover.append(c)

        if prevCor > corTmp : currTol +=1
        else : currTol = 0
        if coverage != None and (cover[-1] > coverage or len(TG_set)>coverNo):
            TF_trimmed.pop()
            break
    seed2weight = seed.loc[TF_trimmed,:].T.to_dict('records')[0]
    _, lst_node_weight = wk.run_exp(seed2weight, TFset, 0.1, normalize=False)
    return spearmanLi, cover, TF_trimmed, lst_node_weight
    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('--tf',help='TF list File')
    parser.add_argument('--nwk',help='Network file')
    parser.add_argument('--lvl',help='DE level of genes')
    parser.add_argument('--deg',help='deg list file')
    parser.add_argument('--out',help='project (output directory) name')
    parser.add_argument('--imround', help="influence maximization round", default=1000)
    parser.add_argument('-p',default='10',type=int,help='# process for multiprocessing')
    parser.add_argument('-c',default='0.5',type=float,help='coverage threshold', choices=range(0, 1))
    parser.add_argument('-coverNo',default='200',type=float)
    args=parser.parse_args()

    if not os.path.exists(os.path.dirname(args.out)): 
        os.makedirs(os.path.dirname(args.out))

    with open(args.tf) as f:
        tfli = f.read().strip().split()

    with open(args.deg) as f:
        degli = set(f.read().strip().split())
   
    #networkFile --> networkx object, expFile --> pd.DataFrame
    print('Number of transcription factors', len(set(tfli)))
    print('Number of differenitial expressed genes', len(set(degli)))
    
    network = nx.read_edgelist(args.nwk, data=(('weight',float),),create_using=nx.DiGraph())

    #IM, NP for each timepoint
    print(network)

    #step0: parsing exp data, network for each timepoint
    start=time.time() 
    
    weight_tp = pd.read_csv(args.lvl, sep = "\t")
    DEG_tp = weight_tp.where(weight_tp["gene"].isin(set(degli)), other = 0)

    DEGliFile = args.deg

    ##Make nwk with DEGs in specific timepoint
    network = network.subgraph(set(degli)|set(tfli))
    DEGnetFile = args.out + ".nwk" 
    nx.write_edgelist(network,DEGnetFile,data=['weight'],delimiter='\t')
    
    ##Assign new weight to timepoint-specific network
    weight_DEG_nwk = weight_tp[weight_tp["gene"].isin(network.nodes())].set_index('gene')
    weight_dict = weight_DEG_nwk.T.to_dict('records')[0]
    nx.set_node_attributes(network,weight_dict,'weight')
    
    #step1: Influence maximization
    start=time.time() 
    ##preprocess
    nodes=pd.Series(list(set([x for x in network.nodes()])))
    TFset=set(nodes[nodes.isin(tfli)].tolist())
    ##influence maximization
    TFrank, infNo = IM(network,TFset, int(args.imround))
    ##Output into file 
    TFrankFile = args.out + '.tf.rank'

    with open(TFrankFile,'w') as TFrankF:
        for TF in TFrank:
            TFrankF.write(TF+'\n')#IM result--> cold.D2.TF_rank.1

    print('TF ranking done', time.strftime("%H:%M:%S",time.gmtime(time.time()-start)), '\n')
    
    #step2: NP based on TF rank
    start=time.time() 
    ##TF adding NP
    seed = weight_tp.set_index('gene')
    spearmanLi,coverage,TF_trimmed,lst_node_weight = TF_adding_NP(degli, degli, tfli, TFrankFile, DEGnetFile, seed, args.coverNo, coverage=args.c) 
    
    ##Save Output as file 
    TFtrimF = TFrankFile + '.trim'

    with open(TFtrimF,'w') as f:
        f.write('\n'.join(TF_trimmed))

    print(DEGliFile, TFrankFile, TFtrimF)
    print("main process end")



    
    