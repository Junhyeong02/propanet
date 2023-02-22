import pandas as pd
import numpy as np
import os
import argparse
import sys
import networkx as nx
import time
import utils.propagation as NP
import utils.target_gene as TG 
import multiprocessing
from glob import glob
from tqdm import tqdm

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
    nodes=pd.Series(list(set([x for x in nwk.nodes()])))
    infNo = {}
    for n in TFset :
        infNo[n]=0.0
    for i in range(repeat) :
        # Produce G'
        g = nwk.copy()
        #for (a,b) in network.edges() :
        #    if np.random.randint(1,1000)>abs(g[a][b]['weight']*1000) :
        #        g.remove_edge(a, b)
                #Calculate influcence (profit)
        for TF in TFset :
            profit, lenInf, lenInfTG = infByNode(TF, g, TFset)
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
    corPre=-np.inf
    corr=[]
    cover=[]
    TF_trimmed=[]
    TG_set=set()
    prevCor = 0
    currTol = 0
    spearmanLi = []
    for iTF,TF in enumerate(TFli_rank):
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
        # print '-----TF\t',TF_trimmed[-1]
        # print 'correlation\t%.3f' % corTmp
        # print 'len TG_set\t', len(TG_set)
        c=float(len(TG_set))/nodeCnt
        # print 'coverage\t%.3f' % c
        cover.append(c)
        if prevCor > corTmp : currTol +=1
        else : currTol = 0
        if coverage != None and (cover[-1] > coverage or len(TG_set)>coverNo):
            TF_trimmed.pop()
            break
        corPre = corTmp
    seed2weight = seed.loc[TF_trimmed,:].T.to_dict('records')[0]
    _, lst_node_weight = wk.run_exp(seed2weight, TFset, 0.1, normalize=False)
    return spearmanLi, cover, TF_trimmed, lst_node_weight
    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('-TFliFile',help='TF list File')
    parser.add_argument('-nwkFile',help='Network file')
    parser.add_argument('-expFile',help='gene expression file with z-value,seed file for NP')
    parser.add_argument('-binFile',help='binary file of condition from DEGlimma')
    parser.add_argument('-geneSet',help='gene set by user') 
    parser.add_argument('-cond',required=True,help='condition name') 
    parser.add_argument('-outD',help='project (output directory) name')
    parser.add_argument('-p',default='10',type=int,help='# process for multiprocessing')
    parser.add_argument('-c',default='0.5',type=float,help='coverage threshold', choices=range(0, 1))
    parser.add_argument('-coverNo',default='200',type=float)
    parser.add_argument('-N', required=True, help= "number of samples", type = int)
    args=parser.parse_args()

    if args.outD != None:
        if not os.path.exists(args.outD) : os.mkdir(args.outD)
    exp = pd.read_csv(args.expFile, sep='\t',index_col=0)
    with open(args.TFliFile) as tfF:
        TFli=tfF.read().strip().split()
   
    #networkFile --> networkx object, expFile --> pd.DataFrame
    print('TF length', len(set(TFli)))
    
    n = args.N
    assert n == len(exp.columns)
    
    networks = []
    for nwkFile in tqdm(glob(args.nwkFile + "*.tsv")):
        network = nx.read_edgelist(nwkFile,data=(('weight',float),),create_using=nx.DiGraph()) 
        networks.append(network)

    print(networks)

    bins = pd.read_csv(args.binFile, sep='\t',index_col=0)
    weight = (exp*(bins.applymap(lambda x:abs(x))))
    #IM, NP for each timepoint

    def main_fxn(i):

        network = networks[i-1]
        print(network)

        #step0: parsing exp data, network for each timepoint
        start=time.time() 
        weight_tp = exp.iloc[:,i-1].reset_index()
        weight_tp.columns = ['gene','weight']

        DEG_tp =  weight.iloc[:,i-1].reset_index()
        DEG_tp.columns = ['gene','weight']
        DEGli = DEG_tp[DEG_tp['weight']!=0]['gene'].tolist()
         
        DEGliFile =args.outD+'/'+ args.cond+'.DEGli.t'+str(i)
        
        with open(DEGliFile,'w') as f:
            f.write('\n'.join(DEGli))
        if args.geneSet != None:
            geneSet = set(np.genfromtxt(args.geneSet, dtype=np.str))&DEGli
        else:
            geneSet = set(DEGli)

        ##Make nwk with DEGs in specific timepoint

        nwk_tp = network.subgraph(geneSet|set(TFli))
        DEGnetFile = args.outD + '/' + args.cond + '.nwk.t'+ str(i) 
        
        nx.write_edgelist(nwk_tp,DEGnetFile,data=['weight'],delimiter='\t')
        ##Assign new weight to timepoint-specific network
        weight_DEG_nwk = weight_tp[weight_tp.gene.isin(nwk_tp.nodes())].set_index('gene')
        weight_dict = weight_DEG_nwk.T.to_dict('records')[0]
        nx.set_node_attributes(nwk_tp,weight_dict,'weight')
        
        #step1: Influence maximization
        start=time.time() 
        ##preprocess
        nodes=pd.Series(list(set([x for x in nwk_tp.nodes()])))
        TFset=set(nodes[nodes.isin(TFli)].tolist())
        ##influence maximization
        TFrank, infNo = IM(nwk_tp,TFset,1)
        ##Output into file 
        TFrankFile =args.outD+'/'+ args.cond+'.TF_rank.t'+str(i)  


        with open(TFrankFile,'w') as TFrankF:
            for TF in TFrank:
                TFrankF.write(TF+'\n')#IM result--> cold.D2.TF_rank.1

        print(sys.stderr,'time point', i,'TF ranking done', time.strftime("%H:%M:%S",time.gmtime(time.time()-start)), '\n')
    
        
        #step2: NP based on TF rank
        start=time.time() 
        ##TF adding NP
        seed = weight_tp.set_index('gene')
        spearmanLi,coverage,TF_trimmed,lst_node_weight = TF_adding_NP(DEGli, geneSet, TFli, TFrankFile, DEGnetFile, seed, args.coverNo, coverage=args.c) 
        
        ##Save Output as file 
        TFtrimF = TFrankFile + '.trim'

        with open(TFtrimF,'w') as f:
            f.write('\n'.join(TF_trimmed))
        
        dic_node2weights = {}
        column_name = exp.columns.tolist()
        lst_seed = exp.index.values.tolist()
        wk = NP.Walker(DEGnetFile, absWeight=True)
        set_nodes = set(wk.OG.nodes())
        set_tmpNodes = set()
        for node, weight_np, all_weight in lst_node_weight:
            if node not in dic_node2weights:
                dic_node2weights[node]=[]
            dic_node2weights[node].append(weight_np)
            set_tmpNodes.add(node)
        
        for node in set_nodes-set_tmpNodes:
            if node not in dic_node2weights:
                dic_node2weights[node]=[]
            dic_node2weights[node].append(0.0)

        print(DEGliFile, TFrankFile, TFtrimF)

        return 0
    
    """
    pool = multiprocessing.Pool(args.p)
    nums=range(1,n+1)
    res = pool.imap_unordered(main_fxn,nums)
    
    time.sleep(10)

    pool.close()
    pool.join()    
    """

    for i in tqdm(range(1, n+1)):
        main_fxn(i)

    print("main process end")



    
    