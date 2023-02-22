import argparse
import sys
import os

import pandas as pd
import numpy as np

from tqdm import tqdm

from scipy.stats import pearsonr
from multiprocessing import Pool, Manager

def calculate_corr(x):
    g1, g2 = x

    try:
        exp1, exp2 = gene_exp.loc[g1, :], gene_exp.loc[g2, :]
    except KeyError:
        print(g1, g2)

    corr, _ = corr_function(exp1, exp2)
    
    q.put((g1, g2, corr))
    return  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='python %(prog)s -nwk nwkFile -exp expFile -o out')
    parser.add_argument('-nwk', required=True, help='Network file')
    parser.add_argument('-exp', required=True, help='gene expression File')
    parser.add_argument('-o', required=True, help = "output file")
    parser.add_argument('-p', required=True, help="number of threads")
    args = parser.parse_args()

    with open(args.nwk) as f:
        edge_list = list(map(lambda x: x.strip().split(), f.readlines()))

    m = Manager()
    q = m.Queue()
    
    gene_exp = pd.read_csv(args.exp, sep = "\t", index_col=0)
    corr_function = pearsonr

    print("calculate correlation coefficent...")
    print(len(edge_list))

    with Pool(processes=int(args.p)) as pool:
        with tqdm(total=len(edge_list)) as pbar:
            for _ in pool.imap_unordered(calculate_corr, edge_list):
                pbar.update()

    with open(args.o, 'w') as fw:
        while not q.empty():
            start, end, corr = q.get()
            fw.write('\t'.join([start, end, str(corr)])+'\n')

    print('end...')
