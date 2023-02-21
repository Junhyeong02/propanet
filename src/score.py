import argparse
from typing import List, Any

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-true', required= True)
    parser.add_argument('-tgt', required = True)

    args = parser.parse_args()

    with open(args.true) as f:
        true_gene_set = set(gene.strip() for gene in f.readlines())
    
    with open(args.tgt) as f:
        pred_gene_set = set(gene.strip() for gene in f.readlines())

    y_true = []
    y_pred = []

    for true_gene in true_gene_set:
        y_true.append(1)
        
        if true_gene in pred_gene_set:
            y_pred.append(1)
            pred_gene_set.remove(true_gene)

        else:
            y_pred.append(0)

    for pred_gene in pred_gene_set:
        y_true.append(0)
        y_pred.append(1)

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    score = f1_score(y_true, y_pred)

    print(score)

    

    