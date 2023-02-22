import argparse

import numpy as np
import pandas as pd

def log2fold():
    pass

def zscore():
    raise NotImplemented


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", default="log2fold")
    parser.add_argument("-exp", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-i", required=True)
    args = parser.parse_args()

    if args.m == "log2fold":
        log2foldlog2fold()

    else:
        zscore()