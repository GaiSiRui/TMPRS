#!/home/zhenghoufengLab/gaisirui/miniconda3/bin/python

import os
import numpy as np
import pandas as pd
import scipy.io as scio
import torch
from scipy import stats
import matplotlib.pyplot as plt

def ReadData(b):
    PHENO = pd.read_csv(b, sep = "\t| ", header = None, engine='python')
    plt.figure()
    plt.scatter(PHENO.iloc[:, 0], PHENO.iloc[:, 6], marker = "+", color = "Purple")
    plt.savefig("./SCORE" + str(b) + ".png", bbox_inches="tight")
    plt.close()
    rho, pval = stats.pearsonr(PHENO.iloc[:, 0], PHENO.iloc[:, 6])
    return(rho)

#with open("./all_type.answer", "a") as f:
    #b = "time.pearson"
    
    #if os.path.exists(b):
        #ReadData(b)
