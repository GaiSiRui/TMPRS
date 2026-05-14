#!/home/zhenghoufengLab/gaisirui/miniconda3/bin/python

import os
import numpy as np
import pandas as pd
import scipy.io as scio
import torch
from scipy import stats
import matplotlib.pyplot as plt

def ReadData(a, b):
    plt.figure()
    plt.scatter(a, b, marker = "+", color = "Purple")
    plt.savefig("./SCORE" + str(a) + "_" + str(b) + ".png", bbox_inches="tight")
    plt.close()
    rho, pval = stats.spearmanr(a, b)
    #print(rho)
    return(rho)

#with open("./all_type.answer", "a") as f:
    #b = "time.spearman"

    #if os.path.exists(b):
        #ReadData(b)
