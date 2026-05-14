#!/home/zhenghoufengLab/gaisirui/miniconda3/bin/python

import os
import numpy as np
import pandas as pd
import scipy.io as scio
import torch
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import metrics

def ReadData(b):
    PHENO = pd.read_csv(b, sep = "\t| ", header = None, engine='python')
    plt.figure()
    plt.scatter(PHENO.iloc[:, 0], PHENO.iloc[:, 6], marker = "+", color = "Purple")
    plt.savefig("./SCORE" + str(b) + ".png", bbox_inches="tight")
    plt.close()
    auc_roc_sum = 0
    fpr, tpr, thresholds = metrics.roc_curve(PHENO.iloc[:, 6], PHENO.iloc[:, 0], pos_label = 2)
    return(metrics.auc(fpr, tpr))

#with open("./all_type.answer", "a") as f:
    #b = "time.spearman"
    
    #if os.path.exists(b):
        #ReadData(b)
