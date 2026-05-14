#!/home/zhenghoufengLab/gaisirui/miniconda3/bin/python

import os
import numpy as np
import pandas as pd
import scipy.io as scio
import torch
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

def ReadData(b):
    PHENO = pd.read_csv(b, sep = "\t| ", header = None, engine='python')
    plt.figure()
    plt.scatter(PHENO.iloc[:, 0], PHENO.iloc[:, 6], marker = "+", color = "Purple")
    plt.savefig("./SCORE" + str(b) + ".png", bbox_inches="tight")
    plt.close()
    auc_roc_sum = 0
    scaler_ss = StandardScaler()
    # 训练接操作
    new_train_x = scaler_ss.fit_transform(PHENO.iloc[:, 6].values.reshape(-1,1))
    # 测试集操作
    new_test_x = scaler_ss.fit_transform(PHENO.iloc[:, 0].values.reshape(-1,1))
    r2 = metrics.r2_score(new_train_x, new_test_x)
    return r2

with open("./all_type.answer", "a") as f:
    b = "time.spearman"
    
    if os.path.exists(b):
        ReadData(b)
