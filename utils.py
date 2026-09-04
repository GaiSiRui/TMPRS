import numpy as np
import sys
from torch.autograd import Variable
import torch.utils.data as data
import math
import time
import numpy as np
import random
import os
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from math import sqrt
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.optim.lr_scheduler import LambdaLR
from sklearn.model_selection import KFold
from scipy import stats
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_curve, precision_recall_curve, average_precision_score, roc_auc_score
from imblearn.under_sampling import ClusterCentroids
from sklearn import preprocessing
import filecmp
from itertools import islice
from pyplink import PyPlink
from sklearn import metrics

def same_fam(inputfile_list):
    with open(inputfile_list, "r") as fileHandler:
        # 读取文件中的下一行
        line = fileHandler.readline().replace('\r','').replace('\n','')
        # 检查行是否为空
        while line:
            #print(line.strip())
            old_line = line
            line = fileHandler.readline().replace('\r','').replace('\n','')

            if line:
                if not filecmp.cmp(line+".fam", old_line+".fam"):
                    print("Different!")
                    return 1

    return 0

def file_line(file, lines, names, number = 1):
    with open(file, 'r') as file_a:
        num_lines = sum(1 for line in file_a)

    with open(file, 'r') as file_a:
        i = 0
        df = pd.DataFrame(data=None,columns=names)

        for line in file_a.readlines():
            i += 1
            if i <= num_lines * number:
                #print(list(np.array(line.split()).take(lines)))
                #print(type(list(np.array(line.split()).take(lines))))
                df.loc[line.split()[1]] = list(np.array(line.split()).take(lines))

            if i % 100 == 99:
                print(i)

    return df

def max_with_source(x, y):
    # x 和 y 是元组，第一个元素是值，第二个元素是来源标识
    if x[0] > y[0]:
        return x
    else:
        return y

def workers(name, pool_list, assoc_combine, threshold, bim, fam, args, origin_PRS, the_type):
    #print(origin_PRS, "origin_PRS")
    PRS = np.zeros_like(fam.loc[:, "status"])
    for row in assoc_combine.itertuples():
        #print(row, "row")
        if row.non_additive_point * 100 > pool_list:
            if row.non_additive == "rec":
                if row.P_REC <= threshold and row.P_REC > threshold - args.interval:
                    #print(row.SNP, row.P_REC, row.non_additive_point * 100, pool_list, "rec")
                    if the_type == "test":
                        with PyPlink(args.test_file) as bed:
                            geno_data = bed.get_geno_marker(row.SNP)

                            for individual in range(len(origin_PRS)):
                                if geno_data[individual] == 2:
                                    PRS[individual] = origin_PRS[individual] + row.beta_REC * HAM(row.P_REC)
                                    #PRS[individual] = origin_PRS[individual] + row.beta_REC
                                    origin_PRS[individual] = PRS[individual]

                    else:
                        with PyPlink(args.validation_file) as bed:
                            geno_data = bed.get_geno_marker(row.SNP)

                            for individual in range(len(origin_PRS)):
                                if geno_data[individual] == 2:
                                    PRS[individual] = origin_PRS[individual] + row.beta_REC * HAM(row.P_REC)
                                    #PRS[individual] = origin_PRS[individual] + row.beta_REC
                                    origin_PRS[individual] = PRS[individual]

            if row.non_additive == "dom":
                if row.P_DOM <= threshold and row.P_DOM > threshold - args.interval:
                    #print(row.SNP, row.P_DOM, row.non_additive_point * 100, pool_list)
                    if the_type == "test":
                        with PyPlink(args.test_file) as bed:
                            geno_data = bed.get_geno_marker(row.SNP)

                            for individual in range(len(origin_PRS)):
                                if geno_data[individual] == 2 or geno_data[individual] == 1:
                                    PRS[individual] = origin_PRS[individual] + row.beta_DOM * HAM(row.P_DOM)
                                    #PRS[individual] = origin_PRS[individual] + row.beta_DOM
                                    origin_PRS[individual] = PRS[individual]

                    else:
                        with PyPlink(args.validation_file) as bed:
                            geno_data = bed.get_geno_marker(row.SNP)

                            for individual in range(len(origin_PRS)):
                                if geno_data[individual] == 2 or geno_data[individual] == 1:
                                    PRS[individual] = origin_PRS[individual] + row.beta_DOM * HAM(row.P_DOM)
                                    #PRS[individual] = origin_PRS[individual] + row.beta_DOM
                                    origin_PRS[individual] = PRS[individual]

        else:
            if row.P_ADD <= threshold and row.P_ADD > threshold - args.interval:
                if the_type == "test":
                    with PyPlink(args.test_file) as bed:
                        geno_data = bed.get_geno_marker(row.SNP)

                        for individual in range(len(origin_PRS)):
                            if geno_data[individual] == 2:
                                PRS[individual] = origin_PRS[individual] + row.beta_ADD * 2 * HAM(row.P_ADD)
                                #PRS[individual] = origin_PRS[individual] + row.beta_ADD * 2
                                origin_PRS[individual] = PRS[individual]

                            elif geno_data[individual] == 1:
                                PRS[individual] = origin_PRS[individual] + row.beta_ADD * HAM(row.P_ADD)
                                #PRS[individual] = origin_PRS[individual] + row.beta_ADD
                                origin_PRS[individual] = PRS[individual]

                else:
                    with PyPlink(args.validation_file) as bed:
                        geno_data = bed.get_geno_marker(row.SNP)

                        for individual in range(len(origin_PRS)):
                            if geno_data[individual] == 2:
                                PRS[individual] = origin_PRS[individual] + row.beta_ADD * 2 * HAM(row.P_ADD)
                                #PRS[individual] = origin_PRS[individual] + row.beta_ADD * 2
                                origin_PRS[individual] = PRS[individual]

                            elif geno_data[individual] == 1:
                                PRS[individual] = origin_PRS[individual] + row.beta_ADD * HAM(row.P_ADD)
                                #PRS[individual] = origin_PRS[individual] + row.beta_ADD
                                origin_PRS[individual] = PRS[individual]

    return PRS

def spearman(a, b):
    rho, pval = stats.spearmanr(a, b)
    return(rho)

def Pearson(b):
    PHENO = pd.read_csv(b, sep = "\t| ", header = None, engine='python')
    rho, pval = stats.pearsonr(PHENO.iloc[:, 0], PHENO.iloc[:, 6])
    return(rho)

def R2(b):
    PHENO = pd.read_csv(b, sep = "\t| ", header = None, engine='python')
    auc_roc_sum = 0
    scaler_ss = StandardScaler()
    # 训练接操作
    new_train_x = scaler_ss.fit_transform(PHENO.iloc[:, 6].values.reshape(-1,1))
    # 测试集操作
    new_test_x = scaler_ss.fit_transform(PHENO.iloc[:, 0].values.reshape(-1,1))
    r2 = metrics.r2_score(new_train_x, new_test_x)
    return r2

def bce(a, b):
    fpr, tpr, thresholds = metrics.roc_curve(a, b, pos_label = 2)
    return(metrics.auc(fpr, tpr))

def preprocess_for_roc(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    mask = ~(np.isnan(a) | np.isnan(b))
    a_clean = a[mask]
    b_clean = b[mask]
    return a_clean, b_clean

def HAM(p_value):
    z_value = stats.norm.ppf(1 - p_value)
    out = 1-(z_value/1+z_value)**2
    #print(out, "HAM", p_value, "p_value")
    return out

def mean(bfile, args):
    adjusted_pheno = pd.read_csv(args.input_file + ".pheno", delim_whitespace=True)

    # 2. 准备结果存储
    results = []

    # 3. 使用PyPlink遍历SNPs
    with PyPlink(args.input_file) as bed:  # your_data是PLINK文件前缀
        sample_pheno = adjusted_pheno.set_index(['FID', 'IID'])['PHENO']
        bim = bed.get_bim()
        fam = bed.get_fam()
        # Iterating over all loci
        for loci_name, genotypes in bed:
            pass

        sample_pheno = adjusted_pheno.set_index(['FID', 'IID'])['PHENO']
        # Getting the genotypes of a single marker (numpy.ndarray)
        for idx, (snp_id, genotypes) in enumerate(bed):
            print(type(snp_info))          # 查看类型
            print(snp_info)               # 查看内容
            print(len(snp_info))          # 查看长度
            print(type(genotypes))          # 查看类型
            print(genotypes)               # 查看内容
            print(len(genotypes))          # 查看长度
            chrom = bim.iloc[idx]['chrom']
            cm_pos = bim.iloc[idx]['cm']
            phys_pos = bim.iloc[idx]['pos']
            a1 = bim.iloc[idx]['a1']
            a2 = bim.iloc[idx]['a2']
            #break
            snp_info = bim_indexed.loc[snp_id]
            chrom = snp_info['chrom']
            cm_pos = snp_info['cm']
            phys_pos = snp_info['pos']
            a1 = snp_info['a1']
            a2 = snp_info['a2']
            #for (chrom, snp_id, cm_pos, phys_pos, a1, a2), genotypes in bed:
            df = pd.DataFrame({
                'FID': bed.fid,
                'IID': bed.iid,
                'genotype': genotypes
            })
            df = df.set_index(['FID', 'IID'])
            df['adjusted_pheno'] = sample_pheno
            df = df[df['genotype'] != -1].copy()
            df['genotype'] = df['genotype'].round().astype(int)
            means = df.groupby('genotype')['adjusted_pheno'].mean()
            results.append({
                'SNP': snp_id,
                'chrom': chrom,
                'position': phys_pos,
                'A1': a1,
                'A2': a2,
                'aa_mean': means.get(0, np.nan),
                'Aa_mean': means.get(1, np.nan),
                'AA_mean': means.get(2, np.nan),
                'N_aa': len(df[df['genotype'] == 0]),
                'N_Aa': len(df[df['genotype'] == 1]),
                'N_AA': len(df[df['genotype'] == 2])
            })
    result_df = pd.DataFrame(results)
    result_df.to_csv("snp_genotype_means.csv", index=False)
    print("分析完成，结果已保存到 snp_genotype_means.csv")
