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
import utils
from pyplink import PyPlink
from multiprocessing import Pool
import statsmodels.api as sm
import csv
import gc

def main(args):
    assoc_combine = before_summary(args)
    #assoc_combine = pd.read_csv('example_big.csv', sep="\t")
    print(assoc_combine)
    utils.mean(args.input_file)
    after_summary(args, assoc_combine)

def before_summary(args):
    file_name = args.test_file + ".summary_and_bim"

    with open(file_name,'a+',encoding='utf-8') as file:
        file.truncate(0)

    file_name = args.input_file + ".fam"

    if args.test_file == "False" and args.validation_file != "False":
        fam_num = utils.file_line(file_name, [0, 1], ['FID', 'IID'], 9 / 10)

    elif args.test_file != "False" and args.validation_file != "False":
        print(file_name)
        fam_num = utils.file_line(file_name, [0, 1], ['FID', 'IID'])

    elif args.test_file != "False" and args.validation_file == "False":
        fam_num = utils.file_line(file_name, [0, 1], ['FID', 'IID'], 9 / 10)

    else:
        fam_num = utils.file_line(file_name, [0, 1], ['FID', 'IID'], 4 / 5)


    file=args.input_file + '.pheno'
    utils.file_line(file_name, [1, 5], ['IID', 'PHENO']).to_csv(file, sep = "\t")

    with open(file, 'r+') as file:
        content = file.read()  # 读取文件内容
        file.seek(0, 0)  # 将文件指针移动到文件开头
        file.write('FID' + content)

    file = args.input_file + ".sex"
    utils.file_line(file_name, [1, 4], ['IID', 'SEX']).to_csv(file, sep = "\t")

    with open(file, 'r+') as file:
        content = file.read()  # 读取文件内容
        file.seek(0, 0)  # 将文件指针移动到文件开头
        file.write('FID' + content)

    file=args.input_file + '.list'
    utils.file_line(file_name, [1], ['IID']).to_csv(file, sep = "\t")

    with open(file, 'r+') as file:
            content = file.read()  # 读取文件内容
            file.seek(0, 0)  # 将文件指针移动到文件开头
            file.write('FID' + content)

    if args.beta == True:
        print("492.1")
        if args.calculate_mod == "origin":
            """
            os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic hide-covar --covar args.covar keep-pheno-on-missing-cov")
            print("492")
            """
            assoc_logistic = pd.read_csv(args.input_file + ".assoc.logistic", sep=r'\s+')
            assoc_add = assoc_logistic[assoc_logistic['TEST'] == 'ADD']
            file = args.input_file + ".assoc.add"
            assoc_add.to_csv(file, sep = "\t", index=False)

            if args.clump == False:
                    os.system("plink --bfile " + args.input_file + " --keep " + args.input_file + ".list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump " + args.input_file + ".assoc.add --clump-snp-field SNP --clump-field P")

            plink_clumped = pd.read_csv("plink.clumped", sep=r'\s+')
            assoc_add = assoc_add[assoc_add['P'] < 0.5]

            # 获取筛选后的 SNP 列表
            clumped_snps = plink_clumped['SNP'].tolist()

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc = assoc_add[assoc_add['SNP'].isin(clumped_snps)]
            snps = assoc.iloc[1, :]

        else:
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic dominant keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic dominant hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_logistic = pd.read_csv(args.input_file + ".assoc.logistic", sep=r'\s+')
            assoc_dom = assoc_logistic[assoc_logistic['TEST'] == 'DOM']
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic recessive keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic recesive hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_logistic = pd.read_csv(args.input_file + ".assoc.logistic", sep=r'\s+')
            assoc_rec = assoc_logistic[assoc_logistic['TEST'] == 'REC']
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --logistic hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_logistic = pd.read_csv(args.input_file + ".assoc.logistic", sep=r'\s+')
            assoc_add = assoc_logistic[assoc_logistic['TEST'] == 'ADD']
            assoc_logistic = pd.concat([assoc_dom, assoc_rec, assoc_add])

            if args.clump == False:
                    os.system("plink --bfile " + args.input_file + " --keep " + args.input_file + ".list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump " + args.input_file + ".assoc.logistic  --clump-snp-field SNP --clump-field P")

            plink_clumped = pd.read_csv("plink.clumped", sep=r'\s+')
            clumped_snps = plink_clumped['SNP'].tolist()
            assoc_add = assoc_logistic[assoc_logistic['TEST'] == 'ADD']
            assoc_add = assoc_add[assoc_add['P'] != "NA"]

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc_add = assoc_add[assoc_add['SNP'].isin(clumped_snps)]
            snps_add = assoc_add.iloc[1, :]
            assoc_add = assoc_add.reset_index(drop=True)

            assoc_dom = assoc_logistic[assoc_logistic['TEST'] == 'DOM']
            assoc_dom = assoc_dom[assoc_dom['SNP'].isin(clumped_snps)]
            snps_dom = assoc_dom.iloc[1, :]
            assoc_dom = assoc_dom.reset_index(drop=True)

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc_rec = assoc_logistic[assoc_logistic['TEST'] == 'REC']
            assoc_rec = assoc_rec[assoc_rec['SNP'].isin(clumped_snps)]
            snps_rec = assoc_rec.iloc[1, :]
            assoc_rec = assoc_rec.reset_index(drop=True)


        print(assoc_add['OR'])
        try:
            or_value = float(assoc_add['OR'])          # 尝试转换为浮点数
            assoc_add['beta'] = np.log(or_value)       # 转换成功则取对数
        except (ValueError, TypeError):                # 转换失败（非数值类型或无法解析的字符串）
            assoc_add['beta'] = assoc_add['OR']        # 保留原始值，不进行 log 操作
        print(assoc_dom['OR'])
        try:
            or_value = float(assoc_dom['OR'])          # 尝试转换为浮点数
            assoc_dom['beta'] = np.log(or_value)       # 转换成功则取对数
        except (ValueError, TypeError):                # 转换失败（非数值类型或无法解析的字符串）
            assoc_dom['beta'] = assoc_dom['OR']        # 保留原始值，不进行 log 操作
        print(assoc_rec['OR'])
        try:
            or_value = float(assoc_rec['OR'])          # 尝试转换为浮点数
            assoc_rec['beta'] = np.log(or_value)       # 转换成功则取对数
        except (ValueError, TypeError):                # 转换失败（非数值类型或无法解析的字符串）
            assoc_rec['beta'] = assoc_rec['OR']        # 保留原始值，不进行 log 操作

    else:
        if args.calculate_mod == "origin":
            """
            os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear hide-covar --covar args.covar keep-pheno-on-missing-cov")
            print("492")
            """
            assoc_linear = pd.read_csv(args.input_file + ".assoc.linear", sep=r'\s+')
            assoc_add = assoc_linear[assoc_linear['TEST'] == 'ADD']
            file = args.input_file + ".assoc.add"
            assoc_add.to_csv(file, sep = "\t", index=False)

            if args.clump == False:
                    os.system("plink --bfile " + args.input_file + " --keep " + args.input_file + ".list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump " + args.input_file + ".assoc.add --clump-snp-field SNP --clump-field P")

            plink_clumped = pd.read_csv("plink.clumped", sep=r'\s+')
            assoc_add = assoc_add[assoc_add['P'] < 0.5]

            # 获取筛选后的 SNP 列表
            clumped_snps = plink_clumped['SNP'].tolist()

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc = assoc_add[assoc_add['SNP'].isin(clumped_snps)]
            snps = assoc.iloc[1, :]

        else:
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear dominant keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear dominant hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_linear = pd.read_csv(args.input_file + ".assoc.linear", sep=r'\s+')
            assoc_dom = assoc_linear[assoc_linear['TEST'] == 'DOM']
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear recessive hide-covar --covar args.covar keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear recessive hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_linear = pd.read_csv(args.input_file + ".assoc.linear", sep=r'\s+')
            assoc_rec = assoc_linear[assoc_linear['TEST'] == 'REC']
            if args.covar == False:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear keep-pheno-on-missing-cov")
            else:
                os.system("plink --bfile " + args.input_file + "  --pheno " + args.input_file + ".pheno --keep " + args.input_file + ".list --geno 0.05 --maf 0.05 --out " + args.input_file + " --linear hide-covar --covar " + args.covar + " keep-pheno-on-missing-cov")
            assoc_linear = pd.read_csv(args.input_file + ".assoc.linear", sep=r'\s+')
            assoc_add = assoc_linear[assoc_linear['TEST'] == 'ADD']
            assoc_linear = pd.concat([assoc_dom, assoc_rec, assoc_add])

            if args.clump == False:
                os.system("plink --bfile " + args.input_file + " --keep " + args.input_file + ".list --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump " + args.input_file + ".assoc.linear  --clump-snp-field SNP --clump-field P")

            plink_clumped = pd.read_csv("plink.clumped", sep=r'\s+')
            clumped_snps = plink_clumped['SNP'].tolist()
            assoc_add = assoc_linear[assoc_linear['TEST'] == 'ADD']
            assoc_add = assoc_add[assoc_add['P'] != "NA"]

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc_add = assoc_add[assoc_add['SNP'].isin(clumped_snps)]
            snps_add = assoc_add.iloc[1, :]

            assoc_dom = assoc_linear[assoc_linear['TEST'] == 'DOM']
            assoc_dom = assoc_dom[assoc_dom['SNP'].isin(clumped_snps)]
            snps_dom = assoc_dom.iloc[1, :]
            assoc_dom.reset_index(drop=True)

            # 筛选 assoc_add 中 SNP 列在 valid_snps 列表中的行
            assoc_rec = assoc_linear[assoc_linear['TEST'] == 'REC']
            assoc_rec = assoc_rec[assoc_rec['SNP'].isin(clumped_snps)]
            snps_rec = assoc_rec.iloc[1, :]
            assoc_rec.reset_index(drop=True)


    if args.calculate_mod == "t":
        try:
            log_dom = np.log(assoc_dom.iloc[:, 7])
        except TypeError:
            log_dom = None

        try:
            log_rec = np.log(assoc_rec.iloc[:, 7])
        except TypeError:
            log_rec = None

        if log_dom is not None and log_rec is not None:
        # 两列均为数值，执行原始逻辑
            max_abs = log_dom.abs().where(log_dom.abs() > log_rec.abs(), log_rec.abs())
        elif log_dom is not None:
            # 仅 dom 有效
            max_abs = log_dom.abs()
        elif log_rec is not None:
            # 仅 rec 有效
            max_abs = log_rec.abs()
        else:
            # 两列均为字符串
            max_abs = pd.Series(["NA"] * len(assoc_dom.iloc[:, 7]), index=assoc_dom.iloc[:, 7].index)

        # 计算来源 DataFrame
        dom_vals = pd.to_numeric(assoc_dom.iloc[:, 7], errors='coerce')
        rec_vals = pd.to_numeric(assoc_rec.iloc[:, 7], errors='coerce')

        # 计算对数绝对值，非数值或非正数会得到 NaN 或 Inf（用 isfinite 过滤）
        log_dom_abs = np.abs(np.log(dom_vals))
        log_rec_abs = np.abs(np.log(rec_vals))

        # 判断有效值（有限且非 NaN）
        dom_valid = np.isfinite(log_dom_abs)
        rec_valid = np.isfinite(log_rec_abs)

        # 初始化结果，默认 'dom'（对应两列都无效的情况）
        max_sources = pd.Series('dom', index=assoc_dom.index)

        # 两者都有效 → 比较绝对值大小
        both_valid = dom_valid & rec_valid
        max_sources[both_valid] = np.where(log_dom_abs[both_valid] > log_rec_abs[both_valid], 'dom', 'rec')

        # 仅 rec 有效 → 输出 'rec'
        only_rec = (~dom_valid) & rec_valid
        max_sources[only_rec] = 'rec'

        # 仅 dom 有效或都无效 → 默认 'dom' 已满足
        # 若需要 numpy 数组，可转为 .values
        max_sources = max_sources.values

        # 2. 计算 A：较大者除以 add 中第七列的绝对值
        add_vals = pd.to_numeric(assoc_add.iloc[:, 7], errors='coerce')
        # 计算对数绝对值（负数或 NaN 会产生 NaN，+inf 亦处理）
        log_add_abs = np.abs(np.log(add_vals))
        valid_add = np.isfinite(log_add_abs)          # 判断对数绝对值是否有限且非 NaN

        # 将 max_abs 也转为数值（可能包含 "NA" 字符串）
        max_abs_num = pd.to_numeric(max_abs, errors='coerce')
        valid_max = np.isfinite(max_abs_num)

        # 两者均有效时执行除法，否则输出 "NA"
        A = pd.Series("NA", index=max_abs.index)      # 默认全部为 "NA"
        both_valid = valid_add & valid_max
        A[both_valid] = max_abs_num[both_valid] / log_add_abs[both_valid]
        # 3. 创建新的 DataFrame summary
        summary = pd.DataFrame({
                'SNP': assoc_add.iloc[:, 1],  # add 的第二列
                'non_additive_point': A,                       # 计算得到的 A
                "non_additive": max_sources
        })

    if args.calculate_mod == "beta":
        try:
            log_dom = np.log(assoc_dom.iloc[:, 6])
        except TypeError:
            log_dom = None

        try:
            log_rec = np.log(assoc_rec.iloc[:, 6])
        except TypeError:
            log_rec = None

        if log_dom is not None and log_rec is not None:
        # 两列均为数值，执行原始逻辑
            max_abs = log_dom.abs().where(log_dom.abs() > log_rec.abs(), log_rec.abs())
        elif log_dom is not None:
            # 仅 dom 有效
            max_abs = log_dom.abs()
        elif log_rec is not None:
            # 仅 rec 有效
            max_abs = log_rec.abs()
        else:
            # 两列均为字符串
            max_abs = pd.Series(["NA"] * len(assoc_dom.iloc[:, 6]), index=assoc_dom.iloc[:, 6].index)

        # 计算来源 DataFrame
        dom_vals = pd.to_numeric(assoc_dom.iloc[:, 6], errors='coerce')
        rec_vals = pd.to_numeric(assoc_rec.iloc[:, 6], errors='coerce')

        # 计算对数绝对值，非数值或非正数会得到 NaN 或 Inf（用 isfinite 过滤）
        log_dom_abs = np.abs(np.log(dom_vals))
        log_rec_abs = np.abs(np.log(rec_vals))

        # 判断有效值（有限且非 NaN）
        dom_valid = np.isfinite(log_dom_abs)
        rec_valid = np.isfinite(log_rec_abs)

        # 初始化结果，默认 'dom'（对应两列都无效的情况）
        max_sources = pd.Series('dom', index=assoc_dom.index)

        # 两者都有效 → 比较绝对值大小
        both_valid = dom_valid & rec_valid
        max_sources[both_valid] = np.where(log_dom_abs[both_valid] > log_rec_abs[both_valid], 'dom', 'rec')

        # 仅 rec 有效 → 输出 'rec'
        only_rec = (~dom_valid) & rec_valid
        max_sources[only_rec] = 'rec'

        # 仅 dom 有效或都无效 → 默认 'dom' 已满足
        # 若需要 numpy 数组，可转为 .values
        max_sources = max_sources.values

        # 2. 计算 A：较大者除以 add 中第七列的绝对值
        add_vals = pd.to_numeric(assoc_add.iloc[:, 6], errors='coerce')
        # 计算对数绝对值（负数或 NaN 会产生 NaN，+inf 亦处理）
        log_add_abs = np.abs(np.log(add_vals))
        valid_add = np.isfinite(log_add_abs)          # 判断对数绝对值是否有限且非 NaN

        # 将 max_abs 也转为数值（可能包含 "NA" 字符串）
        max_abs_num = pd.to_numeric(max_abs, errors='coerce')
        valid_max = np.isfinite(max_abs_num)

        # 两者均有效时执行除法，否则输出 "NA"
        A = pd.Series("NA", index=max_abs.index)      # 默认全部为 "NA"
        both_valid = valid_add & valid_max
        A[both_valid] = max_abs_num[both_valid] / log_add_abs[both_valid]
        summary = pd.DataFrame({
                'SNP': assoc_add.iloc[:, 1],  # add 的第二列
                'non_additive_point': A,                       # 计算得到的 A
                "non_additive": max_sources
        })

    if args.calculate_mod == "mean":
        pheno = pd.read_csv(args.pheno, sep="\s+")
        if args.covar != False:
            covar = pd.read_csv(args.covar, sep="\s+")
        geno_name = pd.read_csv(args.input_file + ".bim", header = None , sep="\s+", usecols=[1])
        os.system("plink --bfile " + args.input_file + " --recode --out " + args.input_file)
        
        with open('mean.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['rs_name', 'aa_mean_adj', 'Aa_mean_adj', 'AA_mean_adj'])
            file.close()
        
        for rs_number, rs_name in geno_name.iterrows():
            if args.covar == False:
                data = pheno
            else:
                data = pd.merge(pheno, covar, on=["FID", "IID"])
            target_snp = rs_name
            columns_to_load = ["FID", "IID", target_snp.iloc[0]]
        
            try:
                geno = pd.read_csv(args.input_file + ".ped", sep="\s+", usecols=columns_to_load)
        
            except Exception as e:
                print("error", e)
        
            data = pd.merge(data, geno[["FID", "IID", rs_name.iloc[0]]], on=["FID", "IID"])
        
            X = data["sex"].values.reshape(-1, 1)
            X = sm.add_constant(X)
            y = data["PHENO"]
            model = sm.OLS(y, X).fit()
            data["RESIDUAL"] = model.resid
        
            aa_mean_adj = data.loc[np.isclose(data[rs_name], 0.0), 'RESIDUAL'].mean() + y.mean()
            Aa_mean_adj = data.loc[np.isclose(data[rs_name], 1.0), 'RESIDUAL'].mean() + y.mean()
            AA_mean_adj = data.loc[np.isclose(data[rs_name], 2.0), 'RESIDUAL'].mean() + y.mean()
        
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')

        # 2. 计算 A：较大者除以 add 中第七列的绝对值
        max_sources = np.where(Aa_mean_adj > (AA_mean_adj + aa_mean_adj) / 2, 'dom', 'rec')
        A = (Aa_mean_adj - (AA_mean_adj + aa_mean_adj) / 2).abs() / ((AA_mean_adj - aa_mean_adj) / 2).abs()

        # 3. 创建新的 DataFrame summary
        summary = pd.DataFrame({
                'SNP': assoc_add.iloc[:, 1],  # add 的第二列
                'non_additive_point': A,                       # 计算得到的 A
                "non_additive": max_sources
        })


    assoc_add = assoc_add[['SNP', 'A1', 'STAT', 'P', "beta"]]

    # 合并 A 和 B 的选定列，基于 A 的第一列和 B 的第二列
    assoc_combine_1 = pd.merge(summary, assoc_add, left_on='SNP', right_on='SNP')
    assoc_combine_1.rename(columns={'STAT': 'STAT_ADD', 'P': 'P_ADD', "beta": "beta_ADD"}, inplace=True)
    assoc_dom = assoc_dom[['SNP', 'STAT', 'P', "beta"]]
    assoc_combine_2 = pd.merge(assoc_combine_1, assoc_dom, left_on='SNP', right_on='SNP')
    assoc_combine_2.rename(columns={'STAT': 'STAT_DOM', 'P': 'P_DOM', "beta": "beta_DOM"}, inplace=True)
    assoc_rec = assoc_rec[['SNP', 'STAT', 'P', "beta"]]
    assoc_combine = pd.merge(assoc_combine_2, assoc_rec, left_on='SNP', right_on='SNP')
    assoc_combine.rename(columns={'STAT': 'STAT_REC', 'P': 'P_REC', "beta": "beta_REC"}, inplace=True)

    # 删除多余的 B_col2 列（如果需要）
    #assoc_combine_1 = assoc_combine_1.drop(columns=['SNP'])

    return(assoc_combine)



def after_summary(args, assoc_combine):
    pedfile = PyPlink(args.input_file)

    with PyPlink(args.test_file) as bed:
        with PyPlink(args.validation_file) as bed_validation:
            bim = bed.get_bim()
            fam = bed.get_fam()
            bim_validation = bed_validation.get_bim()
            fam_validation = bed_validation.get_fam()

            for loci_name, genotypes in bed:
                for loci_name, genotypes in bed_validation:
                    pass

            #pool_list = []
            #assoc_combine.to_csv("example.csv", index=False, encoding="utf-8", sep="\t")
            #bim.to_csv("example_bim.csv", index=False, encoding="utf-8", sep="\t")
            origin_PRS = []
            origin_PRS_validation = []
            for k in np.arange(args.begin_site, args.end_site, args.site_step):
                origin_PRS.append(np.zeros_like(fam.loc[:, "status"]))
                origin_PRS_validation.append(np.zeros_like(fam_validation.loc[:, "status"]))

            origin_result = 0.5

            for i in np.arange(args.lower, args.threshold_t, args.interval):
                with Pool(processes=int((args.end_site - args.begin_site) / args.site_step)) as pool:
                    pool_list = []
                    pool_list_validation = []
                    j = 0

                    for k in np.arange(args.begin_site, args.end_site, args.site_step):
                        pool_list.append(("thread " + str(k), k, assoc_combine, i, bim, fam, args, origin_PRS[j], "test"))
                        pool_list_validation.append(("thread " + str(k), k, assoc_combine, i, bim_validation, fam_validation, args, origin_PRS_validation[j], "validation"))
                        j += 1

                    #pool.starmap(utils.workers, pool_list)
                    results = pool.starmap(utils.workers, pool_list)
                    results_validation = pool.starmap(utils.workers, pool_list_validation)
                    j = 0

                    for k in np.arange(args.begin_site, args.end_site, args.site_step):
                        origin_PRS[j] = results[j]
                        origin_PRS_validation[j] = results_validation[j]
                        j+=1

                j = 0
                for k in np.arange(args.begin_site, args.end_site, args.site_step):
                    a_origin, b_origin = utils.preprocess_for_roc(fam.loc[:, "status"], origin_PRS[j])
                    a_origin_validation, b_origin_validation = utils.preprocess_for_roc(fam_validation.loc[:, "status"], origin_PRS_validation[j])
                    #print(utils.bce(a_origin, b_origin), k, i, a_origin, b_origin)
                    #print(utils.bce(a_origin_validation, b_origin_validation), k, i, a_origin_validation, b_origin_validation, "validation")
                    if utils.bce(a_origin_validation, b_origin_validation) > origin_result:
                        print("result:", utils.bce(a_origin, b_origin))
                        print("result:", utils.bce(a_origin_validation, b_origin_validation))
                        origin_result = utils.bce(a_origin_validation, b_origin_validation)

                    else:
                        print("result:", utils.bce(a_origin, b_origin), "non")
                        print("result:", utils.bce(a_origin_validation, b_origin_validation), "non")

                    j+=1

    return origin_PRS
