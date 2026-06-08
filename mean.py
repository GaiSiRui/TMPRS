import pandas as pd
import statsmodels.api as sm
import numpy as np
import csv

pheno = pd.read_csv("alzheimer_train.pheno", sep="\s+")
covar = pd.read_csv("side_information_including_BMI.tab", sep="\s+")
geno_name = pd.read_csv("alzheimer_train.name", header = None , sep="\s+")
print(geno_name)

with open('mean.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['rs_name', 'aa_mean_adj', 'Aa_mean_adj', 'AA_mean_adj'])
    file.close()

for rs_number, rs_name in geno_name.iterrows():
    data = pd.merge(pheno, covar, on=["FID", "IID"])
    target_snp = rs_name
    columns_to_load = ["FID", "IID", target_snp.iloc[0]]

    try:
        geno = pd.read_csv("alzheimer_train.raw", sep="\s+", usecols=columns_to_load)

    except Exception as e:
        print("error", e)

    data = pd.merge(data, geno[["FID", "IID", rs_name.iloc[0]]], on=["FID", "IID"])

    X = data["sex"].values.reshape(-1, 1)
    X = sm.add_constant(X)
    y = data["PHENO"]
    model = sm.OLS(y, X).fit()
    data["RESIDUAL"] = model.resid

    with open('mean.csv', mode='a', newline='') as file:
        aa_mean_adj = data.loc[np.isclose(data[rs_name], 0.0), 'RESIDUAL'].mean() + y.mean()
        Aa_mean_adj = data.loc[np.isclose(data[rs_name], 1.0), 'RESIDUAL'].mean() + y.mean()
        AA_mean_adj = data.loc[np.isclose(data[rs_name], 2.0), 'RESIDUAL'].mean() + y.mean()
        writer = csv.writer(file)
        writer.writerow([rs_name, aa_mean_adj, Aa_mean_adj, AA_mean_adj])
