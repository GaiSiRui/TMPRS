# TMPRS

TMPRS is a software package for calculating, applying, and evaluating polygenic risk score (PRS) results. The most distinctive feature of TMPRS is its ability to incorporate genetic information from non-additive inheritance models and produce results that take both additive and non-additive genetics into account, along with corresponding analytical visualizations. It can also prune SNPs based on linkage disequilibrium and p-value (i.e., "clumping") and include desired covariates.

TMPRS is a software package written in Python or shell script and runs on top of PLINK, one of the most widely used bioinformatics tools. TMPRS operates as a command-line program with a variety of user options and is freely available for download, compatible with Unix/Linux systems.

# Prerequisite

Plink(>=1.90), matplotlib(>=3.6.2), pandas(>=2.0.3), scikit-learn(>=1.0.2), statsmodels(>=0.14.1), numpy(>=1.24.4), pandas(>=2.0.3), pyplink(>=1.3.7), scipy(>=1.9.3)

# Example

## python

python main.py -i (training set) -v (validation set) -c t -t (test set) -b -o (covariant file) -K (keep list)

## shell

sbatch file_get_feature_01_main.sh -i (training set) -v (validation set) -c False -t (test set) --with_origin -K (keep list) --beta
