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
import model

num=0
logistic="False"
clumped="False"
calculate_mod="False"
filed="False"
testmod="False"
time_point=time.time()
prsice="False"
group="False"
test_group="False"
conv="False"
validationfile="validationfile"
outpercent=1
validationpercent=1
threshold_t=0.05
begin_site=10
end_site=200
step=10
minsite=0.00000005
takestep=0.00005
symbol=0
testfile="testfile"
merge="False"
summary="False"
validationmod="False"
with_origin="False"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-file", help="The name of input file, without suffix.", type=str, required=True)
parser.add_argument("-v", "--validation-file", help="The name of validation file, if did not use this label, the system will use the last ten percent of the input dataset as the validation set.", type=str, default = "False")
parser.add_argument("-B", "--begin-site", help="The lower threshold of the split position of the system. The default value is 0. The value is no less than 0. Notice that the number can set a little higher when useing the t-split mode.", type=int, default=10)
parser.add_argument("-E", "--end-site", help="The higher threshold of the split position of the system. The default value is 200. Notice that the number can set a little lower when useing the t-split mode.", type=int, default=200)
parser.add_argument("-s", "--site-step", help="The threshold's step, too small value may cause the code's running step slow down. The value should no lower than (end_site - begin_site). The dafault value is 10.", type=int, default=10)
parser.add_argument("-b", "--beta", help="logistic association file.", action='store_true')
parser.add_argument("-k", "--keep-t-value", help="Use additive t value to rank the SNPs.", action='store_true')
parser.add_argument("-c", "--calculate-mod", help="The calculate-mods. The chosenable label include: Origin Mean-value T-compare Beta-compare", type=str, required=True, choices=['origin', 'mean', 't', 'beta'])
parser.add_argument("-t", "--test-file", help="The name of test file, if did not use this label, the system will use the last ten percent of the input dataset as the test set. If both of this label and the \"-v\" label have not used, the system will use the last ten percent of the input dataset as the test set, and next last ten percent of the input dataset as the validation set.", type=str, default = "False")
parser.add_argument("-H", "--threshold-t", help="The highest p value of the snps included into the system. The default is 0.05.", type=float, default=0.05)
parser.add_argument("-g", "--group-input", help="Input list mode.", action='store_true')
parser.add_argument("-G", "--group-test", help="Test list mode, the system cannot handle the test file with differ fam file, please try to merge the bfile by plink first, sorry.", action='store_true')
parser.add_argument("-V", "--group-validation", help="Test list mode, the system cannot handle the validation file with differ fam file, please try to merge the bfile by plink first, sorry.", action='store_true')
parser.add_argument("-I", "--interval", help="The step size of the threshold. The default is 5e-05.", type=float, default=0.00005)
parser.add_argument("-L", "--lower", help="The starting p-value threshold. The default is 5e-08.", type=float, default=0.00000005)
parser.add_argument("-T", "--thread", help="The number of thread, recommended no larger than the number of (end_site - begin_site) / site_step.", type=int, default=26)
parser.add_argument("-m", "--merge", help="Merge the result of additive and non-additive calculate model.", action='store_true')
parser.add_argument("-S", "--summary", help="Calculate from summary data.", action='store_true')
parser.add_argument("-w", "--with_origin", help="Calculate origin model in the same time calculating the other model.", action='store_true')
parser.add_argument("-K", "--keep-list", help="The list of people keep in all the process. The file need two column, the first one is family ID and the second is individual ID.", type=str)
parser.add_argument("-C", "--clump", help="Have clumped file.", action = "store_true")
args = parser.parse_args()
print(args)
if "." in args.input_file:
    print("The input file cannot include the character \".\" .")
    sys.exit(1)

if "." in args.validation_file:
    print("The validation file cannot include the character \".\" .")
    sys.exit(1)

if args.begin_site < 0:
    print("The begin site cannot smaller than 0.")
    sys.exit(1)

if args.end_site <= begin_site:
    print("The end site cannot smaller than or equal to begin site.")
    sys.exit(1)

if "." in args.test_file:
    print("The test file cannot include the character \".\" .")
    sys.exit(1)

if args.threshold_t > 1 or args.threshold_t < 0:
    print("The threshold of t can only between 0 and 1.")
    sys.exit(1)

if args.group_input:
    if utils.same_fam(args.input_file):
        print("Please make sure the list of input file have the same fam file, the system cannot handle the input file with differ fam file, please try to merge the bfile by plink first, sorry.")
        sys.exit(1)

if args.group_test:
    if utils.same_fam(args.test_file):
        print("Please make sure the list of test file have the same fam file, the system cannot handle the test file with differ fam file, please try to merge the bfile by plink first, sorry.")
        sys.exit(1)

if args.group_validation:
    if utils.same_fam(args.validation_file):
        print("Please make sure the list of validation file have the same fam file, the system cannot handle the validation file with differ fam file, please try to merge the bfile by plink first, sorry.")
        sys.exit(1)

if args.interval > 1 or args.interval < 0:
    print("The threshold of t can only between 0 and 1.")
    sys.exit(1)

if args.lower > 1 or args.lower < 0:
    print("The threshold of t can only between 0 and 1.")
    sys.exit(1)

if "." in args.keep_list:
    print("The keep list cannot include the character \".\" .")
    sys.exit(1)

model.main(args)
