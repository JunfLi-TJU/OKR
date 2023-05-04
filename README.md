# OKR
Source codes of algorithms and datasets for our paper "Nearly Optimal Algorithms with Sublinear Computational Complexity for Online Kernel Regression", 
accepted in ICML 2023.

We implement all algorithms with R on a Windows machine with 2.8 GHz Core(TM) i7-1165G7 CPU, 
execute each experiment 10 times with random permutation of all datasets and average all of the results.

The default path of codes is "D:/experiment/Conference Paper/ICML/ICML2023/code/". 

The path of datasets is "D:/experiment/online learning dataset/regression/". 

The store path is "D:/experiment/Conference Paper/ICML/ICML2023/result/". 

You can also change all of the default paths.

The baseline algorithms include: FOGD, NOGD, PROS-N-KONS and CON-KONS. 
Our algorithms include: AOGD-ALD and NONS-ALD.

The datasets are downloaded from: https://archive.ics.uci.edu/ml/index.php.

regression datasets:
parkinson (Num:5875, Fea: 16), elevators (Num:16599, Fea: 18), cpusmall (Num:8192, Fea: 12), bank (Num:8192, Fea: 32),
ailerons (Num:13750, Fea: 40), calhousing (Num:14000, Fea: 8), Year (Num:51630, Fea: 90), TomsHardware (Num:28179, Fea: 96)
