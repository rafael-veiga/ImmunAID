# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:12:11 2024

@author: rafae
"""

import pyreadr as py
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc,roc_auc_score
from joblib import Parallel, delayed
n_jobs = 15
seed = 12345

def compute_auc(x, y, LR, skf):
    aucs = np.empty(skf.n_splits)
    for cv, (train_index, test_index) in enumerate(skf.split(x, y)):
        LR.fit(x[train_index], y[train_index])
        yhat = LR.predict(x[test_index])
        aucs[cv] = roc_auc_score(y[test_index], yhat)
    return np.mean(aucs)

dis = {"Inflammation of unknown origin":40,"AOSD":82,"FMF":55,"Behcet":61}

def calculate_auc_vs_n():
    df = py.read_r("./pos_data/data_nor.rds")[None]
    df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
    catego = list(df["disease"].cat.categories)
    catego.remove("Healthy")
    del df["id"]
    del df["batch"]
    dis = {"Inflammation of unknown origin": 15,"AOSD":26,"FMF":7,"Behcet":15}
    skf = StratifiedKFold(n_splits=5,random_state=seed,shuffle=True)
    LR = LogisticRegression(C=0.1,penalty="l2",solver="liblinear",max_iter=1000)
    data = pd.DataFrame()
    for d in catego:
        res = pd.read_csv("./result/table/marks/"+ d +".csv")["var"]
        auc_list = np.empty(len(res)-1)
        auc_list[:] = np.nan
        df1 = df.copy()
        df1["y"] = np.nan
        df1.loc[df1["disease"]=="Healthy","y"] = 0
        df1.loc[df1["disease"]==d,"y"] = 1
        df1.dropna(subset=['y'], inplace=True)
        y = df1["y"]
        y = np.array(y, order="C")
        del df1["y"]
        del df1["disease"]
        for n in range(1, len(res)):
            print("n " + str(n) + " disease " + d)
            if d =="FMF":
                vari = ["age"] + list(res[range(n)])
            else:
                vari = ["sex", "age"] + list(res[range(n)])
            x = df1[vari].copy()
            x = np.array(x, order="C")
            
            # Paralelizar esta parte
            aucs = Parallel(n_jobs=n_jobs)(delayed(compute_auc)(x, y, LR, skf) for i in range(1000))
            
            aucs = np.array(aucs)
            auc_list[n-1] = np.mean(aucs)
                        
        data_aux = pd.DataFrame({"num":range(1,len(res)),"auc":auc_list})
        data_aux["disease"] = d
        data = pd.concat([data,data_aux ], ignore_index=True)
    
    data.to_csv("./result/table/auc_n.csv",index=False)

#calculate_auc_vs_n()



df = py.read_r("./pos_data/data_nor.rds")[None]
df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
catego = list(df["disease"].cat.categories)
catego.remove("Healthy")
del df["id"]
del df["batch"]
skf = StratifiedKFold(n_splits=5,random_state=seed,shuffle=True)
LR = LogisticRegression(C=0.1,penalty="l2",solver="liblinear",max_iter=1000)
data = pd.DataFrame()
for d in catego:
    aux=df.copy()
    res = pd.read_csv("./result/table/marks/"+ d +".csv")
    res = res["var"][0:dis[d]]
    if d=="FMF":
        res = ["age","disease"]+list(res)
    else:
        res = ["sex","age","disease"]+list(res)
    aux = aux[res]
    aux["y"] = np.nan
    aux.loc[aux["disease"]=="Healthy","y"] = 0
    aux.loc[aux["disease"]==d,"y"] = 1
    aux.dropna(subset=['y'], inplace = True)
    y = aux["y"]
    x = aux
    del x["y"]
    del x["disease"]
    y = np.array(y,order="C")
    x = np.array(x,order="C")
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    for i in range(200):
        print(i)
        for cv, (train_index, test_index) in enumerate(skf.split(x, y)):
            LR.fit(x[train_index], y[train_index])
            yhat = LR.predict(x[test_index])
            fpr, tpr, thresholds = roc_curve(y[test_index],yhat,pos_label=1)
            
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    se = 1.96*np.std(tprs, axis=0)/np.sqrt(10)
    tpr_l = mean_tpr - se
    tpr_h = mean_tpr + se
    tpr_l = [0 if v<0 else v for v in tpr_l]
    tpr_h = [1 if v>1 else v for v in tpr_h]
    auc_l = auc(mean_fpr, tpr_l)
    auc_h = auc(mean_fpr, tpr_h)
    data_aux = pd.DataFrame({"tpr":mean_tpr,"fpr":mean_fpr,"tpr_l":tpr_l,"tpr_h":tpr_h})
    data_aux["auc"] = auc(mean_fpr, mean_tpr)
    data_aux["auc_l"] = auc(mean_fpr, tpr_l)
    data_aux["auc_h"] = auc(mean_fpr, tpr_h)
    data_aux["disease"] = d
    data = pd.concat([data, data_aux],ignore_index=True)
data.to_csv("./result/table/auc_curv/auc_curv.csv",index=False)
    
            
        
    
    
    