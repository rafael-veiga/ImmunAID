# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 02:17:09 2024

@author: rafae
"""

import pyreadr as py
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc,roc_auc_score
from joblib import Parallel, delayed
n_jobs = 4
seed = 12345



class AUC_CV:    
    def __init__(self,fpr_tam=100):
        self.fpr_mean = np.linspace(0, 1, fpr_tam)
        self.size = 0
        self.fprs=[]
        self.tprs=[]

    def add_points(self,fpr,tpr):
        self.size = self.size+1
        self.fprs.append(fpr)
        self.tprs.append(tpr)
        
    def get_cur(self,cv):
        return (self.fprs[cv],self.tprs[cv])
    
    def get_mean(self):
        tpr_aux = np.ones((self.size,len(self.fpr_mean)))
        for cv in range(self.size):
            aux = np.interp(self.fpr_mean, self.fprs[cv], self.tprs[cv])
            aux[0] = 0.0
            tpr_aux[cv,:] = aux
        tpr_mean = np.ones(len(self.fpr_mean))
        tpr_mean[-1] = 1.0
        tpr_l = np.ones(len(self.fpr_mean))
        tpr_h = np.ones(len(self.fpr_mean))
        for i in range(len(self.fpr_mean)):
            tpr_mean[i] = np.mean(tpr_aux[:,i])
            se = (1.96*np.std(tpr_aux[:,i]))/np.sqrt(self.size)
            tpr_l[i] = tpr_mean[i]-se
            tpr_h[i] = tpr_mean[i]+se
        data = pd.DataFrame({"fpr":self.fpr_mean,"tpr":tpr_mean,"tpr_l":tpr_l,"tpr_h":tpr_h})
        data.loc[data.tpr_l<0,"tpr_l"] = 0
        data.loc[data.tpr_h>1,"tpr_h"] = 1
        data["auc"] = auc(self.fpr_mean,tpr_mean)
        data["auc_l"] = auc(self.fpr_mean,tpr_l)
        data["auc_h"] = auc(self.fpr_mean,tpr_h)
        return data

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



def lr_auc():
    df = py.read_r("./pos_data/data_nor.rds")[None]
    df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
    catego = list(df["disease"].cat.categories)
    catego.remove("Healthy")
    del df["id"]
    del df["batch"]
    skf = StratifiedKFold(n_splits=10,random_state=seed,shuffle=True)
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
        auc_cv = AUC_CV()
        for cv, (train_index, test_index) in enumerate(skf.split(x, y)):
            LR.fit(x[train_index], y[train_index])
            yhat = LR.predict(x[test_index])
            fpr, tpr, thresholds = roc_curve(y[test_index],yhat,pos_label=1)
            auc_cv.add_points(fpr, tpr)
        data_aux = auc_cv.get_mean()
        data_aux["disease"] = d
        data = pd.concat([data, data_aux],ignore_index=True)
    data.to_csv("./result/table/auc_curv/auc_curv.csv",index=False)

    
def rf_res():
    df = py.read_r("./pos_data/data_nor.rds")[None]
    df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
    catego = list(df["disease"].cat.categories)
    #catego.remove("Healthy")
    #df = df.loc[df.disease!="Healthy"]
    del df["id"]
    del df["batch"]
    skf = StratifiedKFold(n_splits=10,shuffle=True,random_state=seed)
    rf1 = RandomForestClassifier(max_depth=3,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed)
    rf2 = RandomForestClassifier(max_depth=4,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed+1)
    rf3 = RandomForestClassifier(max_depth=5,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed+2)
    rf4 = RandomForestClassifier(max_depth=8,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed+3)
    sk = StandardScaler()
    data = pd.DataFrame()
    for d1 in range(len(catego)):
        for d2 in range(len(catego)):
            if d1 != d2:
                aux = df.copy()
                aux["out"] = np.nan
                aux.loc[aux["disease"]==catego[d2],"out"] = 0
                aux.loc[aux["disease"]==catego[d1],"out"] = 1
                del aux["disease"]
                aux.dropna(how='any',inplace=True)
                y = np.array(aux["out"])
                x = aux
                del x["out"]
                x = sk.fit_transform(x)
                for cv, (train_index, test_index) in enumerate(skf.split(x, y)):
                    print("d1 "+str(d1)+" d2 "+str(d2)+" cv"+str(cv))
                    rf1.fit(x[train_index], y[train_index])
                    rf2.fit(x[train_index], y[train_index])
                    rf3.fit(x[train_index], y[train_index])
                    rf4.fit(x[train_index], y[train_index])
                    par = []
                    yhat = 0
                    if rf4.oob_score_ >= max(rf2.oob_score_,rf3.oob_score_,rf1.oob_score_):
                        yhat = rf4.predict_proba(x[test_index])[:,1]
                        par.append(8)
                    if rf3.oob_score_ >= max(rf2.oob_score_,rf1.oob_score_,rf4.oob_score_):
                        yhat = rf3.predict_proba(x[test_index])[:,1]
                        par.append(5)
                    if rf2.oob_score_ >= max(rf1.oob_score_,rf3.oob_score_,rf4.oob_score_):
                        yhat = rf2.predict_proba(x[test_index])[:,1]
                        par.append(4)
                    if rf1.oob_score_ >= max(rf2.oob_score_,rf3.oob_score_,rf4.oob_score_):
                        yhat = rf1.predict_proba(x[test_index])[:,1]
                        par.append(3)
                    data_aux = pd.DataFrame({"y_real":y[test_index],"y_hat":yhat})
                    data_aux["size_tree"] = par[0]
                    data_aux["cv"] = cv
                    data_aux["d1"] = catego[d1]
                    data_aux["d2"] = catego[d2]
                    data = pd.concat([data,data_aux ], ignore_index=True)
    data.to_csv("./result/RF_auc1.csv",index=False)


    
    
            
    
    

def rf_auc():
    df = pd.read_csv("./result/RF_auc1.csv")
    dis = pd.unique(list(df["d1"])+list(df["d2"]))
    data = pd.DataFrame()
    for d1 in dis:
        for d2 in dis:
            if d1!=d2:
                par = []
                auc_cv = AUC_CV()
                for cv in range(10):
                    print("d1 "+str(d1)+" d2 "+str(d2)+" cv"+str(cv))
                    aux = df.loc[((df["cv"]==cv)&(df["d1"]==d1))&(df["d2"]==d2)]
                    fpr, tpr, thresholds = roc_curve(aux["y_real"],aux["y_hat"],pos_label=1)
                    auc_cv.add_points(fpr,tpr)
                    par.append(aux["size_tree"].iloc[0])
                par = pd.Series(par)
                tb = par.value_counts().sort_values(ascending=False)
                data_aux = auc_cv.get_mean()
                data_aux["d1"] = d1
                data_aux["d2"] = d2
                data_aux["par"]= tb.index[0]
                data = pd.concat([data, data_aux],ignore_index=True)
    data.to_csv("./result/RF_auc2.csv",index=False)

def rf_imp1():
    res = pd.read_csv("./result/RF_auc2.csv")
    df = py.read_r("./pos_data/data_nor.rds")[None]
    df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
    catego = list(df["disease"].cat.categories)
    #catego.remove("Healthy")
    #df = df.loc[df.disease!="Healthy"]
    del df["id"]
    del df["batch"]
    sk = StandardScaler()
    data = pd.DataFrame()
    for d1 in range(len(catego)):
        for d2 in range(len(catego)):
            if d1 != d2:
                print("d1 "+str(d1)+" d2 "+str(d2))
                aux = df.copy()
                aux["out"] = np.nan
                aux.loc[aux["disease"]==catego[d2],"out"] = 0
                aux.loc[aux["disease"]==catego[d1],"out"] = 1
                del aux["disease"]
                aux.dropna(how='any',inplace=True)
                y = np.array(aux["out"])
                x = aux
                del x["out"]
                col = x.columns
                x = sk.fit_transform(x)
                par = res.loc[(res["d1"]==catego[d1]) & (res["d2"]==catego[d2]),"par"].iloc[0]
                rf = RandomForestClassifier(max_depth=par,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed)
                rf.fit(x,y)
                aux = pd.DataFrame({"vari":col,"imp":rf.feature_importances_})
                aux = aux.loc[aux["vari"]!="sex"]
                aux = aux.loc[aux["vari"]!="age"]
                data_aux = pd.DataFrame({"vari":list(aux["vari"]),"imp":list(aux["imp"])})
                data_aux["d1"] = catego[d1]
                data_aux["d2"] = catego[d2]
                data = pd.concat([data, data_aux],ignore_index=True)
    
    data.to_csv("./result/RF_imp.csv",index=False)

lr_auc()
#rf_res()
#rf_auc()
#rf_imp1()


# df = py.read_r("./pos_data/data_nor.rds")[None]
# df["sex"]=pd.get_dummies(df["sex"]).iloc[:,0]
# catego = list(df["disease"].cat.categories)
# catego.remove("Healthy")
# df = df.loc[df.disease!="Healthy"]
# del df["id"]
# del df["batch"]          
# y = df["disease"]
# x = df
# del x["disease"]
# col = x.columns
# rf = RandomForestClassifier(max_depth=8,n_jobs=n_jobs,n_estimators=5001,oob_score=True,random_state=seed)
# rf.fit(x,y)
# data = pd.DataFrame({"vari":col,"imp":rf.feature_importances_})
# data.loc[data["vari"]!="sex"]
# data.loc[data["vari"]!="age"]
# data.sort_values("imp",ascending=False,inplace=True)
# data.to_csv("./result/RF_imp2.csv",index=False)




















