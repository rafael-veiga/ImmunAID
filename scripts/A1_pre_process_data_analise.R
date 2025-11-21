library(tidyverse)
library(VIM)
library(forecast)
data = readRDS(file = "./data_raw.rds")

remove_missing <- function(df, per,colunas){
  total = nrow(df)
  total = round(total*per)
  filt = rep(TRUE,ncol(df))
  for(i in 1:ncol(df)){
    if(colnames(df)[i] %in% colunas){
      col = sum(is.na(df[[i]]))
      if(col>=total){
        filt[i] = FALSE
      }
    }
  }
  return(df[,filt])
}

imputation <- function(df,nk,colunas){
  df_aux = df[,colunas]
  df_base = df
  df_aux = kNN(df_aux,k = nk,imp_var = FALSE)
  df_base[,colunas] = df_aux
  return(df_base)
}

get_fil <- function(col,type="all",pop="all",mark="all"){
  fil_t = rep(FALSE,length(col))
  for(t in type){
    if(t=="all"){
      fil_t=rep(TRUE,length(col))
      break
    }
    if(t=="other"){
      fil_t = (fil_t) | (col=="sample") | (col=="batch") | (col=="id") | (col=="type") | (col=="age") | (col=="sex") | (col=="disease") | (col=="treated") | (col=="fil_dis")
    }
    if(t=="pop_in_pop"){
      fil_t = (fil_t) | (grepl(" in ",col) & (!grepl(" median ",col)) )
    }
    if(t=="mark_in_pop"){
      fil_t = (fil_t) | grepl(" median in ",col)
    }
  }
  fil_p = rep(FALSE,length(col))
  for(p in pop){
    if(p=="all"){
      fil_p=rep(TRUE,length(col))
      break
    }
    f = as.logical(sapply(col,function(s){
      if(grepl(" in ",s)){
        return (unlist(str_split(s," in "))[2])
      }
      return(s)
    })==p)
    ini = as.logical(sapply(col,function(s){
      if(grepl(" in ",s) & (!grepl(" median ",s))){
        return (unlist(str_split(s," in "))[1])
      }
      return(s)
    })==p)
    f_ini = f | ini
    fil_p = fil_p | f_ini
  }
  fil_m = rep(FALSE,length(col))
  for(m in mark){
    if(m=="all"){
      fil_m=rep(TRUE,length(col))
      break
    }
    f = as.logical(sapply(col,function(s){return (unlist(str_split(s," median in "))[1])})==m)
    fil_m = fil_m | f
  }
  return(fil_m & fil_p & fil_t)
}

box_cox_trans <- function(df,fil){
  df1 = df[,fil]
  df_min = rep(NA,length(colnames(df1)))
  lampda = rep(NA,length(colnames(df1)))
  std_mean =  rep(NA,length(colnames(df1)))
  std_sd =  rep(NA,length(colnames(df1)))
  for(i in 1:length(colnames(df1))){
    df_min[i] = min(df1[,i],na.rm = TRUE)-1
    df1[,i] = df1[,i] - df_min[i]
    lampda[i] = BoxCox.lambda(df1[,i][!is.na(df1[,i])], method = c("loglik"), lower = -5, upper = 5)
    if(lampda[i]==0){
      df1[,i] = log(df1[,i])
    }else{
      df1[,i] = (sign(df1[,i])*((abs(df1[,i]))^lampda[i])-1)/lampda[i]
    }
    #standat
    std_mean[i] = mean(df1[,i],na.rm = TRUE)
    std_sd[i] = sd(df1[,i],na.rm = TRUE)
    df1[,i] = (df1[,i]-std_mean[i])/std_sd[i]
  }
  df[,fil] = df1
  trans = data.frame(lampda=lampda,std_mean=std_mean,std_sd = std_sd,df_min = df_min)
  return(list(df=df,trans=trans))
}

inv_box_cox_trans <- function(df,trans,fil){
  df1 = df[,fil]
  lambda = trans$lampda
  std_mean =  trans$std_mean
  std_sd =  trans$std_sd
  df_min = trans$df_min
  for(i in 1:length(colnames(df1))){
    #inv standart
    df1[,i] = df1[,i] * std_sd[i]+std_mean[i]
    #inv boxcox
    if (lambda[i] == 0) {
      df1[,i] = exp(df1[,i])
    }else{
      df1[,i] =  (df1[,i] * lambda[i] + 1)^(1/lambda[i])
    }
    #add min
    df1[,i] = df1[,i] + df_min[i]
  }
  df[,fil] = df1
  return(df)
}
##################################################################################################
df_p = data[["df_p"]]
df_f = data[["df_f"]]

df_p = df_p %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
    ))
df_f = df_f %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
  ))

### remove MFI

df_f = df_f[,get_fil(colnames(df_f),type = c("other","pop_in_pop"))]
col_p = colnames(df_p)[6:ncol(df_p)]
col_f = colnames(df_f)[6:ncol(df_f)]

### remove missing > 20%
df_p = remove_missing(df_p,per=0.2,col_p)
df_f = remove_missing(df_f,per=0.2,col_f)

### save
data[["df_p"]] = df_p
data[["df_f"]] = df_f
saveRDS(data, file = "./data.rds")

### box-cox
res_f = box_cox_trans(df_f,colnames(df_f) %in% col_f)
df_f = res_f$df

### imputation
df_f = imputation(df_f,3,colnames(df_f) %in% col_f)
df_p = imputation(df_p,3,colnames(df_p) %in% col_p)

#save
data[["df_f"]] = df_f
data[["df_p"]] = df_p
saveRDS(data, file = "./data_nor.rds")


## inverse box-cox
df_f = inv_box_cox_trans(df_f,res_f$trans,colnames(df_f) %in% col_f)

#save
data[["df_f"]] = df_f
saveRDS(data, file = "./data_imp.rds")
