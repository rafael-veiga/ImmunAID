library(tidyverse)

df = read_csv("./data/Patient_overview.csv")
df$Disease = as.factor(df$Disease)
colnames(df)[1] = "sample"
colnames(df)[2] = "disease"
colnames(df)[3] = "age"
colnames(df)[4] = "sex"
dis = levels(df$disease) 
dis[dis=="AOSD"] = "Still's disease"
dis[dis=="Negative control"] = "Healthy"
dis[dis=="Behcet (Vasculitis)"] = "Behcet"
dis[dis=="Inflammation of unknown origin"] = "Autoinflammation of unknown origin"
levels(df$disease) = dis


df$fil_dis = FALSE
df$fil_dis[df$disease %in% c("Still's disease","Autoinflammation of unknown origin","FMF","Healthy","Behcet")] = TRUE
df = df[df$fil_dis,]
df = df[df$fil_dis,]
df$fil_dis = NULL

############################################################################################################################
# remove sample 01002 and 27004
df = df[df$sample!="01002",]
df = df[df$sample!="27004",]
############################################################################################################################
aux = read.table("./data/Patient_stratification.csv",sep=",")
colnames(aux) = aux[1,]
aux = aux[2:nrow(aux),]
df$treated = NA

treated = aux$Treated[aux$Treated!=""]
naive = aux$Naive[aux$Naive!=""]
df$treated[df$sample %in% naive] = 0
df$treated[df$sample %in% treated] = 1
############################################################################################################################

df_p = read_csv("./data/immunAIDproteomics.csv")
protein_lab = df_p[,c(1,2,3)]
df_p = df_p[,c(c(1),4:ncol(df_p))]
col = df_p[,1][[1]]

df_p = df_p[,2:ncol(df_p)]
sample = colnames(df_p)
df_p = t(df_p)
df_p = as_tibble(df_p)
colnames(df_p) = col
df_p$sample = sample
df_p = merge.data.frame(df,df_p,by = "sample",all.x = FALSE, all.y = FALSE)
data = list(df_p=df_p,lab=protein_lab)
meta = df
############################################################################################################################
library(future.apply)
library(ggpubr)
#plan(multisession)
open_file_p <- function(file){
  df = read_csv(paste0("./data/AutoGate_results_patients/",file))
  col = colnames(df)[2:length(colnames(df))]
  df <- df %>% mutate_at(col, as.double)
  df$batch =  unlist(strsplit(file,"_statistics.csv"))[1]
  return(df)
}

open_file_s <- function(file){
  df = read_csv(paste0("./data/AutoGate_results_IC/",file))
  col = colnames(df)[2:length(colnames(df))]
  df <- df %>% mutate_at(col, as.double)
  df$batch =  unlist(strsplit(file,"_IC_statistics.csv"))[1]
  return(df)
}

get_marks <- function(col){
  marks = list()
  pop = list()
  fil_other = (!grepl(" in ",col))
  fil_pop = (grepl(" in ",col) & (!grepl(" median ",col)))
  fil_mark = grepl(" median in ",col)
  other = col[fil_other]
  aux = col[fil_mark]
  aux = str_split(aux," median in ")
  for(i in aux){
    marks = append(marks,i[1])
    pop = append(pop,i[2])
  }
  aux = col[fil_pop]
  aux = str_split(aux," in ")
  for(i in aux){
    pop = append(pop,i[1])
    pop = append(pop,i[2])
  }
  
  type = c("other","pop_in_pop","mark_in_pop")
  res = list()
  res[["mark"]]= unique(unlist(marks))
  res[["pop"]] = unique(unlist(pop))
  res[["other"]] = other
  res[["type"]] = type
  
  return(res)
}

get_fil <- function(col,type="all",pop="all",mark="all"){
  fil_t = rep(FALSE,length(col))
  for(t in type){
    if(t=="all"){
      fil_t=rep(TRUE,length(col))
      break
    }
    if(t=="other"){
      fil_t = (fil_t) | (col=="Sample") | (col=="batch") | (col=="id") | (col=="type") | (col=="age") | (col=="sex") | (col=="disease") | (col=="treated") | (col=="fil_dis")
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
#source("functions_pre.R")
# read files
files = list.files("./data/AutoGate_results_patients")

df = future_sapply(files,open_file_p)

df_p = bind_rows(df, .id = "column_label")

files = list.files("./data/AutoGate_results_IC")

df = future_sapply(files,open_file_s)

df_s = bind_rows(df, .id = "column_label")
df_p <- df_p %>% select(!contains(" iqr "))
df_s <- df_s %>% select(!contains(" iqr "))
df_p$column_label = NULL
df_s$column_label = NULL
df_s = df_s[,colnames(df_p)]
df_s$id = NA
df_s$type = NA

df_p$id = NA
df_p$type = NA
df_p$id = as.character(sapply(df_p$Sample, function(s){substr(s,1,5)}))
df_p$type = as.factor(as.character(sapply(df_p$Sample, function(s){substr(s,6,6)})))
rm(files)
df_p$Sample = NULL
df_p$disease = NA
df_p$age = NA
df_p$sex = NA
df_s$disease = NA
df_s$age = NA
df_s$sex = NA

df = read_csv("./data/Patient_overview.csv")
df$ImmunAID_identifier = as.character(sprintf("%05i",as.numeric( df$ImmunAID_identifier)))
id = unique(df_p$id)
for(i in id){
  fil = df_p$id==i
  fil2 = df$ImmunAID_identifier==i
  df_p$disease[fil] = df$Disease[fil2]
  df_p$age[fil] = df$Age[fil2]
  df_p$sex[fil] = df$Sex[fil2]
}

df_p$sex = as.factor(df_p$sex)
df_s$id = df_s$Sample
df_s$Sample = NULL

res = get_marks(colnames(df_p))

# get only A group
df_p = df_p[df_p$type=="A",]
df_p$type = NULL

# get remove CD80 all
fil = get_fil(colnames(df_p),type = "mark_in_pop",mark = "CD80")
df_p = df_p[,!fil]
df_p <- df_p %>% select(-contains("CD80"))

# edit disease
df_p$disease = as.factor(df_p$disease)
dis = levels(df_p$disease) 
dis[dis=="Negative control"] = "Healthy"
dis[dis=="Behcet (Vasculitis)"] = "Behcet"
dis[dis=="Inflammation of unknown origin"] = "Inflammation of unknown origin"
levels(df_p$disease) = dis


df_p$fil_dis = FALSE
df_p$fil_dis[df_p$disease %in% c("AOSD","Inflammation of unknown origin","FMF","Healthy","Behcet")] = TRUE
df_p = df_p[df_p$fil_dis,]
df_p$fil_dis = NULL
df_p$sample = df_p$id
df_p$id = NULL
df_p = df_p[df_p$sample!="01002",]
df_p = df_p[df_p$sample!="27004",]
df_p$disease= NULL
df_p$age = NULL
df_p$sex = NULL

df_p = merge.data.frame(meta,df_p,by="sample",all.x = FALSE, all.y = FALSE)
df_p$batch = NULL

data[["df_f"]] = df_p
############################################################################################################################
saveRDS(data,"./data_raw.rds")


