
source("functions_pre.R")
df_raw_t = readRDS(file = "./pos_data/data_raw.rds")
df_raw_t = df_raw_t %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
    ))



### select disease

df_raw = df_raw_t[df_raw_t$fil_dis,]
df_raw$fil_dis = NULL
df_raw$disease = as.character(df_raw$disease)
df_raw$disease = as.factor(df_raw$disease)

df_raw_t$fil_dis = NULL
df_raw_t$disease = as.character(df_raw_t$disease)
df_raw_t$disease = as.factor(df_raw_t$disease)

### remove missing > 20%
df_raw = remove_missing(df_raw,per=0.2)

### remove MFI

df_raw = df_raw[,get_fil(colnames(df_raw),type = c("other","pop_in_pop"))]

### save
df_raw_t = df_raw_t[,colnames(df_raw)]
saveRDS(df_raw, file = "./pos_data/data.rds")
saveRDS(df_raw_t, file = "./pos_data/data_t.rds")

### box-cox
fil1 = get_fil(colnames(df_raw),type = c("pop_in_pop"))
res = box_cox_trans(df_raw,fil=fil1)
df_raw = res$df
trans = res$trans
df_raw_t = box_cox_trans_add(df_raw_t,trans,fil1)

### imputation
df_raw = imputation(df_raw,fil=fil1,k=3)
df_raw_t = imputation(df_raw_t,fil=fil1,k=3)

#save

saveRDS(df_raw, file = "./pos_data/data_nor.rds")
saveRDS(df_raw_t, file = "./pos_data/data_nor_t.rds")

## inverse box-cox
df_raw = inv_box_cox_trans(df_raw,trans,fil1)
df_raw_t = inv_box_cox_trans(df_raw_t,trans,fil1)

#save
saveRDS(df_raw, file = "./pos_data/data_imp.rds")
saveRDS(df_raw_t, file = "./pos_data/data_imp_t.rds")


