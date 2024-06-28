source("functions_pre.R")
df_raw = readRDS(file = "./pos_data/data.rds")
col = colnames(df_raw)
aux = col[get_fil(col,type="pop_in_pop")]
aux = c("id","disease","age","sex",aux)
df_raw = df_raw[,aux]
write.csv(df_raw,"./pos_data/ImmunAID.csv",row.names = FALSE,quote = TRUE)

### box-cox
fil1 = get_fil(colnames(df_raw),type = c("pop_in_pop"))
res = box_cox_trans(df_raw,fil=fil1)
df_raw = res$df
trans = res$trans

### imputation
df_raw = imputation(df_raw,fil=fil1,k=3)

#save
saveRDS(df_raw, file = "./pos_data/data_nor.rds")


## inverse box-cox
df_raw = inv_box_cox_trans(df_raw,trans,fil1)

#save
saveRDS(df_raw, file = "./pos_data/data_imp.rds")


