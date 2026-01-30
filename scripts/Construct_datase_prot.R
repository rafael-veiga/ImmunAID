library(tidyverse)
library(readxl)
df_f = read_excel("./data/Supplementary Spreadsheet 3.xlsx")
df_f$disease = as.factor(df_f$disease)
df_f$sex = as.factor(df_f$sex)
df_f$treated = as.factor(df_f$treated)
colnames(df_f)[1] = "sample"
df_f[, 6:ncol(df_f)] <- lapply(df_f[, 6:ncol(df_f)], as.numeric)

############################################################################################################################
df_p = read_excel("./data/Supplementary Spreadsheet 4.xlsx")
protein_lab = df_p[,c(1,2,3)]
df_p = df_p[,c(c(1),4:ncol(df_p))]
col = df_p[,1][[1]]
df_p = df_p[,2:ncol(df_p)]
sample = colnames(df_p)
df_p = t(df_p)
df_p = as_tibble(df_p)
colnames(df_p) = col
df_p$sample = sample
df_p = merge.data.frame(df_f[,1:5],df_p,by = "sample",all.x = FALSE, all.y = FALSE)
df_f = as.data.frame(df_f)
df_p = as.data.frame(df_p)
protein_lab = as.data.frame(protein_lab)
data = list(df_f = df_f,df_p=df_p,lab=protein_lab)
saveRDS(data,"./data_raw.rds")
