library(future.apply)
library(ggpubr)
plan(multisession)

source("functions_pre.R")


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
df$ImmunAID_identifier = as.character(sprintf("%05i", df$ImmunAID_identifier))
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
dis[dis=="Negative_control"] = "Healthy"
dis[dis=="Behcet_(Vasculitis)"] = "Behcet"
dis[dis=="Inflammation_of_unknown_origin"] = "Inflammation of unknown origin"
levels(df_p$disease) = dis


df_p$fil_dis = FALSE
df_p$fil_dis[df_p$disease %in% c("AOSD","Inflammation of unknown origin","FMF","Healthy","Behcet")] = TRUE
############################################################################################################################
# save
data_raw = df_p
saveRDS(data_raw, file = "./pos_data/data_raw.rds")
saveRDS(df_s, file = "./pos_data/data_s.rds")


