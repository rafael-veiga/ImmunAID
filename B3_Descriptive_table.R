source("functions_pre.R")


#table 1
df = readRDS(file = "./pos_data/data.rds")
tab = as.data.frame(table(df$disease))
outcome = tab$Var1
N = rep(NA,length(outcome))
rm(tab)
df = df[,get_fil(colnames(df),type = "other")]
for(i in 1:length(outcome)){
  N[i] = length(unique(df$id[df$disease==outcome[i]]))
}
tabela = data.frame(disease = outcome,N=N)
tabela = tabela[order(N,decreasing = TRUE),]

tabela$gender_me = NA
tabela$gender_p = NA
tabela$age_median = NA
tabela$age_iqr = NA
tabela$age_p = NA
tabela$min = NA
tabela$max = NA


df2 = df[,get_fil(colnames(df),type = "other")]
sample = unique(df2$id)
df = data.frame()
for(si in sample){
  aux = df2[df2$id==si,]
  if(nrow(aux)>1){
    a = aux
    if(a$type[1]==a$type[2]){
      print("erro")
      print(a)
    }
    if(a$sex[1]!=a$sex[2]){
      print("erro")
      print(a)
    }
    if(a$age[1]!=a$age[2]){
      print("erro")
      print(a)
    }
    if(a$disease[1]!=a$disease[2]){
      print("erro")
      print(a)
    }
    a = row.names(a)[2]
    df2 = df2[row.names(df2)!=a,]
  }
}

df = df2
#Healthy
out = "Healthy"
a = data.frame(table(df$sex[df$disease==out]))
tabela$gender_me[1] = paste0(a$Freq[a$Var1=="Female"],"(",sprintf("%0.1f",100*a$Freq[a$Var1=="Female"]/tabela$N[1]),"%)")

tabela$gender_p[1] = "---"
age = df$age[df$disease==out]
tabela$age_median[1] = sprintf("%.1f",median(age))
tabela$age_iqr[1] = sprintf("%.1f",IQR(age))
tabela$age_p[1] = "---"
tabela$min[1] = sprintf("%.1f",min(age))
tabela$max[1] = sprintf("%.1f",max(age))
for(i in 2:length(tabela$disease)){
  out = tabela$disease[i]
  print(out)
  df$out = NA
  df$out[df$disease=="Healthy"] = 0
  df$out[df$disease==out] = 1
  df$out = as.factor(df$out)
  a = data.frame(table(df$sex[df$disease==out]))
  tabela$gender_me[i] = paste0(a$Freq[a$Var1=="Female"],"(",sprintf("%0.1f",100*a$Freq[a$Var1=="Female"]/tabela$N[i]),"%)")
  tabela$gender_p[i] = sprintf("%0.3f",chisq.test(x=df$sex,y=df$out,simulate.p.value = TRUE)$p.value)
  age = df$age[df$disease==out]
  tabela$age_median[i] = sprintf("%.1f",median(age))
  tabela$age_iqr[i] = sprintf("%.1f",IQR(age))
  tabela$age_p[i] = sprintf("%.3f",wilcox.test(age~out,data = df)$p.value)
  tabela$min[i] = sprintf("%.1f",min(age))
  tabela$max[i] = sprintf("%.1f",max(age))
}
write.csv(tabela,"./result/table/desc.csv",row.names = FALSE,quote = FALSE)
