library(tidyverse)

color_eu = c("Healthy" ="#66CD00",
             "Autoinflammation of unknown origin" = "#009ACD",
             "Still's disease"="#FF4500",
             "FMF"="#9A32CD",
             "Behcet"="#F5E400")
#####################################################
dis = list(var = c("Healthy","Autoinflammation of unknown origin","Still's disease","FMF","Behcet"),
           size = c(NA,40,50,75,54),
           color = color_eu[1:5],
           ref1 = c(NA,"CD38+ CD8+ T Cells in CD8+ T Cells",
                    "CD38+ CD8+ T Cells in CD8+ T Cells",
                    "NK Cells in Live Cells",
                    "NK Cells in Live Cells"),
           ref2 = c(NA,"IgM+ IgD+ Memory B Cells in Memory B Cells",
                    "BAFF-Rhigh Naive B Cells in Naive B Cells",
                    "CD40+ CD21- B Cells in CD21- B Cells",
                    "CD94+ CD16+ NK Cells in CD16+ NK Cells"),
           pos_refx1 = c(NA,0.5,0.4,1,0.5),
           pos_refy1 =c(NA,0,1.7,1,1),
           pos_refx2 = c(NA,0.8,0.4,0.9,0.25),
           pos_refy2 =c(NA,1.5,0,-0.2,0))

print("OR")
###########################################################
data = readRDS(file = "./data_nor.rds")
df_f = data$df_f
df_p = data$df_p
col_f = colnames(df_f)[6:ncol(df_f)]
col_p = colnames(df_p)[6:ncol(df_p)]
#facs
for(d in 2:5){
  df_f$out = NA
  df_f$out[df_f$disease=="Healthy"] = 0
  df_f$out[df_f$disease==dis$var[d]] = 1
  colu = c("sex","age","out",col_f)
  aux = df_f[!is.na(df_f$out),colu]
  colnames(aux) = c("sex","age","out",paste0("v",1:length(col_f)))
  or = rep(NA,length(col_f))
  or_l = rep(NA,length(col_f))
  or_h = rep(NA,length(col_f))
  p_value =rep(NA,length(col_f))
  for(i in 1:length(col_f)){
    print(paste0(d," ",i))
    if(d==4){
      model = glm(paste0("out~v",i," + age"),data = aux,family = binomial()) 
    }else{
      model = glm(paste0("out~v",i," + sex + age"),data = aux,family = binomial())
    }
    or[i] = exp(coef(model))[paste0("v",i)][[1]]
    ci = exp(confint(model))[paste0("v",i),]
    or_l[i] = ci[[1]]
    or_h[i] = ci[[2]]
    p_value[i] = coef(summary(model))[paste0("v",i),"Pr(>|z|)"]
  }
  tabela = data.frame(col_f,or,or_l,or_h,p_value,ind=1:length(col_f))
  tabela$mag = NA
  tabela$mag[tabela$or>=1] = tabela$or[tabela$or>=1]
  tabela$mag[tabela$or<1] = 1/tabela$or[tabela$or<1]
  tab1 = tabela[tabela$p_value<=0.05,]
  tab2 = tabela[tabela$p_value>0.05,]
  tab1 = tab1[order(tab1$mag,decreasing = TRUE),]
  tab2 = tab2[order(tab2$mag,decreasing = TRUE),]
  tabela = rbind(tab1,tab2)
  tabela$p_value_adj = p.adjust(tabela$p_value, method = "BH")
  ####################################################
  write.csv(tabela,file=paste0("./mark_f/",dis$var[d] ,".csv"),row.names = FALSE,quote = FALSE)
}
#protein
for(d in 2:5){
  df_p$out = NA
  df_p$out[df_p$disease=="Healthy"] = 0
  df_p$out[df_p$disease==dis$var[d]] = 1
  colu = c("sex","age","out",col_p)
  aux = df_p[!is.na(df_p$out),colu]
  colnames(aux) = c("sex","age","out",paste0("v",1:length(col_p)))
  or = rep(NA,length(col_p))
  or_l = rep(NA,length(col_p))
  or_h = rep(NA,length(col_p))
  p_value =rep(NA,length(col_p))
  for(i in 1:length(col_p)){
    print(paste0(d," ",i))
    if(d==4){
      model = glm(paste0("out~v",i," + age"),data = aux,family = binomial()) 
    }else{
      model = glm(paste0("out~v",i," + sex + age"),data = aux,family = binomial())
    }
    or[i] = exp(coef(model))[paste0("v",i)][[1]]
    ci = exp(confint(model))[paste0("v",i),]
    or_l[i] = ci[[1]]
    or_h[i] = ci[[2]]
    p_value[i] = coef(summary(model))[paste0("v",i),"Pr(>|z|)"]
  }
  tabela = data.frame(col_p,or,or_l,or_h,p_value,ind=1:length(col_p))
  tabela$mag = NA
  tabela$mag[tabela$or>=1] = tabela$or[tabela$or>=1]
  tabela$mag[tabela$or<1] = 1/tabela$or[tabela$or<1]
  tab1 = tabela[tabela$p_value<=0.05,]
  tab2 = tabela[tabela$p_value>0.05,]
  tab1 = tab1[order(tab1$mag,decreasing = TRUE),]
  tab2 = tab2[order(tab2$mag,decreasing = TRUE),]
  tabela = rbind(tab1,tab2)
  tabela$p_value_adj = p.adjust(tabela$p_value, method = "BH")
  ########################################################
  write.csv(tabela,file=paste0("./mark_p/",dis$var[d] ,".csv"),row.names = FALSE,quote = FALSE)
}

######################################################################################
saveRDS(data$df_f, file = "./data_nor_f.rds")
saveRDS(data$df_p, file = "./data_nor_p.rds")
