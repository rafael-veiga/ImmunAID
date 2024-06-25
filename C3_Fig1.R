source("functions_pre.R")
library(ape)
library(WCluster)
library(ggExtra)
library(tidyverse)
library("ggpubr")
library(ggtree)
library(patchwork)
#library(factoextra)
library(randomForest)
#library(latex2exp)
library(ggtext)
library(extrafont)
library("RColorBrewer")
set.seed(12345)
#font_import(paths = "c:/windows/Fonts/")
#loadfonts(device = "win")
#windowsFonts()
font ="serif"

#####################################################################################################################
color_eu = c("Healthy" ="#66CD00",
             "Inflammation of unknown origin" = "#009ACD",
             "AOSD"="#FF4500",
             "FMF"="#9A32CD",
             "Behcet"="#F5E400")
#####################################################
dis = list(var = c("Healthy","Inflammation of unknown origin","AOSD","FMF","Behcet"),
           size = c(NA,40,82,55,61),
           color = color_eu[1:5],
           ref1 = c(NA,"CD38+ CD8+ T Cells in CD8+ T Cells",
                    "CD38+ CD8+ T Cells in CD8+ T Cells",
                    "NK Cells in Live Cells",
                    "NK Cells in Live Cells"),
           ref2 = c(NA,"IgM+ IgD+ Memory B Cells in Memory B Cells",
                    "HLA-DR+ CD4+ T Cells in CD4+ T Cells",
                    "CD40+ CD21- B Cells in CD21- B Cells",
                    "CD94+ CD16+ NK Cells in CD16+ NK Cells"),
           pos_refx1 = c(NA,0.5,0.4,1,0.5),
           pos_refy1 =c(NA,0,1.7,1,1),
           pos_refx2 = c(NA,0.8,0.4,0.9,0.25),
           pos_refy2 =c(NA,1.5,0,-0.2,0))
new_labels = function(vet){
  aux = sapply(as.character(vet), function(x) {
    x = gsub("\\+", "<sup>+</sup>", x)
    x = gsub("\\-", "<sup>-</sup>", x)
    x = gsub("Non<sup>-</sup>Switched", "Non-Switched", x)
    x = gsub("low", "<sup>low</sup>", x)
    x = gsub("high", "<sup>high</sup>", x)
    x = gsub("dim", "<sup>dim</sup>", x)
    x = gsub("bright", "<sup>bright</sup>", x)
    return(x)
  })
  return(aux)
}
get_peach_color <- function(value, min_val, max_val) {
  num_colors <- 9
  colors <- brewer.pal(num_colors, "YlOrRd")
  normalized_value <- (value - min_val) / (max_val - min_val)
  color_index <- round(normalized_value * (num_colors - 1)) + 1
  return(colors[color_index])
}
#####################################################################################

# process
## Odds Ration
print("OR")
df = readRDS(file = "./pos_data/data_nor.rds")
col = colnames(df)
var = col[get_fil(col,type = "pop_in_pop")]
for(d in 2:5){
  df$out = NA
  df$out[df$disease=="Healthy"] = 0
  df$out[df$disease==dis$var[d]] = 1
  colu = c("sex","age","out",var)
  aux = df[!is.na(df$out),colu]
  #standartize trans
  aux[,var] = scale(aux[,var])
  colnames(aux) = c("sex","age","out",paste0("v",1:length(var)))
  or = rep(NA,length(var))
  or_l = rep(NA,length(var))
  or_h = rep(NA,length(var))
  p_value =rep(NA,length(var))
  for(i in 1:length(var)){
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
  tabela = data.frame(var,or,or_l,or_h,p_value,ind=1:length(var))
  tabela$mag = NA
  tabela$mag[tabela$or>=1] = tabela$or[tabela$or>=1]
  tabela$mag[tabela$or<1] = 1/tabela$or[tabela$or<1]
  tab1 = tabela[tabela$p_value<=0.05,]
  tab2 = tabela[tabela$p_value>0.05,]
  tab1 = tab1[order(tab1$mag,decreasing = TRUE),]
  tab2 = tab2[order(tab2$mag,decreasing = TRUE),]
  tabela = rbind(tab1,tab2)
  write.csv(tabela,file=paste0("./result/table/marks/",dis$var[d] ,".csv"),row.names = FALSE,quote = FALSE)
}

######################################################################################

#run "auc_curv.py"

########################################################################################
# Figs


#or
fig_a = list()
fig_a[[1]] = NA
for(i in 2:5){
  res = read.csv(paste0("./result/table/marks/" ,dis$var[i] ,".csv"))
  aux = res[res$p_value<=0.05,]
  aux = aux[1:10,]
  aux = aux %>% mutate(out = factor(var, levels=aux$var))
  fig_a[[i]] = aux %>% ggplot(aes(y=out, x=or, xmin=or_l, xmax=or_h)) +
    geom_pointrange()+
    scale_y_discrete (limits = rev(as.vector(aux$out)),labels=new_labels(rev(as.vector(aux$out))))+
    scale_x_continuous(trans = scales::log_trans(),breaks = c(0.1,0.5,1,5,10,50,100),limits=c(0.1,NA))+
    geom_vline(xintercept=1, lty=2, linewidth =1)+
    geom_errorbar(width=0.5, cex=1)+ # Makes whiskers on the range (more aesthetically pleasing)
    geom_point(shape = 15, size = 2,color=dis$color[i])+
    labs(y=element_blank(),x = "OR")+
    theme_classic(base_family = font)+
    theme(axis.title.x = element_text(size = 10,family = font),
          axis.text.x = element_text(size = 8,family=font),
          axis.text.y = ggtext::element_markdown(size = 10,family=font),
          plot.margin = unit(c(2,2,0,2),"mm"))
}

#violin
fig_c = list()
fig_c[[1]] = NA
df = readRDS(file = "./pos_data/data_raw.rds")
df = df[df$disease %in% dis$var,]
df$disease = factor(df$disease,levels = dis$var)

for(i in 2:5){
  res = read.csv(paste0("./result/table/marks/" ,dis$var[i] ,".csv"))
  res = res[res$p_value<=0.05,]
  res = res$var[1:10]
  
  colu = c("disease",res)
  aux = df[,colu]
  aux = gather(aux, marks,value,res)
  aux = aux[!is.na(aux$value), ]
  aux$marks = new_labels(aux$marks)
  aux$marks = factor(aux$marks, levels=new_labels(res))
  fig_c[[i]]=aux %>% ggplot(aes(x=disease,y=value,fill = disease))+
    facet_wrap( ~ marks,scales = "free_y",ncol = 2 )+
    geom_violin(scale = "width",trim = FALSE,draw_quantiles = c(0.5))+
    geom_jitter(width = 0.2,height = 0,size=1)+
    scale_x_discrete(limits = dis$var)+
    scale_y_continuous(limits = c(0,NA))+
    scale_fill_manual(values=dis$color,limits = dis$var)+
    labs(x=element_blank(),y = "Frequency")+
    theme_bw(base_family =font)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(size=8,family = font),
          axis.title.y = element_text(size=10,family = font),
          legend.title=element_blank(),
          legend.text = element_text(size=11,family = font),
          strip.text = ggtext::element_markdown(size = 10,family =font),
          strip.background = element_rect(fill = "bisque"),
          legend.position = "top",
          legend.justification='left',
          legend.direction='horizontal',
          plot.margin = unit(c(0,2,2,2),"mm"))
}





# AUC vc N
fig_b = list()
fig_b[[1]] = NA
data = read.csv("./result/table/auc_n.csv")
for(d in 2:5){
  fig_b[[d]] = data %>% filter(disease==dis$var[d]) %>%
    ggplot(aes(x=num,y=auc))+
    geom_point(size=1,color="gray")+
    geom_line()+
    #geom_ribbon(aes(ymax = auc_h,ymin = auc_l),fill = dis$color[d],alpha=0.4)+
    geom_smooth(method = "loess",span=0.3,fill=dis$color[d],color=dis$color[d])+
    labs(x = "Markers number",y="AUC")+
    scale_y_continuous(limits = c(0.7,NA)) +
    scale_x_continuous(limits = c(1,NA)) +
    geom_vline(xintercept=dis$size[d],color = "darkred") +
    annotate("text",x=dis$size[d]-4,y=0.73,label=paste0("n = ",dis$size[d]),size=4,family=font,angle = 90,color = "darkred")+
    theme_bw(base_family  = font)+
    theme(axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
          axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
}


# AUC curve
data = read.csv(file = "./result/table/auc_curv/auc_curv.csv")
fig_d = list()
fig_d[[1]] = NA
for(d in 2:5){
  aux  = data %>% filter(disease==dis$var[d])
  label =paste0(sprintf(aux$auc[1], fmt = '%#.2f')," (",sprintf(aux$auc_l[1], fmt = '%#.2f')," : ",sprintf(aux$auc_h[1], fmt = '%#.2f'),")")
fig_d[[d]] = aux %>%  ggplot(aes(x=fpr,y=tpr))+
  geom_line(color = dis$color[d],linewidth=1)+
  geom_ribbon(aes(ymax = tpr_h,ymin = tpr_l),fill = dis$color[d],alpha=0.4)+
  geom_richtext(color=get_peach_color(aux$auc[1],0.50,1),fill="gray95",x=0.7,y=0.1,label=label,size=4,family = font,fontface ="plain")+
  geom_abline(slope = 1,intercept = 0,linetype=2,linewidth=1,color="gray")+
  scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_y_continuous(limits = c(0, 1),expand = c(0,0))+
  ylab("True positive rate")+
  xlab("False positive rate")+
  theme_bw(base_family = font)+
  theme(axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
}

#PCA all
fig_e = list()
fig_e[[1]] = NA
df = readRDS(file = "./pos_data/data_nor.rds")
df$id = NULL
df$sex = NULL
df$age = NULL
df$batch = NULL
col = colnames(df)
col = col[col!="disease"]
df_pca =  prcomp(df[, col], scale. = TRUE)
df_pca <- data.frame(PC1 = df_pca$x[, 1], PC2 = df_pca$x[, 2], Disease = as.character(df$disease))
for(i in 2:5){
  aux = df_pca
  aux = aux[aux$Disease %in% c("Healthy",dis$var[i]),]
  p <- aux %>% ggplot(aes(x=PC1, y=PC2))+
    geom_point(aes(color=Disease),size=2)+
    scale_color_manual(values=color_eu,limits = c("Healthy",dis$var[i]))+
    geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
    theme_bw(base_family = font)+
    theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
          axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
  fig_e[[i]]= ggExtra::ggMarginal(p,groupColour = TRUE, groupFill = TRUE)
}


#PCA best
fig_f = list()
fig_f[[1]] = NA
df = readRDS(file = "./pos_data/data_nor.rds")
for(i in 2:5){
  res = read.csv(paste0("./result/table/marks/" ,dis$var[i] ,".csv"))
  res = res$var[1:dis$size[i]]
  aux = df
  aux = df
  aux$out = NA
  aux$out[aux$disease=="Healthy"] = "Healthy"
  aux$out[aux$disease==dis$var[i]] = dis$var[i]
  filx = !is.na(aux$out)
  y = aux$out[filx]
  aux[,"out"] = NULL
  x = aux[filx,res]
  pc = prcomp(x,scale. = TRUE)
  data = as.data.frame(pc$x[,1:2])
  data$Disease = y
  aux =  as.data.frame(pc$rotation[,1:2])[res[1:10],]
  row.names(aux) = NULL
  aux$Disease = res[1:10]
  aux = aux[aux$Disease %in% c(dis$ref1[i],dis$ref2[i]),]
  p <- data %>% ggplot(aes(x=PC1, y=PC2))+
    geom_point(aes(color=Disease),size=2)+
    scale_color_manual(values=color_eu,limits = c("Healthy",dis$var[i]))+
    geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
    theme_bw(base_family = font)+
    theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
          axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
  fig1_pca <- p + 
      geom_richtext(data=aux, aes(x=PC1*20, y=PC2*20, label=new_labels(Disease)), size = 2.5,hjust=c(dis$pos_refx1[i],dis$pos_refx2[i]),vjust=c(dis$pos_refy1[i],dis$pos_refy2[i]), color="black",family=font,fill=fill_alpha("bisque2",0.7))+
    geom_segment(data=aux,aes(x=0, y=0, xend=PC1*20, yend=PC2*20), arrow=arrow(length=unit(0.3,"cm")), color="black",linewidth=1)+
    theme_bw(base_family = font)+
    theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
          axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font),plot.caption = element_markdown(size = 10))
  
  fig_f[[i]]=ggExtra::ggMarginal(fig1_pca,groupColour = TRUE, groupFill = TRUE)
  
}


### painel
p = list()
p[[1]] = NA
for(i in 2:5){
  p1 = ggarrange(fig_a[[i]],fig_c[[i]],nrow = 2,labels=c("A","C"),heights = c(2, 7))
  p2 = ggarrange(fig_b[[i]],fig_d[[i]],fig_e[[i]],fig_f[[i]],ncol = 1,labels=c("B","D","E","F"),widths = c(1),heights = c(1,1,1,1))
  p2 <- p2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  p[[i]] = ggarrange(p1,p2,ncol=2,widths = c(2.6,1))
}
pdf("fig1.pdf", family=font,width = 13.5,height = 12)
print(p[[2]])
dev.off()

pdf("fig2.pdf", family=font,width = 13.5,height = 12)
print(p[[3]])
dev.off()

pdf("fig3.pdf", family=font,width = 13.5,height = 12)
print(p[[4]])
dev.off()

pdf("fig4.pdf", family=font,width = 13.5,height = 12)
print(p[[5]])
dev.off()






###### fig aux
#or
fig_aux = list()
fig_aux[[1]] = NA
for(i in 2:5){
  res = read.csv(paste0("./result/table/marks/" ,dis$var[i] ,".csv"))
  aux = res[res$p_value<=0.05,]
  aux = aux[1:50,]
  aux = aux[!is.na(aux$or),]
  aux = aux %>% mutate(out = factor(var, levels=aux$var))
  fig_aux[[i]] = aux %>% ggplot(aes(y=out, x=or, xmin=or_l, xmax=or_h)) +
    geom_pointrange()+
    scale_y_discrete (limits = rev(as.vector(aux$out)),labels=new_labels(rev(as.vector(aux$out))))+
    scale_x_continuous(trans = scales::log_trans(),breaks = c(0.1,0.5,1,5,10,50,100),limits=c(0.1,NA))+
    geom_vline(xintercept=1, lty=2, linewidth =1,color="darkgray")+
    geom_errorbar(width=0.5, cex=1)+ # Makes whiskers on the range (more aesthetically pleasing)
    geom_point(shape = 15, size = 2,color=dis$color[i])+
    labs(y=element_blank(),x = "OR",title = dis$var[i])+
    theme_classic(base_family = font)+
    theme(axis.title.x = element_text(size = 10,family = font),
          axis.text.x = element_text(size = 8,family=font),
          axis.text.y = ggtext::element_markdown(size = 10,family=font),
          plot.margin = unit(c(2,2,0,2),"mm"),
          plot.title = element_text(size=18,family = font,color = dis$color[i],hjust = 0.3))
}
pdf("fig_aux1.pdf", family=font,width = 11,height = 12)
print(fig_aux[[2]])
dev.off()

pdf("fig_aux2.pdf", family=font,width = 11,height = 12)
print(fig_aux[[3]])
dev.off()

pdf("fig_aux3.pdf", family=font,width = 11,height = 12)
print(fig_aux[[4]])
dev.off()

pdf("fig_aux4.pdf", family=font,width = 11,height = 12)
print(fig_aux[[5]])
dev.off()

## painel auc RF
RF_auc2 <- read.csv("./result/RF_auc2.csv")
df = readRDS(file = "./pos_data/data_nor.rds")
df$id = NULL
df$sex = NULL
df$age = NULL
df$batch = NULL
col = colnames(df)
col = col[col!="disease"]
df_pca =  prcomp(df[, col], scale. = TRUE)
df_pca <- data.frame(PC1 = df_pca$x[, 1], PC2 = df_pca$x[, 2], Disease = df$disease)
doencas <- unique(c(RF_auc2$d1, RF_auc2$d2))
doencas = unique(c("Healthy",doencas))
plot_list <- vector("list", length(doencas)^2)
names(plot_list) <- paste(rep(doencas, each=length(doencas)), rep(doencas, times=length(doencas)), sep="_")
pc1_max = max(df_pca$PC1)
pc1_min = min(df_pca$PC1)
pc2_max = max(df_pca$PC2)
pc2_min = min(df_pca$PC2)
for (i in 1:length(doencas)) {
  for (j in 1:length(doencas)) {
    if (i < j) { # AUC
      plot_data <- subset(RF_auc2, d1 == doencas[i] & d2 == doencas[j])
      label=paste0(format(round(plot_data$auc[1],2), nsmall = 2)," (",format(round(plot_data$auc_l[1],2), nsmall = 2) ," : ",format(round(plot_data$auc_h[1],2), nsmall = 2),")")
      p <- ggplot(plot_data, aes(x=fpr, y=tpr)) +
        geom_line() +
        theme_bw(base_family = font)+
        geom_ribbon(aes(ymin=tpr_l, ymax=tpr_h), alpha=0.2) +
        geom_richtext(color=get_peach_color(plot_data$auc[1],0.50,1),fill="gray95",x=0.6,y=0.1,label=label,size=5,family = font,fontface ="plain")+
        geom_abline(slope = 1,intercept = 0,linetype=2,linewidth=1,color="lightblue")+
        scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
        scale_y_continuous(limits = c(0, 1),expand = c(0,0))+
        theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
              axis.text.x = element_text(size=10,family = font),axis.text.y = element_text(size=10,family = font),
              legend.title=element_blank(),plot.margin = unit(c(2,3,2,0),"mm"))
      plot_list[[paste(doencas[i], doencas[j], sep="_")]] <- p
    } else if (i > j) { # PCA
      aux <- df_pca[df_pca$Disease %in% c(doencas[i], doencas[j]), ]
      p <- ggplot(aux, aes(x=PC1, y=PC2, color=Disease)) +
        geom_point(aes(color=Disease)) +
        scale_color_manual(values=color_eu) +
        scale_x_continuous(limits = c(pc1_min,pc1_max)) + 
        scale_y_continuous(limits = c(pc2_min,pc2_max)) +
        geom_hline(yintercept=0, size=.2,color="darkgray",linetype=3)+geom_vline(xintercept=0, size=.2,color="darkgray",linetype=3)+
        theme_minimal(base_family  = font)+
        theme(legend.position="none",axis.title.x=element_blank(),
                                           legend.title=element_blank(),axis.title.y=element_blank(),
              axis.text.x=element_text(size=10,family = font),
              axis.text.y=element_text(size=10,family = font),
              plot.margin = unit(c(1,1,1,1),"mm"))
      plot_list[[paste(doencas[i], doencas[j], sep="_")]] <- p
    } else { # Diagonal
      disease_name <- doencas[j]
        p <- ggplot() + 
          theme_void(base_family  = font) +
          scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
          scale_y_continuous(limits = c(0, 1),expand = c(0,0))
          if(j!=1){
            p <- p+geom_segment(aes(x=0.5, y=0.5, xend=0.5, yend=0.9), arrow = arrow(length=unit(.5, 'cm')),color=color_eu[disease_name], lwd=3)+ 
              geom_segment(aes(x=0.5, y=0.5, xend=0.1, yend=0.5), arrow = arrow(length=unit(.5, 'cm')),color=color_eu[disease_name], lwd=3) 
          }
          if(j!=5){
            p <-p + geom_segment(aes(x=0.5, y=0.5, xend=0.5, yend=0.1), arrow = arrow(length=unit(.5, 'cm')),color=color_eu[disease_name], lwd=3)+
              geom_segment(aes(x=0.5, y=0.5, xend=0.9, yend=0.5), arrow = arrow(length=unit(.5, 'cm')),color=color_eu[disease_name], lwd=3)
          }
          
          p = p+annotate("label", x = 0.5, y = 0.5, label = gsub("Inflammation of unknown origin", "Inflammation of\nunknown origin", disease_name), 
                   size = 6, hjust = 0.5, vjust = 0.5, 
                   color = color_eu[disease_name]) # Ajuste o size conforme necess?rio e use a cor correspondente
        plot_list[[paste(disease_name, disease_name, sep="_")]] <- p
    }
  }
}

pdf("fig5.pdf",width = 18,height = 16)
print(do.call(gridExtra::grid.arrange, c(plot_list, ncol = length(doencas))))
dev.off()

## violin B
df = read.csv("./result/RF_imp2.csv")
res = as.character(df$vari)
res = res[res!="age"]
res = res[res!="sex"]
res = res[1:16]
df = readRDS(file = "./pos_data/data_raw.rds")
df = df[df$disease %in% dis$var,]
df$disease = factor(df$disease,levels = dis$var)
df = df %>% select(c("disease",res))
df = gather(df, marks,value,res)
df = df[!is.na(df$value),]


df$marks = new_labels(df$marks)
df$marks = factor(df$marks, levels=new_labels(res))
p2=df %>% ggplot(aes(x=disease,y=value,fill = disease))+
  facet_wrap( ~ marks,scales = "free_y",ncol = 2 )+
  geom_violin(scale = "width",trim = FALSE,draw_quantiles = c(0.5))+
  geom_jitter(width = 0.2,height = 0,size=1)+
  scale_x_discrete(limits = dis$var)+
  scale_y_continuous(limits = c(0,NA))+
  scale_fill_manual(values=dis$color,limits = dis$var)+
  labs(x=element_blank(),y = "Frequency")+
  theme_bw(base_family =font)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=8,family = font),
        axis.title.y = element_text(size=10,family = font),
        legend.title=element_blank(),
        legend.text = element_text(size=11,family = font),
        strip.text = ggtext::element_markdown(size = 10,family =font),
        strip.background = element_rect(fill = "bisque"),
        legend.position = "top",
        legend.justification='left',
        legend.direction='horizontal',
        plot.margin = unit(c(0,2,2,2),"mm"))

pdf("fig6.pdf",width = 9,height = 12)
print(p2)
dev.off()



###### relevant marks
## bars 
df = read.csv("./result/RF_imp2.csv")
df = df[1:50,]
df$vari = factor(df$vari,levels = rev(df$vari))
p1 = ggplot(df, aes(x = imp, y = vari,fill = imp)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_gradient("Importance",high="orange",low="darkred")+
  scale_y_discrete(label = new_labels(rev(df$vari)))+
  scale_x_continuous(limits = c(0, NA),expand = c(0,NA))+
  labs(x = "Importance",y=element_blank())+
  #theme(axis.text.x = element_text( hjust = 1, vjust = 1,size=8,family=font),axis.text.y = element_markdown(size=10,family = font),
  theme(axis.text.x = element_blank(),axis.text.y = element_markdown(size=10,family = font),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=10,family=font))

pdf("fig_aux5.pdf",width = 9,height = 12)
print(p1)
dev.off()


#### tree
df = read.csv("./result/RF_imp2.csv")
res = as.character(df$vari)
res = res[res!="age"]
res = res[res!="sex"]
res = res[1:20]
df = readRDS(file = "./pos_data/data_nor_t.rds")
df$id = NULL
df$sex = NULL
df$age = NULL
df$batch = NULL
tb = table(df$disease)
df = df %>% select(all_of(c("disease",res))) %>% group_by(disease) %>% summarise_all(mean)
df$disease = as.character(df$disease)
x = as.matrix(df[,colnames(df)!="disease"])
rownames(x) <- df$disease
df$w = NA
df$w[df$disease==names(tb)] = tb

#cor_matrix <- cor(t(df[,colnames(df)!="disease"]))
#dist_matrix=as.matrix(dist(df[,colnames(df)!="disease"],method = 'euclidean',diag = TRUE, upper = TRUE))
#rownames(dist_matrix) <- df$disease
#colnames(dist_matrix) <- df$disease
#dist_matrix <- as.dist(1 - cor_matrix)
hc <- Whclust(x,df$w)
tree <- as.phylo(hc)
tree$tip.label = df$disease
labdata = data.frame(label = df$disease)
#labdata$label[labdata$label=="Inflammation of unknown origin"] = "Inflammation of<br>unknown origin"

p2 <- ggtree(tree, layout = "rectangular") %<+% labdata +
  geom_tiplab(aes(label = label), offset = 0) +
  scale_x_continuous(limits = c(0,2050))+
  #geom_text(label="A",x=0.1,y=0.9,size=5,family=font)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.background  = element_rect(color = "black", fill = "lightgray", linewidth = 1))

# Adicionar as legendas
p2 <- add_labels(p2, labdata)
# p2 <- ggtree(tree, layout="rectangular") %<+% labdata
# p2 = p2 + geom_tiplab(aes(label = label, color = color), hjust = 0,offset = 0.5) +
#   scale_color_identity()+
#   theme(plot.margin = unit(c(1,1,1,1),"mm"))
final_plot <- p1 +
  inset_element(p2, left = 0.5, bottom = 0.05, right = 0.95, top = 0.25)

# Exibir o gr?fico final

