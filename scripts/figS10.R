#####################################################################################################################
# BASE

library(ape)
library(WCluster)
library(ggExtra)
library(tidyverse)
library("ggpubr")
library(ggtree)
library(patchwork)
library(gghalves)
library(randomForest)
library(ggrepel)
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
             "Autoinflammation of unknown origin" = "#009ACD",
             "Still's disease"="#FF4500",
             "FMF"="#9A32CD",
             "Behcet"="#F5E400")
#####################################################
dis = list(var = c("Healthy","Autoinflammation of unknown origin","Still's disease","FMF","Behcet"),
           size = c(NA,40,82,55,61),
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

centroid_dist <- function(x, y, groups) {
  stopifnot(length(groups) == 2)
  # Filtrar apenas as instâncias dos dois grupos
  sel <- y %in% groups
  Xsub <- x[sel, , drop = FALSE]
  ysub <- droplevels(y[sel])
  # Centróides de cada grupo
  mu1 <- colMeans(Xsub[ysub == groups[1], , drop = FALSE])
  mu2 <- colMeans(Xsub[ysub == groups[2], , drop = FALSE])
  # Distância euclidiana entre os centróides
  dist <- sqrt(sum((mu1 - mu2)^2))
  return(dist)
}

plot_triangle_dist <- function(labels, d12, d13, d23,
                               colors = c("red","green","blue"),
                               edge_label_digits = 2,
                               point_size = 5, edge_size = 0.8,
                               dashed = TRUE, normalize = TRUE) {
  stopifnot(length(labels) == 3, length(colors) == 3)
  
  # normalizar distâncias se pedido
  if (normalize) {
    s <- d12 + d13 + d23
    d12 <- d12 / s
    d13 <- d13 / s
    d23 <- d23 / s
  }
  
  # checagem de triângulo
  eps <- 1e-9
  if ((d12 + d13 <= d23 + eps) || (d12 + d23 <= d13 + eps) || (d13 + d23 <= d12 + eps)) {
    stop("As distâncias não formam um triângulo (violam desigualdade do triângulo).")
  }
  
  # A=labels[1], B=labels[2], C=labels[3]
  a <- d23; b <- d13; c <- d12
  xC <- (b^2 + c^2 - a^2) / (2*c)
  yC <- sqrt(max(b^2 - xC^2, 0))
  
  nodes <- data.frame(
    label = labels,
    x = c(0, c, xC),
    y = c(0, 0, yC),
    color = colors
  )
  
  edges <- data.frame(
    x    = c(0,     0,     c),
    y    = c(0,     0,     0),
    xend = c(c,     xC,    xC),
    yend = c(0,     yC,    yC),
    dist = c(d12,   d13,   d23)
  )
  edges$xm <- (edges$x + edges$xend)/2
  edges$ym <- (edges$y + edges$yend)/2
  
  # padding
  pad <- 0.15 * max(c(d12, d13, d23))
  xlim <- range(c(nodes$x, nodes$xend)) + c(-pad, pad)
  ylim <- range(c(nodes$y, nodes$yend)) + c(-pad, pad)
  
  ggplot() +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = edge_size,
                 linetype = if (dashed) "dashed" else "solid",
                 color = "grey40") +
    geom_point(data = nodes, aes(x = x, y = y), size = point_size) +
    geom_label(data = nodes, aes(x = x, y = y, label = label, fill = color),
               nudge_y = 0.03*max(c(d12,d13,d23)), label.size = 0.2, alpha = 0.9,
               color = "black") +
    geom_text(data = edges,
              aes(x = xm, y = ym, label = round(dist, edge_label_digits)),
              vjust = -0.6) +
    scale_fill_identity() +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, clip = "off") +
    theme_void(base_size = 12) +    # remove eixos e grids
    theme(legend.position = "none")
}

get_protein_names <- function(ids, prot_df) {
  idx <- match(ids, prot_df$`Majority protein IDs`)
  prot_df$`Protein names`[idx]
}


##################################################################################
#figure 5
##or A
d=3
res = read.csv(paste0("./mark_p/" ,dis$var[d] ,".csv"))
data = readRDS(file = "./data.rds")
prot = data$lab
prot$`Gene names` = NULL
res = merge.data.frame(prot,res,by.x="Majority protein IDs",by.y="col_p",all.x = FALSE,all.y = TRUE)
res$`Majority protein IDs` = NULL
colnames(res)[1] = "col_p"

top_labels <- res %>%
  filter(p_value < 0.05) %>%
  arrange(desc(mag)) %>%
  slice(1:10)

# Volcano plot
fig_a = ggplot(res, aes(x = log2(or), y = -log10(p_value))) +
  geom_point(aes(color = p_value <= 0.05), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("TRUE"="red", "FALSE"="gray50")) +
  geom_text_repel(
    data = top_labels,
    aes(label = col_p),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "black",
    segment.size = 0.6,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed") # seta
  ) +
  theme_minimal(base_size = 12) +
  labs(x = "log2(OR)", y = "-log10(p-value)", 
       color = "p ≤ 0.05")

# figure B auc_N 
data = read.csv("./auc_n_p.csv")
fig_b = data %>% filter(disease==dis$var[d]) %>%
  ggplot(aes(x=num,y=auc))+
  geom_point(size=1,color="gray")+
  geom_line()+
  #geom_ribbon(aes(ymax = auc_h,ymin = auc_l),fill = dis$color[d],alpha=0.4)+
  geom_smooth(method = "loess",span=0.3,fill=dis$color[d],color=dis$color[d])+
  labs(x = "Markers number",y="AUC")+
  scale_y_continuous(limits = c(0.7,NA)) +
  scale_x_continuous(limits = c(1,NA)) +
  geom_vline(xintercept=48,color = "darkred") +
  annotate("text",x=48-4,y=0.73,label=paste0("n = ",48),size=4,family=font,angle = 90,color = "darkred")+
  theme_bw(base_family  = font)+
  theme(axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))

#figure C violin plot
data = readRDS(file = "./data.rds")
df = data$df_p
df = df %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
  ))

df$disease = factor(df$disease,levels = dis$var)
res = read.csv(paste0("./mark_p/" ,dis$var[d] ,".csv"))
res = res[res$p_value<=0.05,]
res = res$col_p[1:10]
colu = c("disease",res)
aux = df[,colu]
colnames(aux)[2:11] = get_protein_names(colnames(aux)[2:11],prot)
res = get_protein_names(res,prot)
aux = gather(aux, marks,value,res)
aux = aux[!is.na(aux$value), ]
aux$marks = factor(aux$marks, levels=res)
fig_c=aux %>% ggplot(aes(x=disease,y=value,fill = disease))+
  facet_wrap( ~ marks,scales = "free_y",ncol = 2 )+
  geom_violin(scale = "width",trim = FALSE,draw_quantiles = c(0.5))+
  geom_jitter(width = 0.2,height = 0,size=1)+
  scale_x_discrete(limits = dis$var)+
  #scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values=dis$color,limits = dis$var)+
  labs(x=element_blank(),y = "Log2 LFQ intensity")+
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

# figure D auc curve
data = read.csv(file = "./auc_curv_p.csv")
aux  = data %>% filter(disease==dis$var[d])
label =paste0(sprintf(aux$auc[1], fmt = '%#.2f')," (",sprintf(aux$auc_l[1], fmt = '%#.2f')," : ",sprintf(aux$auc_h[1], fmt = '%#.2f'),")")
fig_d = aux %>%  ggplot(aes(x=fpr,y=tpr))+
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

# figure E PCA ALL
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
col = colnames(df)
col = col[col!="disease"]
df = df[df$disease %in% c("Healthy",dis$var[d]),]
df_pca =  prcomp(df[, col], scale. = TRUE)
df_pca <- data.frame(PC1 = df_pca$x[, 1], PC2 = df_pca$x[, 2], Disease = as.character(df$disease))
aux = df_pca

p <- aux %>% ggplot(aes(x=PC1, y=PC2))+
  geom_point(aes(color=Disease),size=2)+
  scale_color_manual(values=color_eu,limits = c("Healthy",dis$var[d]))+
  geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
fig_e= ggExtra::ggMarginal(p,groupColour = TRUE, groupFill = TRUE)

# figure f PCA best marks
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
res = read.csv(paste0("./mark_p/" ,dis$var[d] ,".csv"))
res = res$col_p[1:dis$size[d]]
aux = df
aux$out = NA
aux$out[aux$disease=="Healthy"] = "Healthy"
aux$out[aux$disease==dis$var[d]] = dis$var[d]
filx = !is.na(aux$out)
y = aux$out[filx]
aux[,"out"] = NULL
x = aux[filx,res]
pc = prcomp(x,scale. = TRUE)
data = as.data.frame(pc$x[,1:2])
data$Disease = y
aux =  as.data.frame(pc$rotation[,1:2])[res[1:10],]
row.names(aux) = NULL
aux$Disease = get_protein_names(res[1:10],prot)
aux = aux[aux$Disease %in% c(aux$Disease[1],aux$Disease[2]),]
p <- data %>% ggplot(aes(x=PC1, y=PC2))+
  geom_point(aes(color=Disease),size=2)+
  scale_color_manual(values=color_eu,limits = c("Healthy",dis$var[d]))+
  geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
fig1_pca <- p + 
  geom_richtext(data=aux, aes(x=PC1*20, y=PC2*20+c(-0.5,2), label=Disease), size = 2.5,hjust=c(dis$pos_refx1[d],dis$pos_refx2[d]),vjust=c(dis$pos_refy1[d],dis$pos_refy2[d]), color="black",family=font,fill=fill_alpha("bisque2",0.7))+
  geom_segment(data=aux,aes(x=0, y=0, xend=PC1*20, yend=PC2*20), arrow=arrow(length=unit(0.3,"cm")), color="black",linewidth=1)+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font),plot.caption = element_text(size = 10))

fig_f=ggExtra::ggMarginal(fig1_pca,groupColour = TRUE, groupFill = TRUE)

# Painel
p1 = ggarrange(fig_a,fig_c,nrow = 2,labels=c("A","C"),heights = c(2, 7))
p2 = ggarrange(fig_b,fig_d,fig_e,fig_f,ncol = 1,labels=c("B","D","E","F"),widths = c(1),heights = c(1,1,1,1))
p2 <- p2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p = ggarrange(p1,p2,ncol=2,widths = c(2.6,1))

pdf("./figS10.pdf", family=font,width = 13.5,height = 12)
print(p)
dev.off()
