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
# Figure 6
#A violin plot 
data = readRDS(file = "./data.rds")
df = data$df_p
prot = data$lab
df = df %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
  ))
df$out = NA
df$out[df$disease=="Healthy"] = "Healthy"
i=2
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 0] ="untreated" 
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 1] ="treated"

df$out = factor(df$out,levels = c("Healthy","treated","untreated"))

res = read.csv(paste0("./mark_p/" ,dis$var[i] ,".csv"))
res = res[res$p_value<=0.05,]
res = res$col_p[1:10]

colu = c("out",res)
aux = df[,colu]
aux = aux[!is.na(aux$out),]
colnames(aux)[2:11] = get_protein_names(colnames(aux)[2:11],prot)
res = get_protein_names(res,prot)
aux = gather(aux, marks,value,res)
aux = aux[!is.na(aux$value), ]
aux$marks = new_labels(aux$marks)
aux$marks = factor(aux$marks, levels=new_labels(res))
fig_a=aux %>% ggplot(aes(x=out,y=value,fill = out))+
  facet_wrap( ~ marks,scales = "free_y",ncol = 2 )+
  geom_violin(scale = "width",trim = FALSE,draw_quantiles = c(0.5))+
  geom_jitter(width = 0.2,height = 0,size=1)+
  #scale_x_discrete(limits = dis$var)+
  #scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values=c(dis$color[1],"treated"="pink2","untreated"=dis$color[i][[1]]),limits = c("Healthy","treated","untreated"))+
  labs(x=element_blank(),y = "Quantification")+
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

#figure B PCA ALL
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$out = NA
df$out[df$disease=="Healthy"] = "Healthy"
i=2
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 0] ="untreated" 
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 1] ="treated"
df$out = factor(df$out,levels = c("Healthy","treated","untreated"))
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
df$disease = NULL
df = df[!is.na(df$out),]
col = colnames(df)
col = col[col!="out"]
df_pca =  prcomp(df[, col], scale. = TRUE)
df_pca <- data.frame(PC1 = df_pca$x[, 1], PC2 = df_pca$x[, 2], Disease = as.character(df$out))
aux = df_pca

p <- aux %>% ggplot(aes(x=PC1, y=PC2))+
  geom_point(aes(color=Disease),size=2)+
  scale_color_manual(values=c(dis$color[1],"treated"="pink3","untreated"=dis$color[i][[1]]),limits = c("Healthy","treated","untreated"))+
  geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
fig_b= ggExtra::ggMarginal(p,groupColour = TRUE, groupFill = TRUE)

# figure C 
# figure f PCA best marks
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$out = NA
df$out[df$disease=="Healthy"] = "Healthy"
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 0] ="untreated" 
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 1] ="treated"

df$out = factor(df$out,levels = c("Healthy","treated","untreated"))
df = df[!is.na(df$out),]
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
df$disease = NULL
i=2
res = read.csv(paste0("./mark_p/" ,dis$var[i] ,".csv"))
res = res$col_p[1:dis$size[i]]
aux = df
y = aux$out
aux[,"out"] = NULL
x = aux[,res]
pc = prcomp(x,scale. = TRUE)
data = as.data.frame(pc$x[,1:2])
data$Disease = y
aux =  as.data.frame(pc$rotation[,1:2])[res[1:10],]
row.names(aux) = NULL
aux$Disease = get_protein_names(res[1:10],prot)
aux = aux[aux$Disease %in% c(aux$Disease[1],aux$Disease[2]),]
p <- data %>% ggplot(aes(x=PC1, y=PC2))+
  geom_point(aes(color=Disease),size=2)+
  scale_color_manual(values=c(dis$color[1],"treated"="pink3","untreated"=dis$color[i][[1]]),limits = c("Healthy","treated","untreated"))+
  geom_hline(yintercept=0, linetype=2,linewidth=1,color="darkgray")+geom_vline(xintercept=0, linetype=2,linewidth=1,color="darkgray")+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font))
fig1_pca <- p + 
  geom_richtext(data=aux, aes(x=PC1*20+c(0,-2), y=PC2*20+c(-1,1), label=new_labels(Disease)), size = 2.5, color="black",family=font,fill=fill_alpha("bisque2",0.7))+
  geom_segment(data=aux,aes(x=0, y=0, xend=PC1*20, yend=PC2*20), arrow=arrow(length=unit(0.3,"cm")), color="black",linewidth=1)+
  theme_bw(base_family = font)+
  theme(legend.position = "none",axis.title.x = element_text(size = 10,family = font),axis.title.y = element_text(size = 10,family = font),
        axis.text.x = element_text(size=8,family = font),axis.text.y = element_text(size=8,family = font),plot.caption = element_markdown(size = 10))

fig_c=ggExtra::ggMarginal(fig1_pca,groupColour = TRUE, groupFill = TRUE)


# figure D dist all
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$out = NA
df$out[df$disease=="Healthy"] = "Healthy"
i=2
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 0] ="untreated" 
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 1] ="treated"
df$out = factor(df$out,levels = c("Healthy","treated","untreated"))
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
df$disease = NULL
df = df[!is.na(df$out),]
col = colnames(df)
col = col[col!="out"]
x = df[,col]
pc = prcomp(x,scale. = TRUE)
x = as.data.frame(pc$x)

h_t = centroid_dist(x,y,c("Healthy","treated"))
h_nt = centroid_dist(x,y,c("Healthy","untreated"))
t_nt = centroid_dist(x,y,c("treated","untreated"))

fig_d = plot_triangle_dist(
  labels = c("Healthy","treated","untreated"),
  colors = c(dis$color[1][[1]],"pink3",dis$color[i][[1]]),
  d12 = h_t,        # H—T
  d13 = h_nt,       # H—NT
  d23 = t_nt        # T—NT
)

# figure E dist best
i=2
data = readRDS(file = "./data_nor.rds")
df = data$df_p
res = read.csv(paste0("./mark_p/" ,dis$var[i] ,".csv"))
res = res$col_p[1:dis$size[i]]
df$out = NA
df$out[df$disease=="Healthy"] = "Healthy"

df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 0] ="untreated" 
df$out[df$disease == dis$var[i] & !is.na(df$treated) & df$treated == 1] ="treated"
df$out = factor(df$out,levels = c("Healthy","treated","untreated"))
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
df$disease = NULL
df = df[!is.na(df$out),]

x = df[,res]
pc = prcomp(x,scale. = TRUE)
x = as.data.frame(pc$x)

h_t = centroid_dist(x,y,c("Healthy","treated"))
h_nt = centroid_dist(x,y,c("Healthy","untreated"))
t_nt = centroid_dist(x,y,c("treated","untreated"))

fig_e = plot_triangle_dist(
  labels = c("Healthy","treated","untreated"),
  colors = c(dis$color[1][[1]],"pink3",dis$color[i][[1]]),
  d12 = h_t,        # H—T
  d13 = h_nt,       # H—NT
  d23 = t_nt        # T—NT
)

# Painel
p1 = ggarrange(fig_b,fig_d,fig_c,fig_e,ncol = 4,labels=c("B","","C",""),widths = c(1,0.8,1,0.8))
p2 = ggarrange(fig_a,labels=c("A"))
p = ggarrange(p2,p1,nrow = 2,heights = c(4,1.5))
#p2 <- p2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

pdf("./fig6.pdf", family=font,width = 11,height = 12)
print(p)
dev.off()

