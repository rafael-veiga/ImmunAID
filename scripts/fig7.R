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


####################################################################################
# figure 7
## painel auc RF
RF_auc2 <- read.csv("./RF_auc2_p.csv")
data = readRDS(file = "./data_nor.rds")
df = data$df_p
df$sample = NULL
df$sex = NULL
df$age = NULL
df$treated = NULL
col = colnames(df)
col = col[col!="disease"]
df_pca =  prcomp(df[, col], scale. = TRUE)
x_pca = as.data.frame(df_pca$x)
y_pca = df$disease
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
      dist = centroid_dist(x_pca, as.factor(y_pca), c(doencas[i], doencas[j]))
      dist_lab <- sprintf("%.2f", dist)
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
              plot.margin = unit(c(1,1,1,1),"mm"))+
        annotate(
          "label",
          x = pc1_min + 0.03 * (pc1_max - pc1_min),
          y = pc2_min + 0.05 * (pc2_max - pc2_min),
          label = dist_lab,
          size = 5, family = font,
          hjust = 0, vjust = 0,
          label.size = 0.2,       # borda discreta
          fill = "gray60", color = "darkblue"
        )
      plot_list[[paste(doencas[i], doencas[j], sep="_")]] <- p
        } else { # Diagonal
      disease_name <- doencas[j]

      p <- ggplot() + 
        theme_void(base_family = font) +
        scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0,0))

      # setas
      if (j != 1){
        p <- p +
          geom_segment(aes(x=0.5, y=0.5, xend=0.5, yend=0.9),
                      arrow=arrow(length=unit(.5,'cm')),
                      color=color_eu[disease_name], linewidth=3) +
          geom_segment(aes(x=0.5, y=0.5, xend=0.1, yend=0.5),
                      arrow=arrow(length=unit(.5,'cm')),
                      color=color_eu[disease_name], linewidth=3)
      }
      if (j != length(doencas)){
        p <- p +
          geom_segment(aes(x=0.5, y=0.5, xend=0.5, yend=0.1),
                      arrow=arrow(length=unit(.5,'cm')),
                      color=color_eu[disease_name], linewidth=3) +
          geom_segment(aes(x=0.5, y=0.5, xend=0.9, yend=0.5),
                      arrow=arrow(length=unit(.5,'cm')),
                      color=color_eu[disease_name], linewidth=3)
      }

      # texto do rótulo (pré-avaliado!)
      lab_text <- gsub("Inflammation of unknown origin",
                      "Inflammation of\nunknown origin",
                      disease_name)

      # rótulo por cima das setas, sem aes dinâmico
      p <- p +
        geom_label(
          data = data.frame(x=0.5, y=0.5, lab=lab_text),
          aes(x=x, y=y, label=lab),
          inherit.aes = FALSE,
          size = 6, family = font,
          hjust = 0.5, vjust = 0.5,
          color = color_eu[disease_name],
          fill = "white",
          label.size = 0.6
        )

      plot_list[[paste(disease_name, disease_name, sep="_")]] <- p
    }

  }
}

pdf("./fig7.pdf",width = 20,height = 16)
print(do.call(gridExtra::grid.arrange, c(plot_list, ncol = length(doencas))))
dev.off()
