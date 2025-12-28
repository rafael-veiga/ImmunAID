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

wilcox_p <- function(x, g, a, b) {
  xa <- x[g == a]; xb <- x[g == b]
  xa <- xa[is.finite(xa)]; xb <- xb[is.finite(xb)]
  if (length(xa) < 2 || length(xb) < 2) return(NA_real_)
  tryCatch(wilcox.test(xa, xb)$p.value, error = function(e) NA_real_)
}

########################################################################################
########################################################################################
# Figure 8
## violin B (usando p-ajustados da pasta mark_p, estilo Figure 4)
########################################################################################

# df de importância
df  = read.csv("./RF_imp2_p.csv")
res = as.character(df$vari)
res = res[res != "age"]
res = res[res != "sex"]
res = res[1:16]

# dados principais
data = readRDS(file = "./data.rds")
df   = data$df_p
prot = data$lab

# treated = 0 para Healthy (mantém para o violin)
df$treated[df$disease == "Healthy"] = 0

# renomeia doenças
df = df %>% 
  dplyr::mutate(
    disease = dplyr::case_when(
      disease == "Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
      disease == "AOSD"                           ~ "Still's disease",
      TRUE                                        ~ disease
    )
  )

# só doenças de interesse listadas em dis$var
df = df[df$disease %in% dis$var, ]

# renomeia colunas de proteína (para nomes "bonitos")
#colnames(df)[6:ncol(df)] = get_protein_names(colnames(df)[6:ncol(df)], prot)
#res = get_protein_names(res, prot)

# ordenação dos grupos
df$disease = factor(df$disease, levels = dis$var)

# long format
df = df %>% dplyr::select(c("disease", "treated", res))
df = tidyr::gather(df, marks, value, res)
df = df[!is.na(df$value), ]

################################################################################
# ANEXANDO p-ajustados da pasta mark_p (col_p, p_value_adj)
################################################################################

ref    <- "Healthy"
others <- setdiff(dis$var, ref)           # todas as outras doenças
n_others <- length(others)

# vamos assumir que os arquivos em ./mark_p/ têm nomes = nome da doença (como em Figure 4)
# e que a coluna de marcador é "col_p", com p-ajustado em "p_value_adj"

res_aux <- vector("list", length = n_others)

for (a in seq_len(n_others)) {
  aux_p <- read.csv(paste0("./mark_p/", others[a], ".csv"))
  # mantém só os marcadores que estão em 'res'
  aux_p <- aux_p[aux_p$col_p %in% res, c("col_p", "p_value_adj")]
  colnames(aux_p)[2] <- paste0("p", a)    # p1, p2, ..., pN
  res_aux[[a]] <- aux_p
}

# tabela base com todos os marcadores
df_aux <- tibble::tibble(col_p = res)

# mescla p1..pN por marcador
for (i in seq_len(n_others)) {
  df_aux <- merge.data.frame(df_aux, res_aux[[i]], by = "col_p")
}

# junta os p-ajustados no df longo
df <- merge.data.frame(df, df_aux, by.x = "marks", by.y = "col_p")

# só agora aplicamos os rótulos bonitos pros painéis
df$marks = new_labels(df$marks)
df$marks = factor(df$marks, levels = new_labels(res))

aux <- df   # objeto que vai pro ggplot e pro cálculo de posições

################################################################################
# PREPARANDO ANOTAÇÕES: estilo Figure 4 (stat_pvalue_manual)
################################################################################

# mapa de comparações: p1 = Healthy vs others[1], ..., pN
p_cols <- paste0("p", seq_len(n_others))

comparisons_map <- tibble::tibble(
  cmp    = p_cols,
  group1 = ref,
  group2 = others[seq_len(n_others)]
)

# estatística por painel para posicionar as barras
panel_stats8 <- aux %>%
  dplyr::group_by(marks) %>%
  dplyr::summarise(
    y_top = as.numeric(stats::quantile(value, 0.98, na.rm = TRUE)),
    y_low = as.numeric(stats::quantile(value, 0.02, na.rm = TRUE)),
    y_rng = y_top - y_low,
    .groups = "drop"
  ) %>%
  dplyr::mutate(y_rng = ifelse(y_rng <= 0, 1, y_rng))

# df de anotações no formato esperado por stat_pvalue_manual
p_anno8 <- aux %>%
  dplyr::distinct(marks, dplyr::across(dplyr::all_of(p_cols))) %>%  # 1 linha por marks
  tidyr::pivot_longer(
    cols      = dplyr::all_of(p_cols),
    names_to  = "cmp",
    values_to = "p.adj"
  ) %>%
  dplyr::left_join(comparisons_map, by = "cmp") %>%
  dplyr::left_join(panel_stats8,   by = "marks") %>%
  dplyr::group_by(marks) %>%
  dplyr::arrange(marks, group2, .by_group = TRUE) %>%
  dplyr::mutate(
    idx = dplyr::row_number(),
    label = dplyr::if_else(
      is.na(p.adj),
      NA_character_,
      dplyr::if_else(
        p.adj < 0.0001,
        "p < 0.0001",
        paste0("p = ", formatC(p.adj, format = "f", digits = 4))
      )
    ),
    # empilha barras acima do dado (igual “espírito” da Figure 4)
    y.position = y_top + (idx - 1) * 0.30 * y_rng
  ) %>%
  dplyr::ungroup() %>%
  # plota só significantes
  dplyr::filter(!is.na(p.adj), p.adj < 0.05)

aux$marks = get_protein_names(aux$marks, prot)
res = get_protein_names(res, prot)
aux$marks = factor(aux$marks,levels = res)
p_anno8$marks = get_protein_names(p_anno8$marks, prot)
panel_stats8$marks = get_protein_names(panel_stats8$marks, prot)
################################################################################
# PLOT: violinos iguais ao original + barras/p-valor estilo Figure 4
################################################################################

fig_3 <- aux %>% 
  ggplot2::ggplot(ggplot2::aes(x = disease, y = value, fill = disease)) +
  ggplot2::facet_wrap(~ marks, scales = "free_y", ncol = 2) +
  gghalves::geom_half_violin(
    data = aux %>% dplyr::filter(treated == 0), side = "l",
    trim = FALSE, scale = "width", width = 1
  ) +
  gghalves::geom_half_violin(
    data = aux %>% dplyr::filter(treated == 1), side = "r",
    trim = FALSE, scale = "width", width = 1
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(shape = factor(treated), group = factor(treated)),
    position = ggplot2::position_jitterdodge(
      dodge.width   = 0.5,
      jitter.width  = 0.6,
      jitter.height = 0
    ),
    size = 1
  ) +
  # barra de média dashed e fina
  ggplot2::stat_summary(
    ggplot2::aes(group = factor(treated)),
    fun      = mean, fun.min = mean, fun.max = mean,
    geom     = "crossbar", width = 0.3,
    position = ggplot2::position_dodge(width = 0.8),
    color    = "black",
    linetype = "dashed",
    size     = 0.3
  ) +
  ggplot2::scale_shape_manual(
    values = c("0" = 19, "1" = 4),
    labels = c("0" = "untreated", "1" = "treated"),
    name   = "Treated"
  ) +
  ggplot2::scale_fill_manual(values = dis$color, limits = dis$var) +
  ggplot2::scale_x_discrete(limits = dis$var) +
  ggplot2::labs(x = NULL, y = "Log2 LFQ intensity") +
  ggplot2::theme_bw(base_family = font) +
  ggplot2::theme(
    axis.title.x      = ggplot2::element_blank(),
    axis.text.x       = ggplot2::element_blank(),
    axis.ticks.x      = ggplot2::element_blank(),
    axis.text.y       = ggplot2::element_text(size = 8, family = font),
    axis.title.y      = ggplot2::element_text(size = 10, family = font),
    legend.title      = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = 11, family = font),
    strip.text        = ggtext::element_markdown(size = 10, family = font),
    strip.background  = ggplot2::element_rect(fill = "bisque"),
    legend.position   = "top",
    legend.justification = "left",
    legend.direction  = "horizontal",
    plot.margin       = grid::unit(c(0, 2, 2, 2), "mm")
  )

# mesma “parte de cima” da Figure 4: espaço extra + barras e p-valor
fig_3 <- fig_3 +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggpubr::stat_pvalue_manual(
    p_anno8,
    label      = "label",
    tip.length = 0.01,
    size       = 2,
    hide.ns    = TRUE
  ) +
  ggplot2::theme(
    plot.margin   = grid::unit(c(6, 2, 2, 2), "mm"),
    panel.spacing = grid::unit(3, "mm")
  )

pdf("./fig8.pdf", width = 9, height = 12)
print(fig_3)
dev.off()
