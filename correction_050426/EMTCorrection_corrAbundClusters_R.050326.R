### Finalized EMT correction reproduction in R
## 050226

## setup
# clear workspace
rm(list = ls())

# packages and functions
packages <- c("tidyverse", "ggplot2", "Hmisc", "data.table","proxy")
lapply(packages, require, character.only = TRUE)

# ploting
my_theme <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      plot.title = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 22),
      legend.key.size = unit(2.0, "lines"),
      legend.spacing.y = unit(1.0, "lines"),
      strip.text = element_text(size = 22)
    )
}
baseR_theme <- function() {
  par(
    cex.lab = 1.5,    
    cex.axis = 1.5,   
    cex.main = 1.5,   
    cex.sub = 1.2,    
    mar = c(5, 5, 4, 2),  
    las = 1,          
    tcl = -0.3,      
    mgp = c(2.5, 0.7, 0) 
  )
}
baseR_theme()

# load data
EMT_prot_unimp <- as.matrix(read.table(file = "EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.txt", sep = "\t", header = TRUE, stringsAsFactors = F))

EMT_cellIdtoTimepoint <-read.table(file = "cellIDToTimepoint.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

## preprocess
dat <- as.matrix(EMT_prot_unimp)
cell_meta <- data.table::as.data.table(EMT_cellIdtoTimepoint)

stopifnot(all(colnames(dat) %in% cell_meta$id))

cell_meta <- cell_meta[match(colnames(dat), id)]
cell_meta[, day := factor(celltype, levels = c("d0", "d3", "d9"))]
cell_meta[, day_num := as.integer(day)]

prot_ids <- rownames(dat)

## filter

keep_total <- rowSums(!is.na(dat)) > 30

keep_by_day <- sapply(levels(cell_meta$day), function(d) {
  rowSums(!is.na(dat[, cell_meta$day == d, drop = FALSE])) > 4
})

keep_prot <- keep_total & rowSums(keep_by_day) == 3

cat(sprintf("Selected Proteins: %d\n\n", sum(keep_prot)))

dat_f <- dat[keep_prot, , drop = FALSE]
prot_f <- prot_ids[keep_prot]

## corrs

cor_by_day <- lapply(levels(cell_meta$day), function(d) {
  Hmisc::rcorr(t(dat_f[, cell_meta$day == d, drop = FALSE]))$r
})

names(cor_by_day) <- levels(cell_meta$day)

r1 <- cor_by_day$d0
r2 <- cor_by_day$d3
r3 <- cor_by_day$d9

## corrDiffs, selection

rd <- r3 - r1

complete_idx <- which(colMeans(is.na(rd)) == 0)
rd_complete <- rd[complete_idx, complete_idx, drop = FALSE]

sum_r <- colSums(rd_complete^2)
thr50 <- stats::quantile(sum_r, 0.5, na.rm = TRUE)

sel <- sum_r > thr50

cat(sprintf("Number of selected proteins: %d\n\n", sum(sel)))

In <- complete_idx[sel]
prot_sel <- prot_f[In]

## cluster

rd_sel <- rd[In, In, drop = FALSE]

dist_cos <- proxy::dist(rd_sel, method = "cosine")
hc <- stats::hclust(dist_cos, method = "complete")
hc_order <- hc$order

X <- cbind(
  r1[In, In, drop = FALSE],
  r2[In, In, drop = FALSE],
  r3[In, In, drop = FALSE]
)

set.seed(1)
km <- stats::kmeans(X, centers = 3)

cluster <- km$cluster
cluster_order <- order(cluster)

Inx <- In[cluster_order]

## plot em

theme_set(my_theme())

draw_mat <- function(M, title, lim) {
  cols <- grDevices::colorRampPalette(c("blue", "white", "red"))(200)
  
  graphics::image(
    z = t(M[nrow(M):1, , drop = FALSE]),
    zlim = c(-lim, lim),
    col = cols,
    axes = FALSE,
    main = title
  )
}

upper_vec <- function(M) { M[upper.tri(M)] }
mean_offdiag <- function(M) { mean(upper_vec(M), na.rm = TRUE) }
## heatmaps

Lim1 <- 0.7
Lim2 <- 0.5

grDevices::pdf("/EMT_Correct/EMT_Global_correlation-Patterns_Kmeans-Ver2-1_file2.pdf", width = 12, height = 7)

graphics::par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))

draw_mat(r1[Inx, Inx], "Day 0", Lim1)
draw_mat(r2[Inx, Inx], "Day 3", Lim1)
draw_mat(r3[Inx, Inx], "Day 9", Lim1)

draw_mat(r2[Inx, Inx] - r1[Inx, Inx], "Day 3 - Day 0", Lim2)
draw_mat(r3[Inx, Inx] - r2[Inx, Inx], "Day 9 - Day 3", Lim2)
draw_mat(r3[Inx, Inx] - r1[Inx, Inx], "Day 9 - Day 0", Lim2)

grDevices::dev.off()

## Day 0 vs 9, corr scaling
ci <- In[cluster == 1]

corr_scaling_dt <- data.table::data.table(
  d0 = upper_vec(r1[ci, ci, drop = FALSE]),
  d9 = upper_vec(r3[ci, ci, drop = FALSE])
)

stats_r <- stats::cor(
  corr_scaling_dt$d0,
  corr_scaling_dt$d9,
  use = "pairwise.complete.obs"
)

p_corr_scaling <- ggplot2::ggplot(corr_scaling_dt, ggplot2::aes(d0, d9)) +
  ggplot2::geom_point(shape = 1) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 1) +
  ggplot2::annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("r = %.2f", stats_r),
    hjust = 1.2,
    vjust = 1.5,
    size = 5
  ) +
  ggplot2::labs(x = "Day 0", y = "Day 9")

p_corr_scaling

ggplot2::ggsave("/EMT_Correct/Corr-scaling_file2.pdf", p_corr_scaling, width = 9, height = 6)

## clusters corrs across doys

cMean <- sapply(1:3, function(k) {
  ci <- In[cluster == k]
  c(
    d0 = mean_offdiag(r1[ci, ci, drop = FALSE]),
    d3 = mean_offdiag(r2[ci, ci, drop = FALSE]),
    d9 = mean_offdiag(r3[ci, ci, drop = FALSE])
  )
})

colnames(cMean) <- paste0("c", 1:3)

cMean_dt <- data.table::as.data.table(cMean, keep.rownames = "day")
cMean_long <- data.table::melt(
  cMean_dt,
  id.vars = "day",
  variable.name = "cluster",
  value.name = "mean_correlation"
)

cMean_long[, day := factor(day, levels = c("d0", "d3", "d9"),
                           labels = c("Day 0", "Day 3", "Day 9"))]

p_mean_corr <- ggplot2::ggplot(
  cMean_long,
  ggplot2::aes(day, mean_correlation, group = cluster, color = cluster)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(shape = 21, fill = "white", size = 4, stroke = 1.2) +
  ggplot2::labs(x = NULL, y = "Mean correlation", color = NULL)

p_mean_corr

ggplot2::ggsave("/EMT_Correct/Mean-cluster-correlation_file2.pdf", p_mean_corr, width = 9, height = 6)

## cluster abudances over time

aMean <- sapply(1:3, function(k) {
  ci <- In[cluster == k]
  c(
    d0 = mean(rowMeans(dat_f[ci, cell_meta$day == "d0", drop = FALSE], na.rm = TRUE), na.rm = TRUE),
    d3 = mean(rowMeans(dat_f[ci, cell_meta$day == "d3", drop = FALSE], na.rm = TRUE), na.rm = TRUE),
    d9 = mean(rowMeans(dat_f[ci, cell_meta$day == "d9", drop = FALSE], na.rm = TRUE), na.rm = TRUE)
  )
})

colnames(aMean) <- paste0("c", 1:3)

aMean_dt <- data.table::as.data.table(aMean, keep.rownames = "day")
aMean_long <- data.table::melt(
  aMean_dt,
  id.vars = "day",
  variable.name = "cluster",
  value.name = "mean_abundance"
)

aMean_long[, day := factor(day, levels = c("d0", "d3", "d9"),
                           labels = c("Day 0", "Day 3", "Day 9"))]

p_mean_abund <- ggplot2::ggplot(
  aMean_long,
  ggplot2::aes(day, mean_abundance, group = cluster, color = cluster)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(shape = 21, fill = "white", size = 4, stroke = 1.2) +
  ggplot2::labs(x = NULL, y = "Mean abundance", color = NULL)

p_mean_abund

ggplot2::ggsave("/EMT_Correct/Mean-cluster-abundance_file2.pdf", p_mean_abund, width = 9, height = 6)

## cluster corrBundance over time

mean_compare_dt <- merge(
  cMean_long,
  aMean_long,
  by = c("day", "cluster")
)

r_mean_corr_overall <- stats::cor(
  mean_compare_dt$mean_abundance,
  mean_compare_dt$mean_correlation,
  use = "pairwise.complete.obs"
)

p_mean_compare <- ggplot2::ggplot(
  mean_compare_dt,
  ggplot2::aes(mean_correlation, mean_abundance, color = cluster)
) +
  ggplot2::geom_point(shape = 21, fill = "white", size = 5, stroke = 1.2) +
  ggplot2::annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf("overall r = %.2f", r_mean_corr_overall),
    hjust = 1.1,
    vjust = -0.5,
    size = 6
  )+
  ggplot2::labs(
    x = "Mean correlation",
    y = "Mean abundance",
    color = NULL
  )

p_mean_compare

ggplot2::ggsave("/EMT_Correct/Mean-clustercorr-vs-abundance_file2.pdf", p_mean_compare, width = 9, height = 6)
