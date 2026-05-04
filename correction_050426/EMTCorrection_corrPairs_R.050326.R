### Finalized EMT correction reproduction in R, protein pairs
## 050326

### define, source functions
# uniprot to gene name
uniprotToGeneNames<- function(uniProtDB) {
  grabHeads <- getAnnot(uniProtDB)
  unlistedHeads <- unlist(lapply(grabHeads, `[[`, 1))
  splitHeads <- strsplit(unlistedHeads, " +" )
  sp_noGN <- splitHeads[-which(grepl(pattern = "GN=.*", splitHeads))]
  sp_noGN_upIds <- unlist(lapply(strsplit(unlist(lapply(sp_noGN, `[[`, 1)), "\\|"), `[[`, 2))
  splitHeads <- splitHeads[which(grepl(pattern = "GN=.*", splitHeads))]
  
  procMap <- data.frame(
    Protein = unlist(lapply(strsplit(unlist(lapply(splitHeads, `[[`, 1)), "\\|"), `[[`, 2)),
    Gene = gsub("GN=", "", unlist(lapply(splitHeads, grep, pattern = "GN=.*", value = TRUE)))
  )
  
  return(procMap)
}

# correlation convenience
corrVect_fromCorrMat <- function(mat1, mat2) {
  
  diag(mat1) <- NA
  diag(mat2) <- NA
  
  correlation_results <- data.frame(Prot = rownames(mat1), Cor = NA_real_)
  
  for (i in seq_len(nrow(correlation_results))) {
    correlation_results$Cor[i] <- stats::cor(
      mat1[i, ],
      mat2[i, ],
      use = "pairwise.complete.obs"
    )
  }
  
  correlation_results
}

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
plot_pair_by_rank <- function(pair_rank = 1, perm_table = perm_results_update_GN) {
  
  dark_blue <- grDevices::rgb(0, 0.2, 0.6)  # darker than MATLAB default
  
  pair_i <- perm_table %>%
    dplyr::arrange(dplyr::desc(abs_delta_r)) %>%
    dplyr::slice(pair_rank)
  
  right_prot <- pair_i$right_protein
  left_prot  <- pair_i$left_protein
  
  right_lab <- dplyr::coalesce(pair_i$right_Gene, right_prot)
  left_lab  <- dplyr::coalesce(pair_i$left_Gene, left_prot)
  
  plot_df <- tibble::tibble(
    x = as.numeric(dat[right_prot, ]),
    y = as.numeric(dat[left_prot, ]),
    day = cell_meta$day
  ) %>%
    dplyr::filter(!is.na(x), !is.na(y))
  
  ## compute correlation per day
  cor_df <- plot_df %>%
    dplyr::group_by(day) %>%
    dplyr::summarise(
      r = stats::cor(x, y, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(day)
  
  ## build facet labels like MATLAB
  facet_labs <- setNames(
    paste0("Day ", c(0, 3, 9), ", ρ = ", sprintf("%.2f", cor_df$r)),
    cor_df$day
  )
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    
    ## points
    ggplot2::geom_point(
      shape = 16,
      color = dark_blue,
      size = 2.6,
      alpha = 0.9
    ) +
    
    ## facet with custom labels
    ggplot2::facet_wrap(
      ~ day,
      nrow = 1,
      labeller = ggplot2::as_labeller(facet_labs)
    ) +
    
    ## labels
    ggplot2::labs(
      x = right_lab,
      y = left_lab
    ) +
    
    ## styling (no grid, MATLAB-ish)
    my_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      strip.text = ggplot2::element_text(size = 20, face = "bold"),
      axis.title = ggplot2::element_text(size = 20),
      axis.text  = ggplot2::element_text(size = 14)
    ) +
    
    ggplot2::coord_cartesian(expand = TRUE)
}

# significance testing via permutation 
Rcpp::sourceCpp("perm_diff_cor_both_update.cpp")

## load data
# Uniprot to Gene name mapping
uniProtDB <-  read.fasta(file = "uniprot-proteome_UP000005640.fasta")

# data matrix and cell mapping table
EMT_prot_unimp <- as.matrix(read.table(file = "EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.txt", sep = "\t", header = TRUE, stringsAsFactors = F))

EMT_cellIdtoTimepoint <-read.table(file = "cellIDToTimepoint.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

## preprocess
Prot_to_Gene_map <- uniprotToGeneNames(uniProtDB)

dat <- EMT_prot_unimp
cell_meta <- data.table::as.data.table(EMT_cellIdtoTimepoint)

stopifnot(all(colnames(dat) %in% cell_meta$id))

cell_meta <- cell_meta[match(colnames(dat), id)]
cell_meta[, day := factor(celltype, levels = c("d0", "d3", "d9"))]

n_obs <- rowSums(!is.na(dat))
names(n_obs) <- rownames(dat)

## corrs
cor_by_day <- lapply(levels(cell_meta$day), function(d) {
  Hmisc::rcorr(t(dat[, cell_meta$day == d, drop = FALSE]))
})

names(cor_by_day) <- levels(cell_meta$day)

r1 <- cor_by_day$d0$r
r2 <- cor_by_day$d3$r
r3 <- cor_by_day$d9$r


## correlation vector analysis
r1_r3_cva <- corrVect_fromCorrMat(r1, r3) %>%
  dplyr::rename(
    Var1 = Prot,
    correlation_between_r = Cor
  ) %>%
  dplyr::mutate(
    n_obs = n_obs[Var1]
  )


## threshold changed and stable tails
left_tail <- r1_r3_cva %>%
  #  dplyr::filter(correlation_between_r < 0.2, n_obs > 200) %>%
  dplyr::filter(correlation_between_r < 0.1, n_obs > 200) %>%
  dplyr::pull(Var1)

right_tail <- r1_r3_cva %>%
  #  dplyr::filter(correlation_between_r > 0.4, n_obs > 400) %>%
  dplyr::filter(correlation_between_r > 0.35, n_obs > 400) %>%
  dplyr::pull(Var1)


## rank and order proteins in tails by magnitude
r1_left_unord <- r1[left_tail, left_tail, drop = FALSE]
r1_right_unord <- r1[right_tail, right_tail, drop = FALSE]

left_ord <- order(colSums(r1_left_unord^2, na.rm = TRUE), decreasing = TRUE)
right_ord <- order(colSums(r1_right_unord^2, na.rm = TRUE), decreasing = TRUE)

left_tail_ord <- left_tail[left_ord]
right_tail_ord <- right_tail[right_ord]

## matrices for selected proteins
tail_prots <- unique(c(right_tail_ord, left_tail_ord))

d0_tailsProtsMat <- t(dat[, cell_meta$day == "d0", drop = FALSE])[, tail_prots, drop = FALSE]
d9_tailsProtsMat <- t(dat[, cell_meta$day == "d9", drop = FALSE])[, tail_prots, drop = FALSE]

stopifnot(identical(colnames(d0_tailsProtsMat), colnames(d9_tailsProtsMat)))


## relevant computations and structure using tails to set up for permutation testing
cross_pair_df <- tidyr::expand_grid(
  right_protein = right_tail_ord,
  left_protein = left_tail_ord
) %>%
  dplyr::mutate(
    r_day0 = r1[cbind(right_protein, left_protein)],
    r_day9 = r3[cbind(right_protein, left_protein)],
    delta_r = r_day9 - r_day0,
    abs_delta_r = abs(delta_r)
  ) %>%
  dplyr::arrange(dplyr::desc(abs_delta_r))


## optional: permutation test on cross-tail pairs
perm_results_update <- cross_pair_df %>%
  dplyr::mutate(
    res = purrr::map2(
      right_protein,
      left_protein,
      ~ perm_diff_cor_both_cpp(
        x1 = d0_tailsProtsMat[, .x],
        y1 = d0_tailsProtsMat[, .y],
        x2 = d9_tailsProtsMat[, .x],
        y2 = d9_tailsProtsMat[, .y],
        N = 100000L,
        do_trick = TRUE,
        return_nulls = FALSE,
        seed = 123L
      )
    ),
    p_perm_direct = purrr::map_dbl(res, "p_direct"),
    p_perm_trick = purrr::map_dbl(res, "p_trick"),
    D_obs_cpp = purrr::map_dbl(res, "D_obs"),
    n1_eff = purrr::map_dbl(res, "n1_eff"),
    n2_eff = purrr::map_dbl(res, "n2_eff")
  ) %>%
  dplyr::select(-res) %>%
  dplyr::mutate(
    FDR_perm_direct = stats::p.adjust(p_perm_direct, method = "BH"),
    FDR_perm_trick = stats::p.adjust(p_perm_trick, method = "BH")
  )

perm_results_update_GN <- perm_results_update %>% 
  left_join(Prot_to_Gene_map %>% rename(right_Gene = Gene), by = c("right_protein" = "Protein")) %>% 
  left_join(Prot_to_Gene_map %>% rename(left_Gene = Gene), by = c("left_protein" = "Protein"))

perm_results_update_GN <- perm_results_oldcpp_thresh2

### plotit

p1 <- plot_pair_by_rank(1, perm_table = perm_results_update_GN %>% dplyr::filter(FDR_perm_trick < 0.1))
p2 <- plot_pair_by_rank(11, perm_table = perm_results_update_GN %>% dplyr::filter(FDR_perm_trick < 0.1))


ggplot2::ggsave(filename = file.path("/EMTCorrect/EMT_CorrPair_CECR2_TAGLN.pdf"), plot = p1, width = 15, height = 6, units = "in", device = grDevices::cairo_pdf)
ggplot2::ggsave(filename = file.path("/EMT_Correct/EMT_CorrPair_YWHAZ_SHMT2.pdf"), plot = p2, width = 15, height = 6, units = "in", device = grDevices::cairo_pdf)

