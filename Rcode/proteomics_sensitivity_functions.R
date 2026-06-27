# proteomics_sensitivity_functions.R

library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(impute)
library(imputeLCMD)
library(pwr)

calc_freq_safe <- function(genes, df, samples, gene_col="GeneID") {
  expr_mat <- as.matrix(df[, samples, drop = FALSE])
  sapply(genes, function(g) {
    idx <- which(df[[gene_col]] == g)
    if (length(idx) == 1) {
      mean(!is.na(expr_mat[idx, ])) * 100
    } else {
      NA_real_
    }
  })
}

prepare_proteomics_matrix <- function(expr_dat_fil, gene_col="GeneID"){
  mat <- as.data.frame(expr_dat_fil)
  rownames(mat) <- make.unique(mat[[gene_col]])
  mat <- mat[, setdiff(colnames(mat), gene_col)]
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat
}

calculate_missingness <- function(mat_raw, phen_dat,
                                  sample_col="SampleID",
                                  group_col="SubGroup",
                                  control_group="HC",
                                  case_group="PDAC"){
  hc <- phen_dat[[sample_col]][phen_dat[[group_col]]==control_group]
  pd <- phen_dat[[sample_col]][phen_dat[[group_col]]==case_group]
  all <- phen_dat[[sample_col]]

  data.frame(
    Protein=rownames(mat_raw),
    n_HC_detected=rowSums(!is.na(mat_raw[,hc,drop=FALSE])),
    n_PDAC_detected=rowSums(!is.na(mat_raw[,pd,drop=FALSE])),
    n_ALL_detected=rowSums(!is.na(mat_raw[,all,drop=FALSE])),
    missing_HC_percent=rowMeans(is.na(mat_raw[,hc,drop=FALSE]))*100,
    missing_PDAC_percent=rowMeans(is.na(mat_raw[,pd,drop=FALSE]))*100,
    missing_ALL_percent=rowMeans(is.na(mat_raw[,all,drop=FALSE]))*100
  )
}

filter_by_detection_frequency <- function(mat_raw, phen_dat,
                                          sample_col="SampleID",
                                          group_col="SubGroup",
                                          control_group="HC",
                                          case_group="PDAC",
                                          min_detect_freq=0.5){
  hc <- phen_dat[[sample_col]][phen_dat[[group_col]]==control_group]
  pd <- phen_dat[[sample_col]][phen_dat[[group_col]]==case_group]
  keep <- rowSums(!is.na(mat_raw[,hc,drop=FALSE]))>=length(hc)*min_detect_freq |
          rowSums(!is.na(mat_raw[,pd,drop=FALSE]))>=length(pd)*min_detect_freq
  mat_raw[keep,,drop=FALSE]
}

run_limma_dep <- function(expr_mat, phen_dat,
                          sample_col="SampleID",
                          group_col="SubGroup",
                          control_group="HC",
                          case_group="PDAC",
                          method_name="limma"){
  expr_mat <- expr_mat[,phen_dat[[sample_col]],drop=FALSE]
  group <- factor(phen_dat[[group_col]], levels=c(control_group,case_group))
  design <- model.matrix(~0+group)
  colnames(design) <- c(control_group,case_group)
  cont <- makeContrasts(contrasts=paste0(case_group,"-",control_group), levels=design)
  fit <- eBayes(contrasts.fit(lmFit(expr_mat,design),cont))
  topTable(fit,coef=1,n=Inf,adjust.method="BH") %>%
    rownames_to_column("Protein") %>%
    mutate(Method=method_name)
}

run_lcn2_sensitivity <- function(mat_raw, phen_dat,
                                 target = "LCN2",
                                 sample_col = "SampleID",
                                 group_col = "SubGroup",
                                 control_group = "HC",
                                 case_group = "PDAC",
                                 min_detect_freq = 0.5,
                                 knn_k = 5,
                                 minprob_q = 0.01,
                                 log_transform = TRUE,
                                 pseudocount = 1) {
  
  # ----------------------------
  # 1. filter by detection frequency
  # ----------------------------
  hc_samples   <- phen_dat[[sample_col]][phen_dat[[group_col]] == control_group]
  pdac_samples <- phen_dat[[sample_col]][phen_dat[[group_col]] == case_group]
  
  keep <- rowSums(!is.na(mat_raw[, hc_samples, drop = FALSE])) >= length(hc_samples) * min_detect_freq |
    rowSums(!is.na(mat_raw[, pdac_samples, drop = FALSE])) >= length(pdac_samples) * min_detect_freq
  
  mat_fil <- mat_raw[keep, , drop = FALSE]
  
  # ----------------------------
  # 2. log transform
  # ----------------------------
  if (log_transform) {
    mat_log <- log2(mat_fil + pseudocount)
  } else {
    mat_log <- mat_fil
  }
  
  # ----------------------------
  # 3. No imputation
  # ----------------------------
  res_noimp <- run_limma_dep(mat_log, phen_dat, method_name = "No imputation")
  
  # ----------------------------
  # 4. KNN imputation
  # ----------------------------
  knn_mat <- impute::impute.knn(mat_log, k = knn_k)$data
  
  res_knn <- run_limma_dep(knn_mat, phen_dat, method_name = "KNN imputation")
  
  # ----------------------------
  # 5. MinProb (MNAR)
  # ----------------------------
  minprob_mat <- imputeLCMD::impute.MinProb(mat_log, q = minprob_q)
  
  res_minprob <- run_limma_dep(minprob_mat, phen_dat, method_name = "MinProb imputation")
  
  # ----------------------------
  # 6. LCN2 sensitivity table
  # ----------------------------
  lcn2_sensitivity <- bind_rows(res_knn, res_noimp, res_minprob) %>%
    dplyr::filter(Protein == target) %>%
    dplyr::select(Method, Protein, logFC, AveExpr, t, P.Value, adj.P.Val, B)
  
  # ----------------------------
  # 7. raw (complete-case) test
  # ----------------------------
  lcn2_vec <- mat_log[target, ]
  
  df <- data.frame(
    Expression = as.numeric(lcn2_vec),
    Group = phen_dat[[group_col]]
  )
  
  df <- df[!is.na(df$Expression), ]
  
  wt <- wilcox.test(Expression ~ Group, data = df)
  
  lcn2_complete_case <- data.frame(
    Method = "Complete-case Wilcoxon",
    Protein = target,
    logFC = mean(df$Expression[df$Group == case_group]) -
      mean(df$Expression[df$Group == control_group]),
    AveExpr = mean(df$Expression),
    t = NA,
    P.Value = wt$p.value,
    adj.P.Val = NA,
    B = NA
  )
  
  lcn2_sensitivity <- bind_rows(lcn2_sensitivity, lcn2_complete_case)
  
  # ----------------------------
  # 8. return
  # ----------------------------
  list(
    filtered_matrix = mat_fil,
    log_matrix = mat_log,
    knn_mat = knn_mat,
    minprob_mat = minprob_mat,
    res_knn = res_knn,
    res_noimp = res_noimp,
    res_minprob = res_minprob,
    target_sensitivity = lcn2_sensitivity,
    target_raw_df = df
  )
}

count_dep_direction <- function(res, fc_col = "logFC", p_col = "adj.P.Val", p_cut = 0.05, fc_cut = 2) {
  
  res <- res %>%
    mutate(
      Regulation = case_when(
        !!sym(p_col) < p_cut & !!sym(fc_col) > log2(fc_cut) ~ "Up",
        !!sym(p_col) < p_cut & !!sym(fc_col) < -log2(fc_cut) ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  data.frame(
    Up = sum(res$Regulation == "Up", na.rm = TRUE),
    Down = sum(res$Regulation == "Down", na.rm = TRUE),
    Not_significant = sum(res$Regulation == "Not significant", na.rm = TRUE)
  )
}

plot_gene_imputation_comparison <- function(gene,
                                            phen_dat,
                                            mat_list,
                                            out_prefix,
                                            fig_plot,
                                            methods_names = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  
  if (is.null(methods_names)) {
    methods_names <- names(mat_list)
  }
  
  plot_results <- list()
  
  for (i in seq_along(mat_list)) {
    
    mat <- mat_list[[i]]
    
    # -----------------------------
    # extract gene expression
    # -----------------------------
    gene_df <- mat %>%
      data.frame() %>%
      tibble::rownames_to_column("GeneID") %>%
      filter(GeneID == gene)
    
    if (nrow(gene_df) == 0) next
    
    plot_df <- gene_df %>%
      pivot_longer(cols = any_of(phen_dat$SampleID),
                   names_to = "SampleID",
                   values_to = "Expression") %>%
      mutate(Expression = as.numeric(Expression)) %>%
      inner_join(phen_dat, by = "SampleID") %>%
      mutate(Expression = replace_na(Expression, 0),
             SubGroup = factor(SubGroup, levels = c("HC", "PDAC")))
    
    # -----------------------------
    # boxplot
    # -----------------------------
    p_box <- ggplot(plot_df, aes(x = SubGroup, y = Expression, fill = SubGroup)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
      scale_fill_manual(values = c("HC" = "#3C8DBC", "PDAC" = "#D73027")) +
      labs(title = paste0(gene, " (", methods_names[i], ")"),
           y = "Log2 Intensity", x = "") +
      theme_bw() +
      stat_compare_means(method = "wilcox.test",
                         comparisons = list(c("HC", "PDAC")),
                         label = "p.signif") +
      theme(legend.position = "none")
    
    # save boxplot
    pdf(paste0(fig_plot, out_prefix, "_", gene, "_", methods_names[i], "_boxplot.pdf"),
        width = 3, height = 4)
    print(p_box)
    dev.off()
    
    # -----------------------------
    # barplot
    # -----------------------------
    plot_df_bar <- plot_df %>%
      arrange(SubGroup, desc(Expression)) %>%
      mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    
    p_bar <- ggplot(plot_df_bar, aes(x = SampleID, y = Expression, fill = SubGroup)) +
      geom_bar(stat = "identity", width = 0.7) +
      facet_grid(. ~ SubGroup, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = c("HC" = "#3C8DBC", "PDAC" = "#D73027")) +
      labs(title = paste0(gene, " (", methods_names[i], ")"),
           y = "Log10 Intensity", x = "Samples") +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        legend.position = "top"
      )
    
    pdf(paste0(fig_plot, out_prefix, "_", gene, "_", methods_names[i], "_barplot.pdf"),
        width = 6, height = 4)
    print(p_bar)
    dev.off()
    
    plot_results[[methods_names[i]]] <- plot_df
  }
  
  return(plot_results)
}

estimate_target_power <- function(target_raw_df,
                                  group_col = "Group",
                                  control_group = "HC",
                                  case_group = "PDAC",
                                  expression_col = "Expression",
                                  target = "LCN2") {
  
  library(dplyr)
  library(pwr)
  
  # -----------------------------
  # 1. summarize group stats
  # -----------------------------
  stat_df <- target_raw_df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n = n(),
      mean = mean(.data[[expression_col]], na.rm = TRUE),
      sd = sd(.data[[expression_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  n1 <- stat_df$n[stat_df[[group_col]] == control_group]
  n2 <- stat_df$n[stat_df[[group_col]] == case_group]
  
  m1 <- stat_df$mean[stat_df[[group_col]] == control_group]
  m2 <- stat_df$mean[stat_df[[group_col]] == case_group]
  
  s1 <- stat_df$sd[stat_df[[group_col]] == control_group]
  s2 <- stat_df$sd[stat_df[[group_col]] == case_group]
  
  # -----------------------------
  # 2. pooled SD
  # -----------------------------
  pooled_sd <- sqrt(
    ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  )
  
  # -----------------------------
  # 3. Cohen's d
  # -----------------------------
  cohens_d <- abs(m2 - m1) / pooled_sd
  
  # -----------------------------
  # 4. power estimation (two-sample t-test)
  # -----------------------------
  power_res <- pwr.t2n.test(
    n1 = n1,
    n2 = n2,
    d = cohens_d,
    sig.level = 0.05,
    alternative = "two.sided"
  )
  
  # -----------------------------
  # 5. output table
  # -----------------------------
  power_table <- data.frame(
    Protein = target,
    n_Control = n1,
    n_Case = n2,
    mean_Control = m1,
    mean_Case = m2,
    sd_Control = s1,
    sd_Case = s2,
    pooled_sd = pooled_sd,
    cohens_d = cohens_d,
    estimated_power = power_res$power
  )
  
  return(power_table)
}

plot_volcano <- function(deg_table,
                         title = "Volcano Plot",
                         file_path = NULL,
                         fc_cutoff = 2,
                         p_cutoff = 0.05,
                         logfc_col = "logFC",
                         padj_col = "adj.P.Val") {
  
  library(dplyr)
  library(ggplot2)
  
  deg_table <- deg_table %>%
    mutate(
      change = case_when(
        .data[[padj_col]] < p_cutoff & .data[[logfc_col]] > log2(fc_cutoff) ~ "UP",
        .data[[padj_col]] < p_cutoff & .data[[logfc_col]] < -log2(fc_cutoff) ~ "DOWN",
        TRUE ~ "NOT"
      )
    ) %>%
    arrange(desc(.data[[logfc_col]]))
  
  # clip extreme values for visualization
  max_fc <- max(abs(deg_table[[logfc_col]]), na.rm = TRUE)
  max_fc <- min(max_fc, 10)
  
  deg_table[[logfc_col]] <- pmax(
    pmin(deg_table[[logfc_col]], max_fc),
    -max_fc
  )
  
  # counts
  stat <- table(deg_table$change)
  
  legend_label <- paste0(names(stat), " (", as.numeric(stat), ")")
  
  p <- ggplot(deg_table,
              aes(x = .data[[logfc_col]],
                  y = -log10(.data[[padj_col]]),
                  color = change)) +
    
    geom_point(alpha = 0.5, size = 1.2) +
    
    scale_color_manual(
      values = c("UP" = "#d44a3d",
                 "DOWN" = "#387eb8",
                 "NOT" = "grey"),
      labels = legend_label
    ) +
    
    geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
               linetype = 4,
               color = "black") +
    
    geom_hline(yintercept = -log10(p_cutoff),
               linetype = 4,
               color = "black") +
    
    theme_minimal() +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10 adj.P.Val") +
    theme(legend.title = element_blank())
  
  # save if path provided
  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 5, height = 4)
  }
  
  return(p)
}

run_enrichment <- function(genes,
                           db = "GO",
                           org_db = org.Hs.eg.db,
                           p_cutoff = 0.05) {
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  
  # -------------------------
  # 1. clean input
  # -------------------------
  genes <- unique(na.omit(genes))
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  # -------------------------
  # 2. SYMBOL → ENTREZID
  # -------------------------
  ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  )
  
  if (is.null(ids) || nrow(ids) == 0) {
    return(NULL)
  }
  
  # -------------------------
  # 3. GO enrichment
  # -------------------------
  if (db == "GO") {
    
    ego <- enrichGO(
      gene = ids$ENTREZID,
      OrgDb = org_db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = p_cutoff,
      readable = TRUE
    )
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      return(NULL)
    }
    
    return(as.data.frame(ego))
  }
  
  # -------------------------
  # 4. KEGG enrichment
  # -------------------------
  if (db == "KEGG") {
    
    ek <- enrichKEGG(
      gene = ids$ENTREZID,
      organism = "hsa",
      pvalueCutoff = p_cutoff
    )
    
    if (is.null(ek) || nrow(as.data.frame(ek)) == 0) {
      return(NULL)
    }
    
    return(as.data.frame(ek))
  }
  
  return(NULL)
}
