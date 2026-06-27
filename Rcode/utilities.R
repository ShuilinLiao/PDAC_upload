Define_Change <- function(deg_main, sel_Pcol, sel_FCcol, p_thred = 0.05, fc_thred = 2){
  if(!(sel_Pcol %in% colnames(deg_main))) {
    stop(paste("列", sel_Pcol, "不存在于数据框中"))
  }
  if(!(sel_FCcol %in% colnames(deg_main))) {
    stop(paste("列", sel_FCcol, "不存在于数据框中"))
  }
  if("change" %in% colnames(deg_main)) {
    message("检测到已存在的 'change' 列，将被覆盖更新")
  }
  
  deg_main$change <- "NotSig"
  up_idx <- deg_main[[sel_Pcol]] < p_thred & deg_main[[sel_FCcol]] > log2(fc_thred)
  down_idx <- deg_main[[sel_Pcol]] < p_thred & deg_main[[sel_FCcol]] < -log2(fc_thred)
  
  deg_main$change[up_idx] <- "Up"
  deg_main$change[down_idx] <- "Down"
  
  return(deg_main)
}

AUROC <- function(
  DataTest, 
  PredProb = pred_prob, 
  label = group_names[1], 
  DataProf, 
  tit = "") {
  
  # DataTest = y_test
  # PredProb = pred_prob
  # label = group_names[1]
  # DataProf = 41
  
  # ROC object
  rocobj <- roc(DataTest$Group, PredProb[, 1])
  
  
  # AUROC data
  roc <- data.frame(tpr = rocobj$sensitivities,
                    fpr = 1 - rocobj$specificities)
  
  # AUC 95% CI
  rocobj_CI <- roc(DataTest$Group, PredProb[, 1], 
                   ci = TRUE, percent = TRUE)
  roc_CI <- round(as.numeric(rocobj_CI$ci)/100, 3)
  # roc_CI_lab <- paste0(label, 
  #                      " (", "AUC=", roc_CI[2], 
  #                      ", 95%CI ", roc_CI[1], "-", roc_CI[3], 
  #                      ")")
  
  roc_CI_lab <- paste0("AUC=", roc_CI[2], 
                       "\n95% CI: ", roc_CI[1], "-", roc_CI[3])  
  # ROC dataframe
  rocbj_df <- data.frame(threshold = round(rocobj$thresholds, 3),
                         sensitivities = round(rocobj$sensitivities, 3),
                         specificities = round(rocobj$specificities, 3),
                         value = rocobj$sensitivities + 
                           rocobj$specificities)
  max_value_row <- which(max(rocbj_df$value) == rocbj_df$value)
  threshold <- rocbj_df$threshold[max_value_row]
  
  # plot
  pl <- ggplot(data = roc, aes(x = fpr, y = tpr)) +
    geom_path(color = "red", size = 1) +
    geom_abline(intercept = 0, slope = 1, 
                color = "grey", linewidth = 1, linetype = 2) +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate",
         title = paste0("AUROC (", DataProf, " Features)")) +
    annotate("text", 
             x = 1 - rocbj_df$specificities[max_value_row] + 0.15, 
             y = rocbj_df$sensitivities[max_value_row] - 0.05, 
             label = paste0(threshold, " (", 
                            rocbj_df$specificities[max_value_row], ",",
                            rocbj_df$sensitivities[max_value_row], ")"),
             size=5) +
    annotate("point", 
             x = 1 - rocbj_df$specificities[max_value_row], 
             y = rocbj_df$sensitivities[max_value_row], 
             color = "black", size = 2) +    
    annotate("text", 
             x = .75, y = .25, 
             label = roc_CI_lab,
             size = 5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    ggtitle(tit) +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.title = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 11, color = "black", face = "bold"), 
          axis.text = element_text(size= 10, color = "black"),
          axis.ticks.length = unit(0.4, "lines"),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())
  
  res <- list(rocobj = rocobj,
              roc_CI = roc_CI_lab,
              roc_pl = pl)
  
  return(res)
}

multi_ROC <- function(DataTest, PredProb1, PredProb2, PredProb3, tit, cols){
  
  # DataTest = testData
  # PredProb1 = data.frame(pred_prob = pred_prob_lr)
  # PredProb2 = pred_prob_rf
  # PredProb3 = pred_prob_svm
  # tit = ""
  # cols = c("LR" = "#0072B5FF", "RF" = "#BC3C29FF", "SVM" = "#E18727FF")
  
  rocobj1 <- roc(DataTest$Group, PredProb1[, 1])
  rocobj2 <- roc(DataTest$Group, PredProb2[, 1])
  rocobj3 <- roc(DataTest$Group, PredProb3[, 1])
  
  roc1 <- data.frame(tpr = rocobj1$sensitivities, fpr = 1 - rocobj1$specificities, model = "LR")
  roc2 <- data.frame(tpr = rocobj2$sensitivities, fpr = 1 - rocobj2$specificities, model = "RF")
  roc3 <- data.frame(tpr = rocobj3$sensitivities, fpr = 1 - rocobj3$specificities, model = "SVM")
  
  rocobj_CI <- roc(DataTest$Group, PredProb1[, 1],  ci = TRUE, percent = TRUE)
  roc_CI <- round(as.numeric(rocobj_CI$ci)/100, 3)
  roc_CI_lab1 <- paste0("LR, AUC = ", roc_CI[2], " (95% CI: ", roc_CI[1], "-", roc_CI[3], ")") 
  
  rocobj_CI <- roc(DataTest$Group, PredProb2[, 1],  ci = TRUE, percent = TRUE)
  roc_CI <- round(as.numeric(rocobj_CI$ci)/100, 3)
  roc_CI_lab2 <- paste0("RF, AUC = ", roc_CI[2], " (95% CI: ", roc_CI[1], "-", roc_CI[3], ")") 
  
  rocobj_CI <- roc(DataTest$Group, PredProb3[, 1],  ci = TRUE, percent = TRUE)
  roc_CI <- round(as.numeric(rocobj_CI$ci)/100, 3)
  roc_CI_lab3 <- paste0("SVM, AUC = ", roc_CI[2], " (95% CI: ", roc_CI[1], "-", roc_CI[3], ")") 
  
  roc_all = rbind(roc1, roc2, roc3)
  
  # roc_all <- data.frame(roc1_tpr = roc1$tpr, roc1_fpr = roc1$fpr, 
  #                       roc2_tpr = roc2$tpr, roc2_fpr = roc2$fpr,
  #                       roc3_tpr = roc3$tpr, roc3_fpr = roc3$fpr)
  
  
  pl <- ggplot(data = roc_all, aes(x = fpr, y = tpr, color = model)) +
    geom_path(size = 1) +  # Draw each ROC curve with the specified color
    geom_abline(intercept = 0, slope = 1,
                color = "grey", linewidth = 1, linetype = 2) +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate",
         title = paste0("AUROC (Comparison of Models)")) +
    annotate("text", x = 0.65, y = 0.25, size = 4.5, label = roc_CI_lab1, color = cols[1]) +
    annotate("text", x = 0.65, y = 0.15, size = 4.5, label = roc_CI_lab2, color = cols[2]) +
    annotate("text", x = 0.65, y = 0.05, size = 4.5, label = roc_CI_lab3, color = cols[3]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    ggtitle(tit) +
    theme(
          plot.title = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 11, color = "black", face = "bold"), 
          axis.text = element_text(size= 10, color = "black"),
          
          # axis.ticks.length = unit(0.4, "lines"),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.position="none") +
    scale_color_manual(values = cols)
  
  return(pl)
}
