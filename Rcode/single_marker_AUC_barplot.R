GetDat <- function(MergeData, cpar, select_genes = c(genes_PAAD, genes_Healthy)){
  if(cpar[2] == "PAAD"){
    MergeData_merg <- MergeData %>% dplyr::select(one_of(c("Group", select_genes)))
  }else{
    cpar2 = strsplit(cpar[2],"/")[[1]]
    MergeData2 <- MergeData %>% dplyr::select(-one_of("Group"))
    rownames(MergeData2) <- gsub("-", "_", rownames(MergeData2))
    MergeData2 <-  inner_join(phen_all[, c("SampleID", "stage")], 
                              MergeData2 %>% rownames_to_column("SampleID"), by = ("SampleID")) %>%
      column_to_rownames("SampleID") 
    
    MergeData_P1 <-  MergeData2  %>% filter(stage %in% cpar2) %>%
      mutate(Group = paste0(cpar[2])) %>%
      dplyr::select(one_of(c("Group", select_genes)))
    
    MergeData_P2 <- MergeData %>% filter(Group %in% cpar[1])
    rownames(MergeData_P2) <- gsub("-", "_", rownames(MergeData_P2))
    MergeData_P2 <-  MergeData_P2 %>%
      dplyr::select(one_of(c("Group", select_genes)))
    
    MergeData_merg <- rbind(MergeData_P1, MergeData_P2)
  }
  return(MergeData_merg)
}

get_raincloud <- function(
  dat_wide,
  dat,
  group = "Group",
  group_names = "all",
  group_colors,
  nrow_num = 1) {
  
  # dat_wide = MergeData
  # dat = MergeData_long
  # group = "Group"
  # group_names = grp_names
  # group_colors = c("#0073C2FF", "#CD534CFF")
  
  dat_cln2 <- dat
  colnames(dat_cln2)[which(colnames(dat_cln2) == group)] <- "Group_new"
  
  if (group_names[1] == "all") {
    tempdata <- dat_cln2
  } else {
    tempdata <- dat_cln2 %>%
      dplyr::filter(Group_new %in% group_names)
  }
  tempdata$Group_new <- factor(tempdata$Group_new, levels = group_names)
  
  plotdata <- tempdata
  
  cmp <- list()
  num <- utils::combn(length(group_names), 2)
  for (i in 1:ncol(num)) {
    cmp[[i]] <- num[, i]
  }
  
  pl <- ggplot(plotdata, aes(x = Group_new, y = log2(value + 1), fill = Group_new)) +
    # ggdist::stat_halfeye(adjust = 0.5, width = 0.3,
    #                      .width = 0, justification = -0.3, point_colour = NA) +
    geom_violin(trim = F, outlier.shape = NA, color = NA) +
    geom_boxplot(width = 0.1, color = "white", outlier.shape = NA, alpha = 0.7) +
    # stat_boxplot(aes(color = Group_new), geom = "errorbar", width = 0.1) +    
    # geom_boxplot(width = 0.1, outlier.shape = NA) + 
    # gghalves::geom_half_point(side = "l", range_scale = 0.4, alpha = 0.5) +
    stat_summary(geom = "crossbar", width = 0.08, fatten = 0, color = "white", 
                 fun.data = function(x){c(y = median(x), ymin = median(x), ymax = median(x))}) +       
    labs(x = "", y = "log2 (TPM)") + 
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    guides(fill = "none", color = "none") + 
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    # ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),
    #                            method = "wilcox.test",
    #                            comparisons = cmp) +
    geom_signif(comparisons = cmp,
                map_signif_level = TRUE,
                test = "wilcox.test",
                tip_length = 0) +
    facet_wrap(.~ Gene, scales = "free_y", nrow = nrow_num) +     
    theme_minimal() +
    theme(axis.title.y = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold", color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())
  
    pl2 <- ggplot(plotdata, aes(x = Group_new, y = log2(value + 1))) +
      geom_violin(
        aes(color = Group_new),
        fill = NA,              
        trim = FALSE,
        linewidth = 0.9,
        scale = "width", 
        adjust = 4) +
      geom_jitter(
        aes(color = Group_new),
        width = 0.15,
        size = 1.5,       
        alpha = 0.65) +
      stat_summary(
        fun = median,
        geom = "point",
        size = 2.2,
        color = "grey") +
      labs(x = "", y = "log2 (TPM+1)") +
      scale_color_manual(values = group_colors) +
      guides(fill = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      geom_signif(
        comparisons = cmp,
        map_signif_level = TRUE,
        test = "wilcox.test",
        tip_length = 0) +
      
      facet_wrap(. ~ Gene, scales = "free_y", nrow = nrow_num) +
      theme_minimal() +
      theme(
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line()
      )
  return(pl2)
}


SingleGeneROC <- function(MergeData, grps){
  # MergeData = dat
  # grps = c("Healthy", "PAAD")
  
  MergeData$Group <- ifelse(MergeData$Group == grps[1], "Ctrl", "Case")
  MergeData$Group <- factor(MergeData$Group, levels = c("Ctrl", "Case"))
  
  auc_dat <- matrix(NA, nrow = length(MergeData)-1, ncol = 3) %>%
    data.frame()
  colnames(auc_dat) <- c("lower", "AUC", "upper")
  genes = colnames(MergeData)[-1]
  rownames(auc_dat) <- genes
  
  for(i in 1:length(genes)){
    
    ## split data
    set.seed(222)
    trainIndex <- caret::createDataPartition(
      MergeData$Group, 
      p = 0.8, 
      list = FALSE, 
      times = 1)
    
    trainData <- MergeData[trainIndex, c(1, i+1)]
    X_train <- subset(trainData, select = colnames(trainData)[2])
    y_train <- subset(trainData, select = colnames(trainData)[1])
    
    testData <- MergeData[-trainIndex, c(1, i+1)]
    X_test <- subset(testData, select = colnames(testData)[2])
    y_test <- subset(testData, select = colnames(testData)[1])
    
    ## model
    # sigmoid <- function(x) {1 / ( 1 + exp(-x) )} # sigmoid
    # pred_prob <- sigmoid(X_train)
    
    lr_fit_optimal <- glm(Group ~ ., data = trainData, family = binomial) # LR
    pred_prob <- predict(lr_fit_optimal, X_test, type="response")

    rocobj_CI <- roc(testData$Group, pred_prob,  ci = TRUE, percent = TRUE)
    roc_CI <- round(as.numeric(rocobj_CI$ci)/100, 3)

    auc_dat[i, 1] = roc_CI[1]
    auc_dat[i, 2] = roc_CI[2]
    auc_dat[i, 3] = roc_CI[3]
  }

  auc_dat <- auc_dat %>% rownames_to_column("gene") %>%
    arrange(desc(AUC))
  
  return(auc_dat)
}

MultiGeneROC <- function(MergeData, dataset_number = "TCGA", grps){
 
  # MergeData = dat
  # dataset_number = "TCGA"
  # grps = c("Healthy", "PAAD")
  
  MergeData$Group <- ifelse(MergeData$Group == grps[2], "Case", "Ctrl")
  MergeData$Group <- factor(MergeData$Group, levels = c("Ctrl", "Case"))
  set.seed(123)
  trainIndex <- caret::createDataPartition(
    MergeData$Group, 
    p = 0.8, 
    list = FALSE, 
    times = 1)
  
  trainData <- MergeData[trainIndex, ]
  X_train <- trainData[, -1]
  y_train <- trainData[, 1]
  
  testData <- MergeData[-trainIndex, ]
  X_test <- testData[, -1]
  y_test <- testData[, 1]
  
  # models
  ## LR
  lr_fit_optimal <- glm(Group ~ ., data = trainData, family = binomial)
  # summary(lr_fit_optimal)
  # cat("模型的参考水平（基线）:", levels(lr_fit_optimal$model$Group)[1], "\n")
  
  ### Confusion matrix
  pred_prob <- predict(lr_fit_optimal, X_test, type="response")
  pred_prob_lr = pred_prob

  ## RF
  set.seed(123)
  # N-repeat K-fold cross-validation
  myControl <- trainControl(
    method = "cv",
    number = 10,
    search = "random",
    classProbs = TRUE,
    verboseIter = FALSE,
    allowParallel = TRUE)
  tuneGrid <- expand.grid(.mtry = c(1:min(10, ncol(X_train)-1)))
  
  set.seed(123)
  tune_fit <- train(
    Group ~.,
    data = trainData,
    method = "rf",
    trControl = myControl,
    tuneGrid = tuneGrid,
    metric = "Accuracy",
    verbose = FALSE) 
  
  optimalVar <- tune_fit$bestTune

  set.seed(123)
  rf_fit_optimal <- randomForest::randomForest(
    Group ~ ., 
    data = trainData, 
    importance = TRUE,
    ntree = 1000,
    mtry = optimalVar$mtry,
    proximity = ifelse(ncol(X_train) <= 20, TRUE, FALSE)
  )
  pred_raw <- predict(rf_fit_optimal, newdata = X_test, type = "response")
  pred_prob <- predict(rf_fit_optimal, newdata = X_test, type = "prob")
  pred_prob_rf = pred_prob[, "Case"]
  
  ## SVM
  myControl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "random",
    summaryFunction = twoClassSummary,  
    classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE)
  
  tuneGrid <- expand.grid(
    sigma= 2^c(-25, -20, -15,-10, -5, 0), 
    C = 2^c(0:5))
  
  set.seed(123)
  tune_fit <- train(
    Group ~.,
    data = trainData,
    method = "svmRadial",
    trControl = myControl,
    tuneGrid = tuneGrid)
  
  optimalVar <- data.frame(tune_fit$results[which.max(tune_fit$results[, 3]), ])
  
  set.seed(123)
  svm_fit_optimal <- e1071::svm(
    Group ~ ., 
    data = trainData, 
    kernel = "radial",
    metric = "ROC",
    gamma = optimalVar$sigma,
    cost = optimalVar$C,
    summaryFunction = twoClassSummary,  
    probability = TRUE)
  
  pred_raw <- predict(svm_fit_optimal, newdata = X_test, type = "response")
  pred_prob_temp <- predict(svm_fit_optimal, newdata = X_test, probability = TRUE)
  pred_prob <- attr(pred_prob_temp, "probabilities")
  pred_prob_svm = pred_prob[,"Case"]
  
  ## merge one plot
  threeML_plot = multi_ROC(
    DataTest = testData,
    PredProb1 = data.frame(pred_prob_lr),
    PredProb2 = data.frame(pred_prob_rf),
    PredProb3 = data.frame(pred_prob_svm),
    tit = dataset_number,
    cols = c("LR" = "#0072B5FF", "RF" = "#BC3C29FF", "SVM" = "#E18727FF")
  )
  
  return(threeML_plot)
}