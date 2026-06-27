# minimal_mul_ML.R
# Helper functions for minimal multigene classifier optimization,
# compartment-specific ERS model fitting, performance evaluation, and ROC plotting.

load_required_packages <- function() {
  pkgs <- c("tidyverse", "caret", "pROC", "randomForest", "e1071")
  invisible(lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package not installed: ", pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

create_project_dirs <- function(...) {
  dirs <- c(...)
  invisible(lapply(dirs, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)))
}

build_external_paths <- function(external_dir,
                                 cohort_names = c("GSE28735", "GSE62452", "GSE71729")) {
  stats::setNames(file.path(external_dir, paste0(cohort_names, ".RData")), cohort_names)
}

save_csv <- function(x, out_dir, filename, row.names = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file.path(out_dir, filename), row.names = row.names)
}

load_discovery_rdata <- function(rdata_path) {
  if (!file.exists(rdata_path)) stop("File not found: ", rdata_path)
  env <- new.env()
  load(rdata_path, envir = env)
  if (!exists("profile", envir = env)) stop("Object 'profile' was not found in ", rdata_path)
  if (!exists("metadata", envir = env)) stop("Object 'metadata' was not found in ", rdata_path)

  profile_obj <- get("profile", envir = env)
  metadata_obj <- get("metadata", envir = env)

  if (!"Group" %in% colnames(metadata_obj)) {
    if ("Compvar" %in% colnames(metadata_obj)) {
      metadata_obj <- metadata_obj %>% dplyr::rename(Group = Compvar)
    } else {
      stop("Neither 'Group' nor 'Compvar' was found in metadata from ", rdata_path)
    }
  }

  list(profile = profile_obj, metadata = metadata_obj)
}

prepare_model_data <- function(profile, metadata, genes,
                               sample_col = "SampleID",
                               group_col = "Group",
                               negative = "Healthy",
                               positive = "PAAD") {
  available_genes <- intersect(genes, rownames(profile))
  missing_genes <- setdiff(genes, rownames(profile))
  if (length(missing_genes) > 0) {
    warning("Missing genes removed: ", paste(missing_genes, collapse = ", "))
  }
  if (length(available_genes) == 0) stop("No model genes were found in profile.")

  profile_sub <- profile[available_genes, , drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(sample_col)

  dat <- metadata %>%
    dplyr::select(dplyr::all_of(c(sample_col, group_col))) %>%
    dplyr::inner_join(profile_sub, by = sample_col) %>%
    tibble::column_to_rownames(sample_col) %>%
    dplyr::mutate(
      Group = .data[[group_col]],
      Group = dplyr::recode(Group,
                            "Normal" = negative,
                            "Control" = negative,
                            "HC" = negative,
                            "Tumor" = positive,
                            "Cancer" = positive,
                            "PDAC" = positive,
                            .default = as.character(Group)),
      Group = factor(Group, levels = c(negative, positive))
    ) %>%
    dplyr::select(Group, dplyr::all_of(available_genes)) %>%
    tidyr::drop_na(Group)

  dat
}

split_discovery_data <- function(dat, p = 0.8, seed = 123) {
  set.seed(seed)
  train_index <- caret::createDataPartition(dat$Group, p = p, list = FALSE)
  list(
    train = dat[train_index, , drop = FALSE],
    test  = dat[-train_index, , drop = FALSE]
  )
}

get_default_candidate_sets <- function() {
  list(
    "SERF2" = c("SERF2"),
    "LCN2" = c("LCN2"),
    "LRRC42" = c("LRRC42"),
    "FLOT1" = c("FLOT1"),
    "LCN2_LRRC42" = c("LCN2", "LRRC42"),
    "LCN2_FLOT1" = c("LCN2", "FLOT1"),
    "LRRC42_FLOT1" = c("LRRC42", "FLOT1"),
    "SERF2_LCN2" = c("SERF2", "LCN2"),
    "LCN2_LRRC42_FLOT1" = c("LCN2", "LRRC42", "FLOT1"),
    "SERF2_LCN2_LRRC42" = c("SERF2", "LCN2", "LRRC42"),
    "SERF2_LCN2_FLOT1" = c("SERF2", "LCN2", "FLOT1"),
    "SERF2_LRRC42_FLOT1" = c("SERF2", "LRRC42", "FLOT1")
  )
}

filter_candidate_sets <- function(candidate_sets, dat) {
  candidate_sets[sapply(candidate_sets, function(x) all(x %in% colnames(dat)))]
}

get_svm_prob <- function(model, newdata, positive = "PAAD") {
  pred <- predict(model, newdata = newdata, probability = TRUE)
  prob_mat <- attr(pred, "probabilities")
  if (positive %in% colnames(prob_mat)) {
    prob_mat[, positive]
  } else {
    stop("Positive class not found in SVM probability matrix.")
  }
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

fit_ers_models <- function(trainData, genes,
                           positive = "PAAD", seed = 123,
                           tune_svm = TRUE) {
  trainData <- trainData[, c("Group", genes), drop = FALSE]

  lr_fit <- stats::glm(Group ~ ., data = trainData, family = stats::binomial)

  set.seed(seed)
  rf_fit <- randomForest::randomForest(
    Group ~ ., data = trainData, ntree = 500, importance = TRUE, proximity = TRUE
  )

  if (tune_svm) {
    set.seed(seed)
    svm_ctrl <- caret::trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 3,
      classProbs = TRUE,
      summaryFunction = caret::twoClassSummary,
      savePredictions = "final"
    )

    svm_grid <- expand.grid(
      sigma = 2^c(-25, -20, -15, -10, -5, 0),
      C = 2^c(0:5)
    )

    svm_tune <- caret::train(
      Group ~ ., data = trainData,
      method = "svmRadial",
      trControl = svm_ctrl,
      tuneGrid = svm_grid,
      metric = "ROC"
    )

    svm_fit <- e1071::svm(
      Group ~ ., data = trainData,
      kernel = "radial",
      gamma = svm_tune$bestTune$sigma,
      cost = svm_tune$bestTune$C,
      probability = TRUE,
      scale = TRUE
    )
  } else {
    svm_tune <- NULL
    set.seed(seed)
    svm_fit <- e1071::svm(Group ~ ., data = trainData, kernel = "radial", probability = TRUE, scale = TRUE)
  }

  list(
    models = list(LR = lr_fit, RF = rf_fit, SVM = svm_fit),
    svm_tune = svm_tune
  )
}

fit_final_ers_models <- function(trainData, final_ers_genes,
                                 positive = "PAAD", seed = 123,
                                 tune_svm = TRUE) {
  fit_ers_models(trainData = trainData, genes = final_ers_genes,
                 positive = positive, seed = seed, tune_svm = tune_svm)
}

predict_model_prob <- function(model, newdata, model_name, genes,
                               positive = "PAAD") {
  nd <- newdata[, genes, drop = FALSE]
  if (model_name == "LR") {
    prob <- stats::predict(model, newdata = nd, type = "response")
  } else if (model_name == "RF") {
    prob <- stats::predict(model, newdata = nd, type = "prob")[, positive]
  } else if (model_name == "SVM") {
    prob <- get_svm_prob(model, nd, positive = positive)
  } else {
    stop("Unknown model_name: ", model_name)
  }
  clip_prob(prob)
}

evaluate_gene_set <- function(gene_set_name, genes, trainData_all, testData_all,
                              positive = "PAAD", negative = "Healthy", seed = 123) {
  trainData <- trainData_all[, c("Group", genes), drop = FALSE]
  testData  <- testData_all[, c("Group", genes), drop = FALSE]
  y_test <- factor(testData$Group, levels = c(negative, positive))

  lr_fit <- stats::glm(Group ~ ., data = trainData, family = stats::binomial)
  lr_prob <- stats::predict(lr_fit, newdata = testData[, genes, drop = FALSE], type = "response")
  lr_auc <- as.numeric(pROC::auc(pROC::roc(y_test, lr_prob, levels = c(negative, positive), direction = "<", quiet = TRUE)))

  set.seed(seed)
  rf_fit <- randomForest::randomForest(Group ~ ., data = trainData, ntree = 500, importance = TRUE)
  rf_prob <- stats::predict(rf_fit, newdata = testData[, genes, drop = FALSE], type = "prob")[, positive]
  rf_auc <- as.numeric(pROC::auc(pROC::roc(y_test, rf_prob, levels = c(negative, positive), direction = "<", quiet = TRUE)))

  set.seed(seed)
  svm_fit <- e1071::svm(Group ~ ., data = trainData, kernel = "radial", probability = TRUE, scale = TRUE)
  svm_prob <- get_svm_prob(svm_fit, testData[, genes, drop = FALSE], positive = positive)
  svm_auc <- as.numeric(pROC::auc(pROC::roc(y_test, svm_prob, levels = c(negative, positive), direction = "<", quiet = TRUE)))

  tibble::tibble(
    Model = gene_set_name,
    Genes = paste(genes, collapse = " + "),
    N_genes = length(genes),
    LR_AUC = lr_auc,
    RF_AUC = rf_auc,
    SVM_AUC = svm_auc,
    Mean_AUC = mean(c(lr_auc, rf_auc, svm_auc), na.rm = TRUE)
  )
}

run_minimal_optimization <- function(candidate_sets, trainData_all, testData_all,
                                     positive = "PAAD", negative = "Healthy", seed = 123) {
  purrr::imap_dfr(
    candidate_sets,
    ~ evaluate_gene_set(.y, .x, trainData_all, testData_all,
                        positive = positive, negative = negative, seed = seed)
  ) %>%
    dplyr::arrange(dplyr::desc(Mean_AUC))
}

plot_optimization_results <- function(optimization_res, title = "Minimal multigene classifier optimization") {
  plot_df <- optimization_res %>%
    tidyr::pivot_longer(
      cols = c(LR_AUC, RF_AUC, SVM_AUC),
      names_to = "Classifier",
      values_to = "AUC"
    ) %>%
    dplyr::mutate(
      Classifier = stringr::str_replace(Classifier, "_AUC", ""),
      Model = factor(Model, levels = rev(optimization_res$Model))
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = AUC, y = Model, fill = Classifier)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
    ggplot2::geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey50") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "AUC", y = NULL) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
}

load_external_cohort <- function(rdata_path, genes,
                                 negative = "Healthy",
                                 positive = "PAAD") {
  if (!file.exists(rdata_path)) {
    warning("External cohort file not found: ", rdata_path)
    return(NULL)
  }

  raw <- load_discovery_rdata(rdata_path)
  prepare_model_data(
    profile = raw$profile,
    metadata = raw$metadata,
    genes = genes,
    negative = negative,
    positive = positive
  )
}

load_external_cohorts <- function(external_paths, genes,
                                  negative = "Healthy", positive = "PAAD") {
  external_cohorts <- purrr::imap(
    external_paths,
    ~ load_external_cohort(.x, genes = genes, negative = negative, positive = positive)
  )
  external_cohorts[!sapply(external_cohorts, is.null)]
}

get_training_threshold <- function(model, trainData, model_name, genes,
                                   positive = "PAAD", negative = "Healthy") {
  prob <- predict_model_prob(model, trainData, model_name, genes = genes, positive = positive)
  truth <- factor(trainData$Group, levels = c(negative, positive))
  roc_obj <- pROC::roc(truth, prob, levels = c(negative, positive), direction = "<", quiet = TRUE)
  as.numeric(pROC::coords(
    roc_obj, x = "best", best.method = "youden",
    ret = "threshold", transpose = FALSE
  )$threshold[1])
}

binom_ci <- function(x, n) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  as.numeric(stats::binom.test(x, n)$conf.int)
}

calibration_metrics <- function(truth, prob, positive = "PAAD") {
  y <- as.integer(truth == positive)
  prob <- clip_prob(prob)
  brier <- mean((y - prob)^2, na.rm = TRUE)

  cal <- tryCatch({
    fit <- stats::glm(y ~ qlogis(prob), family = stats::binomial)
    c(Calibration_intercept = unname(stats::coef(fit)[1]),
      Calibration_slope = unname(stats::coef(fit)[2]))
  }, error = function(e) {
    c(Calibration_intercept = NA_real_, Calibration_slope = NA_real_)
  })

  hl_p <- if (requireNamespace("ResourceSelection", quietly = TRUE)) {
    tryCatch(ResourceSelection::hoslem.test(y, prob, g = 10)$p.value,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }

  c(Brier_score = brier, cal, HL_p = hl_p)
}

evaluate_classifier <- function(model, newdata, model_name, cohort, threshold, genes,
                                positive = "PAAD", negative = "Healthy",
                                model_domain = NA_character_) {
  dat <- newdata[, c("Group", genes), drop = FALSE] %>% tidyr::drop_na()
  truth <- factor(dat$Group, levels = c(negative, positive))
  prob <- predict_model_prob(model, dat, model_name, genes = genes, positive = positive)

  roc_obj <- pROC::roc(truth, prob, levels = c(negative, positive), direction = "<", quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95))

  pred <- factor(ifelse(prob >= threshold, positive, negative), levels = c(negative, positive))
  tab <- table(Predicted = pred, Reference = truth)

  TP <- ifelse(!is.na(tab[positive, positive]), tab[positive, positive], 0)
  TN <- ifelse(!is.na(tab[negative, negative]), tab[negative, negative], 0)
  FP <- ifelse(!is.na(tab[positive, negative]), tab[positive, negative], 0)
  FN <- ifelse(!is.na(tab[negative, positive]), tab[negative, positive], 0)

  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  sens_ci <- binom_ci(TP, TP + FN)
  spec_ci <- binom_ci(TN, TN + FP)
  cal <- calibration_metrics(truth, prob, positive = positive)

  tibble::tibble(
    Model_domain = model_domain,
    Cohort = cohort,
    Classifier = model_name,
    N = length(truth),
    N_positive = sum(truth == positive),
    N_negative = sum(truth == negative),
    Threshold_from_training = threshold,
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUC_95CI_lower = auc_ci[1],
    AUC_95CI_upper = auc_ci[3],
    Sensitivity = sensitivity,
    Sensitivity_95CI_lower = sens_ci[1],
    Sensitivity_95CI_upper = sens_ci[2],
    Specificity = specificity,
    Specificity_95CI_lower = spec_ci[1],
    Specificity_95CI_upper = spec_ci[2],
    Brier_score = unname(cal["Brier_score"]),
    Calibration_intercept = unname(cal["Calibration_intercept"]),
    Calibration_slope = unname(cal["Calibration_slope"]),
    HL_p = unname(cal["HL_p"])
  )
}

make_performance_table <- function(model_list, all_cohorts, trainData, genes,
                                   positive = "PAAD", negative = "Healthy",
                                   model_domain = NA_character_) {
  thresholds <- purrr::imap_dbl(
    model_list,
    ~ get_training_threshold(.x, trainData, model_name = .y,
                             genes = genes, positive = positive, negative = negative)
  )

  performance_table <- purrr::imap_dfr(model_list, function(model, model_name) {
    purrr::imap_dfr(all_cohorts, function(dat, cohort_name) {
      evaluate_classifier(
        model = model,
        newdata = dat,
        model_name = model_name,
        cohort = cohort_name,
        threshold = thresholds[[model_name]],
        genes = genes,
        positive = positive,
        negative = negative,
        model_domain = model_domain
      )
    })
  }) %>%
    dplyr::mutate(
      AUC_95CI = sprintf("%.3f (%.3f-%.3f)", AUC, AUC_95CI_lower, AUC_95CI_upper),
      Sensitivity_95CI = sprintf("%.3f (%.3f-%.3f)", Sensitivity, Sensitivity_95CI_lower, Sensitivity_95CI_upper),
      Specificity_95CI = sprintf("%.3f (%.3f-%.3f)", Specificity, Specificity_95CI_lower, Specificity_95CI_upper)
    )

  manuscript_table <- performance_table %>%
    dplyr::select(
      Model_domain, Cohort, Classifier, N, N_positive, N_negative,
      AUC_95CI, Sensitivity_95CI, Specificity_95CI,
      Brier_score, Calibration_intercept, Calibration_slope, HL_p
    )

  list(full = performance_table, manuscript = manuscript_table, thresholds = thresholds)
}

make_three_classifier_roc <- function(dat, cohort_name, model_list, genes,
                                      positive = "PAAD", negative = "Healthy",
                                      multi_roc_fun = NULL,
                                      plot_title_prefix = "final ERS",
                                      cols = c("LR" = "#b2182b", "RF" = "#006837", "SVM" = "#4575b4")) {
  dat <- dat[, c("Group", genes), drop = FALSE] %>% tidyr::drop_na()
  dat$Group <- factor(dat$Group, levels = c(negative, positive))

  lr_prob  <- predict_model_prob(model_list$LR,  dat, "LR", genes = genes, positive = positive)
  rf_prob  <- predict_model_prob(model_list$RF,  dat, "RF", genes = genes, positive = positive)
  svm_prob <- predict_model_prob(model_list$SVM, dat, "SVM", genes = genes, positive = positive)

  pred_prob_lr <- data.frame(pred_prob = lr_prob)
  pred_prob_rf <- data.frame(Healthy = 1 - rf_prob, PAAD = rf_prob)
  pred_prob_svm <- data.frame(Healthy = 1 - svm_prob, PAAD = svm_prob)

  if (is.null(multi_roc_fun)) {
    multi_roc_fun <- get("multi_ROC", envir = .GlobalEnv)
  }

  multi_roc_fun(
    DataTest = dat,
    PredProb1 = pred_prob_lr,
    PredProb2 = pred_prob_rf,
    PredProb3 = pred_prob_svm,
    tit = paste0(cohort_name, " ", plot_title_prefix),
    cols = cols
  )
}

save_roc_plots <- function(roc_cohorts, model_list, genes, out_dir,
                           prefix = "ROC_final_ERS",
                           positive = "PAAD", negative = "Healthy",
                           plot_title_prefix = "final ERS") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  roc_plot_list <- purrr::imap(
    roc_cohorts,
    ~ make_three_classifier_roc(.x, .y, model_list = model_list, genes = genes,
                                positive = positive, negative = negative,
                                plot_title_prefix = plot_title_prefix)
  )

  purrr::iwalk(roc_plot_list, function(pl, nm) {
    out_name <- gsub("[^A-Za-z0-9_]+", "_", nm)
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(prefix, "_", out_name, ".pdf")),
      plot = pl, width = 6, height = 5
    )
  })

  grDevices::pdf(file.path(out_dir, paste0(prefix, "_all_cohorts.pdf")), width = 6, height = 5)
  for (nm in names(roc_plot_list)) print(roc_plot_list[[nm]])
  grDevices::dev.off()

  roc_plot_list
}
