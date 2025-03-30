library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(pROC)
library(sva)
library(caret)

# Input of data and formatting to tibble
training_data <- as_tibble(as.matrix(readRDS("C:/Users/callu/Desktop/BIOF_520/UROMOL_TaLG.teachingcohort.rds")))
validation_data <- as_tibble(as.matrix(readRDS("C:/Users/callu/Desktop/BIOF_520/knowles_matched_TaLG_final.rds")))

# Remove samples with no recurrence data
training_data <- training_data %>%
  filter(!is.na(Recurrence))
validation_data <- validation_data %>%
  filter(!is.na(Recurrence))

# Select columns that we have data for in both cohorts
white_list_colnames <- intersect(colnames(training_data), colnames(validation_data))

training_whitelist <- training_data %>%
  dplyr::select(white_list_colnames) %>%
  dplyr::select(-Progression, -PFS_time.) 

validation_whitelist <- validation_data %>%
  dplyr::select(white_list_colnames) %>%
  dplyr::select(-Progression, -PFS_time.)

# Identify NAs in training data
NA_counts <- training_whitelist %>% summarise(across(everything(), ~sum(is.na(.))))

# Median imputation for 3 NA age values
training_whitelist$Age[is.na(training_whitelist$Age)] <- median(training_whitelist$Age, na.rm = TRUE)

# Batch effect correct expression data
training_expression <- training_whitelist[, 11:19097]
validation_expression <- validation_whitelist[, 11:19097]
training_expression_t <- t(training_expression)
validation_expression_t <- t(validation_expression)

expression_combined <- cbind(training_expression_t, validation_expression_t)
mode(expression_combined) <- "numeric"

batch <- c(rep(1, ncol(training_expression_t)), rep(2, ncol(validation_expression_t)))
corrected_expression <- ComBat(dat = expression_combined, batch = batch, par.prior = TRUE, prior.plots = FALSE)
corrected_expression_t <- t(corrected_expression)
training_whitelist[, 11:19097] <- corrected_expression_t[1:nrow(training_whitelist), ]
validation_whitelist[, 11:19097] <- corrected_expression_t[(nrow(training_whitelist) + 1):nrow(corrected_expression_t), ]

# Create objects for regression analyses
y <- as.character(training_whitelist$Recurrence)

# Drop unnecessary variables and create matrix for regression analysis
x <- as.data.frame(training_whitelist[, !colnames(training_whitelist) %in% c("RFS_time", "Recurrence",
                                                                         "Tumor.stage", "Tumor.grade")])

numeric_columns <- c(1:2, 7:19093)

x_numeric <- sapply(x[, numeric_columns], as.numeric)
x_numeric <- scale(x_numeric)
x[, numeric_columns] <- x_numeric

validation_x <- as.data.frame(validation_whitelist[, colnames(validation_whitelist) %in% white_list_colnames &
                                          !colnames(validation_whitelist) %in% c("RFS_time", "Recurrence",
                                                                            "Tumor.stage", "Tumor.grade",
                                                                            "Progression","PFS_time.")])

validation_x_numeric <- sapply(validation_x[, numeric_columns], as.numeric)
validation_x_numeric <- scale(validation_x_numeric, center = attr(x_numeric, "scaled:center"), scale = attr(x_numeric, "scaled:scale"))
validation_x[, numeric_columns] <- validation_x_numeric

# Univariate cox regression to select relevant genes
training_expression_scaled <- x[, 7:19093]  
recurrence <- as.numeric(trimws(training_whitelist$Recurrence))
RFS_time <- as.numeric(trimws(training_whitelist$RFS_time))

p_values <- numeric(ncol(training_expression_scaled))

for (i in 1:ncol(training_expression_scaled)) {
  gene_expr <- training_expression_scaled[, i]
  cox_model <- coxph(Surv(RFS_time, recurrence) ~ gene_expr)
  p_values[i] <- summary(cox_model)$coefficients[5]
}

gene_pvalues <- data.frame(gene = colnames(training_expression_scaled), p_value = p_values)

significant_genes <- gene_pvalues[gene_pvalues$p_value < 0.001, ]

# Filter out genes with p value > 0.001 in x
x_clinical_vars <- x[, 1:6]
x_gene_columns <- colnames(x)[7:ncol(x)]
selected_gene_columns <- x_gene_columns[x_gene_columns %in% significant_genes$gene]

x <- cbind(x_clinical_vars, x[, selected_gene_columns])

# Filter out genes with p value > 0.001 in validation_x
validation_x_clinical_vars <- validation_x[, 1:6]
validation_x_gene_columns <- colnames(validation_x)[7:ncol(validation_x)]
selected_gene_columns <- validation_x_gene_columns[validation_x_gene_columns %in% significant_genes$gene]

validation_x <- cbind(validation_x_clinical_vars, validation_x[, selected_gene_columns])

dim(x)
dim(validation_x)

# Dimmensionality reduction with PCA
x_clinical <- x[, 1:6]
x_genes <- x[, 7:ncol(x)]
pca_model <- prcomp(x_genes, center = TRUE, scale. = TRUE)
explained_variance <- summary(pca_model)$importance[2, ]  
num_pcs <- which(cumsum(explained_variance) >= 0.90)[1]
x_pca_scores <- pca_model$x[, 1:num_pcs]
x <- cbind(x_clinical, x_pca_scores)

# Apply PCA to validation data
validation_pca <- predict(pca_model, newdata = validation_x[, 7:ncol(validation_x)])
validation_pca <- validation_pca[, 1:num_pcs]
validation_x <- cbind(validation_x[, 1:6], validation_pca)

# Set seed for reproducible regression analyses
set.seed(154)

#  Set categorical clinical variables as numbers
x$Sex <- ifelse(x$Sex == "M", 1, 0)
x$Concomitant.CIS <- ifelse(x$Concomitant.CIS == "Yes", 1, 0)
x$Class_1 <- as.integer(trimws(x$UROMOL2021.classification) == "Class 1")
x$Class_2a <- as.integer(trimws(x$UROMOL2021.classification) == "Class 2a")
x$Class_2b <- as.integer(trimws(x$UROMOL2021.classification) == "Class 2b")
x$Class_3 <- as.integer(trimws(x$UROMOL2021.classification) == "Class 3")
x <- x %>%
  select(-UROMOL2021.classification, -Class_1)

validation_x$Sex <- ifelse(validation_x$Sex == "M", 1, 0)
validation_x$Concomitant.CIS <- ifelse(validation_x$Concomitant.CIS == "Yes", 1, 0)
validation_x$Class_1 <- as.integer(trimws(validation_x$UROMOL2021.classification) == "Class_1")
validation_x$Class_2a <- as.integer(trimws(validation_x$UROMOL2021.classification) == "Class_2a")
validation_x$Class_2b <- as.integer(trimws(validation_x$UROMOL2021.classification) == "Class_2b")
validation_x$Class_3 <- as.integer(trimws(validation_x$UROMOL2021.classification) == "Class_3")
validation_x <- validation_x %>%
  select(-UROMOL2021.classification, -Class_1)

x <- as.matrix(x)
mode(x) <- "numeric"

validation_x <- as.matrix(validation_x)
mode(validation_x) <- "numeric"

# LASSO with logistic regression
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
plot(cv_lasso)
best_lambda <- cv_lasso$lambda.min
final_model <- glmnet(x, y, family = "binomial", alpha = 1, lambda = 2*best_lambda)
selected_features <- as.matrix(coef(final_model))
selected_features <- selected_features[selected_features != 0, ]
selected_features_df <- as.data.frame(selected_features)
selected_features_df$variables <- rownames(selected_features_df)

training_risk <- as.data.frame(predict(final_model, newx = x,  type = "response"))
training_risk$recurrence <- training_whitelist$Recurrence
training_risk$RFS_time <- training_whitelist$RFS_time
colnames(training_risk) <- c("probability", "recurrence", "RFS_time")

training_roc_curve <- roc(training_risk$recurrence, training_risk$probability)
plot(training_roc_curve, main = "ROC Curve")
auc(training_roc_curve)

optimal_threshold <- coords(training_roc_curve, "best", ret = "threshold")

training_risk <- training_risk %>%
  mutate(risk_group = ifelse(training_risk$probability > median(training_risk$probability), "high_risk", "low_risk"))
table(training_risk$risk_group, str_trim(training_risk$recurrence))

# Validation of model with validation data set
validation_risk <- as.data.frame(predict(final_model, newx = validation_x,  type = "response"))
validation_risk$recurrence <- validation_data$Recurrence
validation_risk$RFS_time <- validation_whitelist$RFS_time
colnames(validation_risk) <- c("probability", "recurrence", "RFS_time")

validation_risk <- validation_risk %>%
  mutate(risk_group = ifelse(validation_risk$probability > as.numeric(optimal_threshold$threshold), "high_risk", "low_risk"))
table(validation_risk$risk_group, str_trim(validation_risk$recurrence))

# ROC curve/AUC performance metric for validation set
validation_roc_curve <- roc(validation_risk$recurrence, validation_risk$probability)
plot(validation_roc_curve, main = "ROC Curve")
auc(validation_roc_curve)

# Kaplan-Meier plot with training data
training_km_fit <- survfit(Surv(as.numeric(RFS_time), as.numeric(recurrence)) ~ risk_group, data = training_risk)

ggsurvplot(training_km_fit, 
           data = training_risk,
           pval = TRUE,
           pval.method = TRUE,              
           risk.table = TRUE,
           ggtheme = theme_bw(),
           title = "Recurrence-free survival of NMIBC patients (Training Set)",
           xlab = "Time (months)", 
           ylab = "Recurrence-free probability",
           legend.title = "Risk Group")

# Kaplan-Meier plot with validation data
validation_km_fit <- survfit(Surv(as.numeric(RFS_time), as.numeric(recurrence)) ~ risk_group, data = validation_risk)

ggsurvplot(validation_km_fit, 
           data = validation_risk,
           pval = TRUE,
           pval.method = TRUE,              
           risk.table = TRUE,
           ggtheme = theme_bw(),
           title = "Recurrence-free survival of NMIBC patients (Validation Set)",
           xlab = "Time (months)", 
           ylab = "Recurrence-free probability",
           legend.title = "Risk Group")

