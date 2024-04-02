#### correlation testing ####
## correlation test between the selected CpGs and to age as well

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)
library(tibble)
library(dplyr)
library(ggplot2)

## loading data 
# Data should be in environment from previous script

#### functions ####

cor.test.age <- function(methyl_values, age, SMR = "not_defined", species = "undefined", method = "pearson") {
  correlation_results <- list()
  print(paste0("Running correlation test against age with ", method, " method. Results are stored in tibble."))
  # Loop through each methylation site
  for (i in 1:ncol(methyl_values)) {
    site_name <- colnames(methyl_values)[i]
    # Perform correlation test with age
    test_result <- cor.test(methyl_values[,i], age, method = method) # Use "spearman" or "kendall" if more appropriate
    
    # Store the results
    correlation_results[[site_name]] <- list(
      correlation_coefficient = test_result$estimate,
      p_value = test_result$p.value
    )
  }
  
  # Optionally, convert the results list to a more convenient format like a dataframe
  correlation_summary <- tibble(
    Site = names(correlation_results),
    Correlation = sapply(correlation_results, function(x) x$correlation_coefficient),
    P_value = sapply(correlation_results, function(x) x$p_value),
    SMR = SMR,
    species = species
  )
  return(correlation_summary)
}

cor.test.age.filter <- function(input, p_value = 0.05) {
  significant_vector <- as.vector(ifelse(input$P_value <= p_value, TRUE, FALSE))
  input$significant <- significant_vector
  return(input)
}


#### Correlation tests ####

AC_cor_age_pearson <- cor.test.age(AC_meth_values, AC_age, AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], species = "AC")
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

AS_cor_age_pearson <- cor.test.age(AS_meth_values, AS_age, AS_methyl_sites$SMR, species = "AS")
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age_pearson, 0.05)

EH_cor_age_pearson <- cor.test.age(EH_meth_values, EH_age, EH_methyl_sites$SMR, species = "EH")
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age_pearson, 0.05)

ZF_cor_age_pearson <- cor.test.age(ZF_meth_values, ZF_age, ZF_methyl_sites$SMR, species = "ZF")
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age_pearson, 0.05)

cor_all <- rbind(AC_cor_age_filtered_pearson,
                 AS_cor_age_filtered_pearson,
                 EH_cor_age_filtered_pearson, 
                 ZF_cor_age_filtered_pearson)

## selecting CpGs

ncol(AC_meth_values[AC_cor_age_filtered_pearson$significant] == TRUE)
ncol(AS_meth_values[AS_cor_age_filtered_pearson$significant] == TRUE)
ncol(EH_meth_values[EH_cor_age_filtered_pearson$significant] == TRUE)
ncol(ZF_meth_values[ZF_cor_age_filtered_pearson$significant] == TRUE)

#### function to choose only the highest correlating CpGs per SMR

select.max.cor <- function(cor_tibble, filter_significant = FALSE) {
  filtered_data <- cor_tibble
  if(filter_significant == TRUE) {filtered_data <- filter(filtered_data, significant)}
  
  filtered_data <- filtered_data %>% 
    group_by(SMR) %>%
    # Add a temporary column for the absolute correlation values
    mutate(abs_correlation = abs(Correlation)) %>%
    # For each group, filter the row with the max absolute correlation
    filter(abs_correlation == max(abs_correlation)) %>%
    # Remove the temporary column
    select(-abs_correlation) %>%
    # Optionally, ensure only one row per group if there are ties
    slice(1)
  return(filtered_data)
}

## selecting only positively correlating samples
AC_pos_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation > 0,])
AS_pos_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation > 0,])
EH_pos_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation > 0,])
ZF_pos_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation > 0,])

all_pos_cor_CpG <- rbind(AC_pos_cor_CpGs, AS_pos_cor_CpGs, EH_pos_cor_CpGs, ZF_pos_cor_CpGs)
all_pos_cor_CpG_common  <- all_pos_cor_CpG %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

## selecting only negatively correlating samples
AC_neg_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation < 0,])
AS_neg_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation < 0,])
EH_neg_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation < 0,])
ZF_neg_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation < 0,])

all_neg_cor_CpG <- rbind(AC_neg_cor_CpGs, AS_neg_cor_CpGs, EH_neg_cor_CpGs, ZF_neg_cor_CpGs)
all_neg_cor_CpG_common  <- all_neg_cor_CpG %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

## selecting a mixture of both
temp_index_vec <- (all_pos_cor_CpG_common$SMR %in% all_neg_cor_CpG_common$SMR) == FALSE
all_mix_cor_CpG_common  <- rbind(all_neg_cor_CpG_common, all_pos_cor_CpG_common[temp_index_vec,]) %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

## selecting significant ones
AC_sig_CpGs <- select.max.cor(AC_cor_age_filtered_pearson, TRUE)
AS_sig_CpGs <- select.max.cor(AS_cor_age_filtered_pearson, TRUE)
EH_sig_CpGs <- select.max.cor(EH_cor_age_filtered_pearson, TRUE)
ZF_sig_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson, TRUE)

all_sig_CpGs <- rbind(AC_sig_CpGs, AS_sig_CpGs, EH_sig_CpGs, ZF_sig_CpGs)

all_sig_CpGs_common <- all_sig_CpGs %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup


#### model creation ####
library(glmnet)
AC_glm <- cv.glmnet(as.matrix(AC_meth_values), AC_age, alpha = 0.5)
coef(AC_glm, s=0.01)
AC_glm_prediction <- predict(AC_glm, newx = as.matrix(AC_meth_values))
plot(AC_glm)

#### ----
### selecting only pos cor ones 
#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_pos_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_pos_cor_CpG_common$SMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_pos_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_pos_cor_CpG_common$SMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_pos_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_pos_cor_CpG_common$SMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values) %in% all_pos_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_pos_cor_CpG_common$SMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$species <- "ZF"

#### ----
### selecting only neg cor ones 
#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_neg_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_neg_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_neg_cor_CpG_common$SMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_neg_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_neg_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_neg_cor_CpG_common$SMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_neg_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_neg_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_neg_cor_CpG_common$SMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values) %in% all_neg_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_neg_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_neg_cor_CpG_common$SMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$species <- "ZF"

#### ----
# selecting a mixture of them 

#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_mix_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_mix_cor_CpG_common$SMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_mix_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_mix_cor_CpG_common$SMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_mix_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_mix_cor_CpG_common$SMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values) %in% all_mix_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_mix_cor_CpG_common$SMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$species <- "ZF"

#### ----
#combining methylation values
all_meth_values_selected <- rbind(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected)

all_meth_values_selected_train <- all_meth_values_selected[-seq(1, nrow(all_meth_values_selected), 4),]
all_meth_values_selected_test <- all_meth_values_selected[seq(1, nrow(all_meth_values_selected), 4),]

ZF_meth_values_selected_train <- ZF_meth_values_selected[-seq(1, nrow(ZF_meth_values_selected), 4),]
ZF_meth_values_selected_test <- ZF_meth_values_selected[seq(1, nrow(ZF_meth_values_selected), 4),]

save(all_meth_values_selected, file = "/workspace/cfngle/results-data/06_model_creation/HS_all_meth_values_selected.RData")
# for all
X <- all_meth_values_selected_train %>% 
  select(-rel_age, -species) %>% 
  as.matrix()

Y <- all_meth_values_selected_train[,"rel_age"]

X_test <- all_meth_values_selected_test %>% 
  select(-rel_age, -species) %>% 
  as.matrix()
Y_test <- all_meth_values_selected_test[,"rel_age"]


#### Testing model CV GLM ####
set.seed(123)

# change age here!
Y <- all_meth_values_selected_train[,"rel_age"]
Y_test <- all_meth_values_selected_test[,"rel_age"]
#changing value
Y <- -log(-log(Y))
Y_test <- -log(-log(Y_test))

set.seed(123)
GLM_test <- cv.glmnet(X, Y, alpha = 0.5)
# plot(glm_test)

coef(GLM_test, s=GLM_test$lambda.min)

### running model on testing data
### plotting
GLM_predictions_test <- predict(GLM_test, X_test, s=GLM_test$lambda.min)
GLM_predictions_cor_test <- cor(GLM_predictions_test, Y_test, method = "pearson") %>% 
  round(4)
GLM_mse_test <- mean((GLM_predictions_test - Y_test)^2) %>% 
  round(4) 
GLM_mae_test <- mean(abs(GLM_predictions_test - Y_test)) %>% 
  round(4)

# for normal Y
GLM_result_df <- data.frame(predictions = GLM_predictions_test,
                            rel_age = Y_test,
                            species = all_meth_values_selected_test$species)
# for transformed Y!

GLM_result_df <- data.frame(predictions = exp(-exp(-GLM_predictions_test)),
                            rel_age = exp(-exp(-Y_test)),
                            species = all_meth_values_selected_test$species)
# for transformed Y!
GLM_mse_test <- mean((exp(-exp(-GLM_predictions_test)) - exp(-exp(-Y_test)))^2) %>% 
  round(4) 
GLM_mae_test <- mean(abs(exp(-exp(-GLM_predictions_test)) - exp(-exp(-Y_test)))) %>% 
  round(4)


colnames(GLM_result_df) <- c("age_predicted", "age", "species")

GLM_plot_test <- ggplot(GLM_result_df, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction GLM (test data set)",
       subtitle = paste0("R=", GLM_predictions_cor_test,
                         " MSE=", GLM_mse_test,
                         " MAE=", GLM_mae_test,
                         " N=", nrow(X_test))) +
  theme_minimal()

### running model on training data
GLM_predictions_train <- predict(GLM_test, newx= X, s=GLM_test$lambda.min)
GLM_predictions_cor_train <- cor(GLM_predictions_train, Y, method = "pearson") %>% 
  round(4)
GLM_mse_train <- mean((GLM_predictions_train - Y)^2) %>% 
  round(4)
GLM_mae_train <- mean(abs(GLM_predictions_train - Y)) %>% 
  round(4)

# for normal Y
GLM_result_df_train <- data.frame(predictions = GLM_predictions_train,
                                  rel_age = Y,
                                  species = all_meth_values_selected_train$species)

# for transformed Y!
GLM_result_df_train <- data.frame(predictions = exp(-exp(-GLM_predictions_train)),
                                  rel_age = exp(-exp(-Y)),
                                  species = all_meth_values_selected_train$species)
# for transformed Y!
GLM_mse_train <- mean((exp(-exp(-GLM_predictions_train)) - exp(-exp(-Y)))^2) %>% 
  round(4) 
GLM_mae_train <- mean(abs(exp(-exp(-GLM_predictions_train)) - exp(-exp(-Y)))) %>% 
  round(4)

colnames(GLM_result_df_train) <- c("age_predicted", "age", "species")

GLM_plot_train <- ggplot(GLM_result_df_train, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction GLM (train data set)",
       subtitle = paste0("R=", GLM_predictions_cor_train, " MSE=", GLM_mse_train, " MAE=", GLM_mae_train, " N=", nrow(X))) +
  theme_minimal()

GLM_plot_train + GLM_plot_test +
  plot_layout(nrow = 1)


#### Testing multiple regression models ####

mlm_test <- lm(-log(-log(Y)) ~ X)
summary(mlm_test)

# selecting only significant values
sign_vec <- as.vector((summary(mlm_test)$coefficients[,4] < 0.05)[-1])

mlm_test_significant <- lm(-log(-log(Y)) ~ X[,sign_vec])
summary(mlm_test_significant)

plot(mlm_test)

#### Testing random forest model ####

### Shortcoming of random forest
# Random Forests aren't good at generalizing cases with completely new data. For example, if I tell you that one ice-cream costs $1, 2 ice-creams cost $2, and 3 ice-creams cost $3, how much do 10 ice-creams cost? A linear regression can easily figure this out, while a Random Forest has no way of finding the answer.
# Random forests are biased towards the categorical variable having multiple levels (categories). It is because feature selection based on impurity reduction is biased towards preferring variables with more categories so variable selection (importance) is not accurate for this type of data.

library(randomForest)
set.seed(123)
# change age here!
Y <- all_meth_values_selected_train[,"rel_age"]
Y_test <- all_meth_values_selected_test[,"rel_age"]
#changing value
Y <- -log(-log(Y))
Y_test <- -log(-log(Y_test))

## actual model
RF_test <- randomForest(Y ~ ., data = X, mtry = 9, ntree = 1500)
plot(RF_test)
varImpPlot(RF_test)
importance(RF_test)

which.min(RF_test$mse)

## tuning model
tuneRF(
  x=X, #define predictor variables
  y=Y, #define response variable
  ntreeTry=500,
  mtryStart=4, 
  stepFactor=1.5,
  improve=0.01,
  trace=TRUE #don't show real-time progress
)

### plotting and evaluating
RF_predictions_test <- predict(RF_test, X_test)
RF_predictions_cor_test <- cor(RF_predictions_test, Y_test, method = "pearson") %>% 
  round(4)
RF_mse_test <- mean((RF_predictions_test - Y_test)^2) %>% 
  round(4) 
RF_mae_test <- mean(abs(RF_predictions_test - Y_test)) %>% 
  round(4)

# for normal Y
RF_result_df <- data.frame(predictions = RF_predictions_test,
                            rel_age = Y_test,
                            species = all_meth_values_selected_test$species)
# for transformed Y!

RF_result_df <- data.frame(predictions = exp(-exp(-RF_predictions_test)),
                            rel_age = exp(-exp(-Y_test)),
                            species = all_meth_values_selected_test$species)
# for transformed Y!
RF_mse_test <- mean((exp(-exp(-RF_predictions_test)) - exp(-exp(-Y_test)))^2) %>% 
  round(4) 
RF_mae_test <- mean(abs(exp(-exp(-RF_predictions_test)) - exp(-exp(-Y_test)))) %>% 
  round(4)


colnames(RF_result_df) <- c("age_predicted", "age", "species")

RF_plot_test <- ggplot(RF_result_df, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction RF (test data set)",
       subtitle = paste0("R=", RF_predictions_cor_test,
                         " MSE=", RF_mse_test,
                         " MAE=", RF_mae_test,
                         " N=", nrow(X_test))) +
  theme_minimal()

### running model on training data
RF_predictions_train <- predict(RF_test, newx= X)
RF_predictions_cor_train <- cor(RF_predictions_train, Y, method = "pearson") %>% 
  round(4)
RF_mse_train <- mean((RF_predictions_train - Y)^2) %>% 
  round(4)
RF_mae_train <- mean(abs(RF_predictions_train - Y)) %>% 
  round(4)

# for normal Y
RF_result_df_train <- data.frame(predictions = RF_predictions_train,
                                  rel_age = Y,
                                  species = all_meth_values_selected_train$species)

# for transformed Y!
RF_result_df_train <- data.frame(predictions = exp(-exp(-RF_predictions_train)),
                                  rel_age = exp(-exp(-Y)),
                                  species = all_meth_values_selected_train$species)
# for transformed Y!
RF_mse_train <- mean((exp(-exp(-RF_predictions_train)) - exp(-exp(-Y)))^2) %>% 
  round(4) 
RF_mae_train <- mean(abs(exp(-exp(-RF_predictions_train)) - exp(-exp(-Y)))) %>% 
  round(4)

colnames(RF_result_df_train) <- c("age_predicted", "age", "species")

RF_plot_train <- ggplot(RF_result_df_train, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction RF (train data set)",
       subtitle = paste0("R=", RF_predictions_cor_train, " MSE=", RF_mse_train, " MAE=", RF_mae_train, " N=", nrow(X))) +
  theme_minimal()

RF_plot_train + RF_plot_test +
  plot_layout(nrow = 1)

#### Testing support vector regression models ####
library(e1071)

# change age here!
Y <- all_meth_values_selected_train[,"rel_age"]
Y_test <- all_meth_values_selected_test[,"rel_age"]
#changing value
Y <- -log(-log(Y))
Y_test <- -log(-log(Y_test))

set.seed(123)
SVM_test <- svm(Y ~ ., data = X, type = "eps-regression")
SVM_test <- svm(Y ~ ., data = X, type = "nu-regression")
SVM_test <- svm(Y ~ ., data = X, type = "nu-regression", kernel = "polynomial")

summary(SVM_test)


### plotting
SVM_predictions_test <- predict(SVM_test, X_test)
SVM_predictions_cor_test <- cor(SVM_predictions_test, Y_test, method = "pearson") %>% 
  round(4)
SVM_mse_test <- mean((SVM_predictions_test - Y_test)^2) %>% 
  round(4) 
SVM_mae_test <- mean(abs(SVM_predictions_test - Y_test)) %>% 
  round(4)

# for normal Y
SVM_result_df <- data.frame(predictions = SVM_predictions_test,
                        rel_age = Y_test,
                        species = all_meth_values_selected_test$species)
# for transformed Y!

SVM_result_df <- data.frame(predictions = exp(-exp(-SVM_predictions_test)),
                            rel_age = exp(-exp(-Y_test)),
                            species = all_meth_values_selected_test$species)
# for transformed Y!
SVM_mse_test <- mean((exp(-exp(-SVM_predictions_test)) - exp(-exp(-Y_test)))^2) %>% 
  round(4) 
SVM_mae_test <- mean(abs(exp(-exp(-SVM_predictions_test)) - exp(-exp(-Y_test)))) %>% 
  round(4)


colnames(SVM_result_df) <- c("age_predicted", "age", "species")

SVM_plot_test <- ggplot(SVM_result_df, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction SVM (test data set)",
       subtitle = paste0("R=", SVM_predictions_cor_test,
                         " MSE=", SVM_mse_test,
                         " MAE=", SVM_mae_test,
                         " N=", nrow(X_test))) +
  theme_minimal()

### running model on training data
SVM_predictions_train <- predict(SVM_test, newx= X)
SVM_predictions_cor_train <- cor(SVM_predictions_train, Y, method = "pearson") %>% 
  round(4)
SVM_mse_train <- mean((SVM_predictions_train - Y)^2) %>% 
  round(4)
SVM_mae_train <- mean(abs(SVM_predictions_train - Y)) %>% 
  round(4)

# for normal Y
SVM_result_df_train <- data.frame(predictions = SVM_predictions_train,
                              rel_age = Y,
                              species = all_meth_values_selected_train$species)

# for transformed Y!
SVM_result_df_train <- data.frame(predictions = exp(-exp(-SVM_predictions_train)),
                            rel_age = exp(-exp(-Y)),
                            species = all_meth_values_selected_train$species)
# for transformed Y!
SVM_mse_train <- mean((exp(-exp(-SVM_predictions_train)) - exp(-exp(-Y)))^2) %>% 
  round(4) 
SVM_mae_train <- mean(abs(exp(-exp(-SVM_predictions_train)) - exp(-exp(-Y)))) %>% 
  round(4)

colnames(SVM_result_df_train) <- c("age_predicted", "age", "species")

SVM_plot_train <- ggplot(SVM_result_df_train, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), cex = 3) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction SVM (train data set)",
       subtitle = paste0("R=", SVM_predictions_cor_train, " MSE=", SVM_mse_train, " MAE=", SVM_mae_train, " N=", nrow(X))) +
  theme_minimal()

SVM_plot_train + SVM_plot_test +
  plot_layout(nrow = 1)

#### Testing Bayesian models ####
install.packages("rstanarm")
library(rstanarm)

#### Testing deep learning models ####
install.packages("keras")
library(keras)
DL_test <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(num_features)) %>%
  layer_dense(units = 1)
DL_test %>% compile(
  optimizer = 'rmsprop',
  loss = 'mse',
  metrics = 'mae'
)
DL_test %>% fit(X, Y, epochs = 100, batch_size = 10, validation_split = 0.2)


#### plotting ####
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>%
  .[c(-1,-9)]
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## all
ggplot(cor_all, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  # facet_row(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(cor_all, aes()) +
  geom_point(aes(x = Site, y = Correlation, color = species, alpha = significant)) +
  # geom_line(aes(x = c(-1,1), y = log2(0.05), color = "#CC79A7")) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  facet_wrap(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only significant
ggplot(subset(cor_all, significant == TRUE), aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  # facet_wrap(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only selected 
ggplot(all_sig_CpGs_common, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  facet_row(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only cor positive 
ggplot(all_pos_cor_CpG, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  facet_row(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.ticks.x = element_blank())

## plotting SMR groups 24 and 28
selected_methyl_values <- subset(all_meth_values_long, Site %in% subset(all_sig_CpGs_common, SMR == "SMR_024" | SMR == "SMR_026")$Site)
selected_methyl_values <- subset(all_meth_values_long, Site %in% all_mix_cor_CpG_common$Site)

ggplot(selected_methyl_values, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species, color = species), alpha = 0.9, outlier.size = 0.1) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI) +
  scale_color_manual(values = colpalOI) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)",
       subtitle = "Selected values are correlating with age")

ggplot(selected_methyl_values, aes(x = species, y = Methylation_Value)) +
  geom_sina(aes(color = rel_age, shape = species)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_manual(aesthetics = "legend") +
  theme_classic() +
  # theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)")