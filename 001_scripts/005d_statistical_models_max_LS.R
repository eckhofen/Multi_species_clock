#### Overview ####
# testing different statistical models

#### Settings ####
# change working directory accordingly extension
# setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
extension <- "_max_LS.pdf"


# setting up color palette 
colpal_CB <- c("#c06d00", "#f9cf6e", "#6a5d00", "#44a02b", "#008649", "#12ebf0", "#65a9ff", "#004588", "#660077", "#ff98f7", "#954674", "#630041")
colpal_CB_01 <- colpal_CB[c(TRUE, FALSE)]
colpal_CB_02 <- colpal_CB[c(FALSE, TRUE)]

colpal_CB_a <- c("#f8cbb1","#006786","#182057","#6b6300","#ff8ab9","#f1aaff","#bb005a","#013aa8","#01ef9a","#fa8200","#ee0028","#26c100")
colpal_CB_a_01 <- colpal_CB_a[1:6]
colpal_CB_a_02 <- colpal_CB_a[7:12]

colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)

color_compare <- c("#005AB5", "#DC3220")

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(caret)
library(glmnet)
library(svglite)

#### Loading data ####

load("000_data/006_model_creation/all_meth_values_selected.RData")

#### modifying age (optional) ####

# AC + 30%
extension <- "_max_LS_AC_p_30.pdf"
save_file_AE <- "_max_LS_AC_p_30"
AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age / 1.3

# # AC + 40%
# extension <- "_max_LS_AC_p_40.pdf"
# save_file_AE <- "_max_LS_AC_p_40"
# AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age / 1.4

# AS+ 30%
# extension <- "_max_LS_AS_p_30.pdf"
# save_file_AE <- "_max_LS_AS_p_30"
# AS_meth_values_selected$rel_age <- AS_meth_values_selected$rel_age / 1.3

# EH + 30%
# extension <- "_max_LS_EH_p_30.pdf"
# save_file_AE <- "_max_LS_EH_p_30"
# EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age / 1.3

# EH - 30%
# extension <- "_max_LS_EH_m_25.pdf"
# save_file_AE <- "_max_LS_EH_m_25"
# EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age / 0.75


all_meth_values_selected <- rbind(AC_meth_values_selected,AS_meth_values_selected,EH_meth_values_selected,ZF_meth_values_selected)

#### Data splitting ####
# defining arguments
ds_breaks <- 3
ds_prop <- 5/6
# using a stratified splitting technique
set.seed(123)
AC_split <- initial_split(AC_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)
set.seed(123)
AS_split <- initial_split(AS_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)
set.seed(123)
EH_split <- initial_split(EH_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)
set.seed(123)
ZF_split <- initial_split(ZF_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)

# combining data into training and testing sets
meth_train <- rbind(training(AC_split), training(AS_split), training(EH_split), training(ZF_split))
meth_test <- rbind(testing(AC_split), testing(AS_split), testing(EH_split), testing(ZF_split))

### modified for LOSO

## NO AC
# extension <- "_no_AC.pdf"
# save_file_name <- "no_AC_SVM_eps_val_REL.Rdata"
# 
# meth_train <- rbind(AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected)
# meth_test <- AC_meth_values_selected
# 
# # checking how many CpGs are present per data set
# nrow(meth_train) #261
# nrow(meth_test) #110

## NO ZF
# extension <- "_no_ZF.pdf"
# save_file_name <- "no_ZF_SVM_eps_val.Rdata"
# 
# meth_train <- rbind(AS_meth_values_selected, EH_meth_values_selected, AC_meth_values_selected)
# meth_test <- ZF_meth_values_selected
# 
# # checking how many CpGs are present per data set
# nrow(meth_train) #275
# nrow(meth_test) #96

## NO AS
# extension <- "_no_AS.pdf"
# save_file_name <- "no_AS_SVM_eps_val.Rdata"
# 
# meth_train <- rbind(ZF_meth_values_selected, EH_meth_values_selected, AC_meth_values_selected)
# meth_test <- AS_meth_values_selected
# 
# # checking how many CpGs are present per data set
# nrow(meth_train) #300
# nrow(meth_test) #71

## NO EH
# extension <- "_no_EH.pdf"
# save_file_name <- "no_EH_SVM_eps_val.Rdata"
# 
# meth_train <- rbind(ZF_meth_values_selected, AS_meth_values_selected, AC_meth_values_selected)
# meth_test <- EH_meth_values_selected
# 
# # checking how many CpGs are present per data set
# nrow(meth_train) #277
# nrow(meth_test) #94

## NO EH only AC training
# extension <- "_no_EH_only_AC.pdf"
# save_file_name <- "no_EH_only_AC_SVM_eps_val.Rdata"
# 
# meth_train <- AC_meth_values_selected
# meth_test <- EH_meth_values_selected
# 
# # checking how many CpGs are present per data set
# nrow(meth_train) #110
# nrow(meth_test) #94


### defining data
# training data
X <- meth_train %>% select(-rel_age, -species)

Y <- meth_train[,"rel_age"]

# testing data
X_test <- meth_test %>% select(-rel_age, -species)

Y_test <- meth_test[,"rel_age"]

#### model evaluation ####

## function for model testing plus calculating metrics (MSE, MAE, R)
evaluate.model <- function(model, X_train, Y_train, X_test, Y_test, species_train,
                           species_test, transform = FALSE, colpalOI = color_species, plot_title = "Model evaluation:", 
                           y_lim = c(0,.3), x_lim = c(-0,.3), CpGs = "not defined", s = NA) {
  # calculate predictions
  if (!is.na(s)) {
    predictions_train <- predict(model, X_train, s = s)
    predictions_test <- predict(model, X_test, s = s)
  } else {
    predictions_train <- predict(model, X_train)
    predictions_test <- predict(model, X_test)
  }
  
  # optional transformation
  if (transform) {
    predictions_train <- exp(-exp(-predictions_train))
    predictions_test <- exp(-exp(-predictions_test))
    Y_train
    Y_test
  }
  
  # calculate metrics
  metrics_train <- data.frame(
    R = round(cor(predictions_train, Y_train, method = "pearson"), 4),
    MSE = round(mean((predictions_train - Y_train)^2), 4),
    MAE = round(mean(abs(predictions_train - Y_train)), 4),
    N = nrow(X_train)
  )
  
  metrics_test <- data.frame(
    R = round(cor(predictions_test, Y_test, method = "pearson"), 4),
    MSE = round(mean((predictions_test - Y_test)^2), 4),
    MAE = round(mean(abs(predictions_test - Y_test)), 4),
    N = nrow(X_test)
  )
  values_AE_train <- data.frame(
    AE = round(abs(predictions_train - Y_train), 4)
  )
  values_AE_test <- data.frame(
    AE = round(abs(predictions_test - Y_test), 4)
  )
  
  # students t-test for training and testing estimation errors
  metrics_ttest <- t.test(values_AE_train, values_AE_test)
  
  # prepare data frames for plotting
  if(is.matrix(predictions_train)) {
    result_df_train <- data.frame(age_predicted = c(predictions_train), age = Y_train, species = species_train)
    result_df_test <- data.frame(age_predicted = c(predictions_test), age = Y_test, species = species_test)
  } else {
    result_df_train <- data.frame(age_predicted = predictions_train, age = Y_train, species = species_train)
    result_df_test <- data.frame(age_predicted = predictions_test, age = Y_test, species = species_test)
  }
  
  
  # create plots
  plot_train <- ggplot(result_df_train, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_smooth(method = "lm", se = FALSE, color = "black") +
    # geom_line(data = result_df_train, aes(x = age, y = age_predicted), color = "black", alpha = 0.7, linetype = "dashed") +
    scale_color_manual(values = colpalOI) +
    ylim(y_lim) +
    xlim(x_lim) +
    labs(title = paste(plot_title, "(Training Set)"), y = "Estimated age", x = "Relative age", color = "Species") +
    # labs(subtitle = paste0("R=", metrics_test$R, "\nMSE=", metrics_test$MSE, "\nMAE=", metrics_test$MAE, "\nN=", nrow(X_test), " CpGs=", CpGs)) +
    theme_bw() +
    annotate("text", x = 0, y = 0.25, label = paste0("R=", metrics_train$R, "\nMSE=", metrics_train$MSE, "\nMAE=", metrics_train$MAE), size = 3.5, hjust = 0, vjust = .1) +
    theme(plot.title = element_text(hjust = .5, face = "bold"))
  
  plot_test <- ggplot(result_df_test, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_manual(values = colpalOI) +
    ylim(y_lim) +
    xlim(x_lim) +
    labs(title = paste(plot_title, "(Testing Set)"), y = "Estimated age", x = "Relative age", color = "Species") +
    # labs(subtitle = paste0("R=", metrics_test$R, "\nMSE=", metrics_test$MSE, "\nMAE=", metrics_test$MAE, "\nN=", nrow(X_test), " CpGs=", CpGs)) +
    theme_bw() +
    annotate("text", x = 0, y = 0.25, label = paste0("R=", metrics_test$R, "\nMSE=", metrics_test$MSE, "\nMAE=", metrics_test$MAE), size = 3.5, hjust = 0, vjust = .1) +
    theme(plot.title = element_text(hjust = .5, face = "bold"))
  
  # return list containing metrics and plots
  return(list(metrics_train = metrics_train, metrics_test = metrics_test, plot_train = plot_train, plot_test = plot_test, values_AE_train = values_AE_train, values_AE_test = values_AE_test, t_test = metrics_ttest))
}


#### CV ENR ####

# change age here!
Y <- meth_train[,"rel_age"]
Y_test <- meth_test[,"rel_age"]

# age transformation
Y_log <- -log(-log(Y))
Y_log_test <- -log(-log(Y_test))

# define alpha for either lasso, rigid or elastic net regression
glm_alpha <- 0.5

### 1) normal model
# setting seed for reproducibility 
set.seed(123)

## model
ENR_test <- cv.glmnet(as.matrix(X), Y, alpha = glm_alpha, family = "gaussian")
plot(ENR_test)
coef(ENR_test, s=ENR_test$lambda.min)

### 2) log transformed
# setting seed for reproducibility 
set.seed(123)

## model
ENR_test_log <- cv.glmnet(as.matrix(X), Y_log, alpha = glm_alpha)

### running models on testing data
ENR_eval_rel <-  evaluate.model(ENR_test, s = ENR_test$lambda.min, as.matrix(X), Y, as.matrix(X_test), Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "ENR", CpGs = "40")

ENR_eval_rel_t <-  evaluate.model(ENR_test_log, s = ENR_test_log$lambda.min, as.matrix(X), Y, as.matrix(X_test), Y_test, meth_train$species, meth_test$species, transform = TRUE, 
                                  colpalOI= color_species, plot_title = "ENR (-log-log(age))", CpGs = "40")
##grouped
# normal age
ENR_eval_rel_plot <- ENR_eval_rel$plot_train + ENR_eval_rel$plot_test  + ENR_eval_rel_t$plot_train + ENR_eval_rel_t$plot_test +
  plot_layout(nrow = 2, guides = "collect")

# saving plots
# ggsave(filename = paste0("002_plots/005_m_ENR_rel_age_TE", extension), ENR_eval_rel$plot_test, width = 8, height = 7)
# ggsave(filename = paste0("002_plots/005_m_ENR_rel_age_TR", extension), ENR_eval_rel$plot_train, width = 8, height = 7)
# 
# ggsave(filename = paste0("002_plots/005_m_ENR_rel_log-age_TE", extension), ENR_eval_rel_t$plot_test, width = 8, height = 7)
# ggsave(filename = paste0("002_plots/005_m_ENR_rel_log-age_TR", extension), ENR_eval_rel_t$plot_train, width = 8, height = 7)

ggsave(filename = paste0("002_plots/005_m_ENR_rel_age_all", extension), ENR_eval_rel_plot, width = 10, height = 7)

#### Testing multivariate linear regression models ####
### pre-testing with base R package
mlm_test <- lm(Y ~., data = X)
summary(mlm_test)
coef(mlm_test)


# with transformed age
mlm_test_t <- lm(Y_log ~., data = X)

## selecting only significant values
sign_vec <- as.vector((summary(mlm_test)$coefficients[,4] < 0.05)[-1])

# stat tests 
plot(density(mlm_test$residuals))
shapiro.test(mlm_test$residuals)

## selecting only significant values
mlm_test_opt <- lm(Y ~., data = X[,sign_vec])
summary(mlm_test_opt)

# with transformed age
mlm_test_opt_t <- lm(Y_log ~.,X[sign_vec])

## using centered and scaled dataset 
# did not change anything!
# mlm_alt <- lm(Y ~., trainingData)
# predict(mlm_alt, testingData)

### testing and plotting model

# testing normal mlm
mlm_eval_rel <-  evaluate.model(mlm_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, 
                                colpalOI= color_species, plot_title = "MLM", CpGs = length(mlm_test$coefficients)-1)
# mlm_alt_eval_rel <-  evaluate.model(mlm_alt, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= colpal_CB, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)

mlm_eval_rel_t <-  evaluate.model(mlm_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, 
                                  colpalOI= color_species, plot_title = "MLM (-log-log(age))", CpGs = length(mlm_test_t$coefficients)-1)

# testing optimized mlm
mlm_eval_rel_opt <-
  evaluate.model(mlm_test_opt, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI = color_species, plot_title = "MLM (sign. CpGs only)", CpGs = length(mlm_test_opt$coefficients) - 1)
mlm_eval_rel_opt_t <-
  evaluate.model(mlm_test_opt_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI = color_species, plot_title = "MLM (sig. only; -log-log(age))", CpGs = length(mlm_test_opt$coefficients) - 1)

## plotting 
# normal age
mlm_eval_rel_plot <- mlm_eval_rel$plot_train + mlm_eval_rel$plot_test  + mlm_eval_rel_opt$plot_train + mlm_eval_rel_opt$plot_test +
  plot_layout(nrow = 2, guides = "collect")
# transformed age
mlm_eval_rel_plot_t <- mlm_eval_rel_t$plot_train + mlm_eval_rel_t$plot_test  + mlm_eval_rel_opt_t$plot_train + mlm_eval_rel_opt_t$plot_test +
  plot_layout(nrow = 2, guides = "collect")

ggsave(filename = paste0("002_plots/005_m_MLM_rel_age_all", extension), mlm_eval_rel_plot, width = 10, height = 7)
ggsave(filename = paste0("002_plots/005_m_MLM_rel_log-age_all", extension), mlm_eval_rel_plot_t, width = 10, height = 7)


#### ENR (caret) ####
## pre processing
# centering means that all the the mean is being subtracted from all values and scale divides them by the standard deviation. 

# training data
preProcValues <- preProcess(X, method = c("center", "scale"))
trainingData <- predict(preProcValues, X)
trainingData$rel_age <- Y
# testing data
preProcValues_test <- preProcess(X_test, method = c("center", "scale"))
testingData <- predict(preProcValues_test, X_test)
testingData$rel_age <- Y_test

# transformation
testingData_t <- testingData
trainingData_t <- trainingData
# -log-log transformation
trainingData_t$rel_age <- -log(-log(trainingData_t$rel_age))
testingData_t$rel_age <- -log(-log(testingData_t$rel_age))

### training model 
set.seed(123)
trainControl <- trainControl(method = "cv", number = 10) # 10-fold CV
MLM_model <- train(rel_age ~ ., data = trainingData, method = "lm", trControl = trainControl) #names(getModelInfo()) for all method variables

#tuning model 
MLM_tuned <- train(rel_age ~ ., data = testingData,  method = "lm", tuneLength = 10, trControl = trainControl)

importance <- varImp(MLM_tuned, scale = FALSE)
plot(importance)

### model transformed
set.seed(123)
trainControl <- trainControl(method = "cv", number = 10) # 10-fold CV
MLM_model_t <- train(rel_age ~ ., data = trainingData_t, method = "lm", trControl = trainControl) #names(getModelInfo()) for all method variables

#tuning model 
MLM_tuned_t <- train(rel_age ~ ., data = testingData_t,  method = "lm", tuneLength = 10, trControl = trainControl)

## evaluation 
# normal models
MLM_eval_rel <-  evaluate.model(MLM_model, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_eval_rel_tuned <-  evaluate.model(MLM_tuned, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "MLM tuned prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_eval_rel$plot_train + MLM_eval_rel$plot_test + MLM_eval_rel_tuned$plot_train + MLM_eval_rel_tuned$plot_test +
  plot_layout(nrow=2, guides = "collect")

# transformed models
MLM_t_eval_rel <-  evaluate.model(MLM_model_t, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "MLM (-log-log(age)) prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_t_eval_rel_tuned <-  evaluate.model(MLM_tuned_t, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "MLM (-log-log(age)) tuned prediction", CpGs = length(mlm_test$coefficients)-1)

### LOOCV

# preparing data 
all_data <- all_meth_values_selected[,-length(all_meth_values_selected)]
# all_data$rel_age <- -log(-log(all_data$rel_age))
# all_age <- all_meth_values_selected$rel_age %>% -log(-log(.))
all_data <- rbind(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected) %>% 
  .[,-length(AC_meth_values_selected)]

all_data_t <- all_data

all_data_t$rel_age <- log(all_data_t$rel_age)
# setting up training method
loocv_train_control <- trainControl(method = "LOOCV")

# run model
set.seed(123)
MLM_LOOCV_model <- train(rel_age ~ ., data = all_data, method = "lm", trControl = loocv_train_control)
MLM_LOOCV_model_t <- train(rel_age ~ ., data = all_data_t, method = "lm", trControl = loocv_train_control)

MLM_LOOCV_model$results
MLM_LOOCV_model_t$results

# evaluation
metrics_LOOCV <- data.frame(
  R = round(cor(MLM_LOOCV_model$pred$pred, MLM_LOOCV_model$pred$obs, method = "pearson"), 4),
  MSE = round(mean((MLM_LOOCV_model$pred$pred - MLM_LOOCV_model$pred$obs)^2), 4),
  MAE = round(mean(abs(MLM_LOOCV_model$pred$pred - MLM_LOOCV_model$pred$obs)), 4),
  N = nrow(all_data))

values_LOOCV_AE <- data.frame(
  AE = round(abs(MLM_LOOCV_model$pred$pred - MLM_LOOCV_model$pred$obs), 4)
)

metrics_LOOCV_t <- data.frame(
  R = round(cor(MLM_LOOCV_model_t$pred$pred, MLM_LOOCV_model_t$pred$obs, method = "pearson"), 4),
  MSE = round(mean((MLM_LOOCV_model_t$pred$pred - MLM_LOOCV_model_t$pred$obs)^2), 4),
  MAE = round(mean(abs(MLM_LOOCV_model_t$pred$pred - MLM_LOOCV_model_t$pred$obs)), 4),
  N = nrow(all_data))

values_LOOCV_t_AE <- data.frame(
  AE = round(abs(MLM_LOOCV_model_t$pred$pred - MLM_LOOCV_model_t$pred$obs), 4)
)

# evaluate model
# MLM_LOOCV_eval_rel <-  evaluate.model(MLM_LOOCV_model, all_data, all_age, testingData, Y_test, all_meth_values_selected$species, meth_test$species, transform = FALSE, colpalOI= colpal_CB_a_01, plot_title = "MLM LOOCV prediction", CpGs = length(mlm_test$coefficients)-1)

#### Testing random forest model ####

### Shortcoming of random forest
# Random Forests aren't good at generalizing cases with completely new data. For example, if I tell you that one ice-cream costs $1, 2 ice-creams cost $2, and 3 ice-creams cost $3, how much do 10 ice-creams cost? A linear regression can easily figure this out, while a Random Forest has no way of finding the answer.
# Random forests are biased towards the categorical variable having multiple levels (categories). It is because feature selection based on impurity reduction is biased towards preferring variables with more categories so variable selection (importance) is not accurate for this type of data.

library(randomForest)
set.seed(123)

## actual model
RF_test <- randomForest(Y ~ ., data = X, mtry = 4, ntree = 1500)
plot(RF_test)
varImpPlot(RF_test)
importance(RF_test)

which.min(RF_test$mse)

set.seed(123)
RF_test_t <- randomForest(Y_log ~ ., data = X, mtry = 4, ntree = 1500)

## tuning model
set.seed(123)
tuneRF(
  x=X, #define predictor variables
  y=Y, #define response variable
  ntreeTry=1500,
  mtryStart=4, 
  stepFactor=1.5,
  improve=0.01,
  trace=TRUE
)

# take suggested mtry
set.seed(123)
RF_test_tuned <- randomForest(Y ~ ., data = X, mtry = 9, ntree = 1500)

## evaluation
RF_eval_rel <-  evaluate.model(RF_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "RF", CpGs = length(mlm_test$coefficients)-1)

RF_eval_rel_t <-  evaluate.model(RF_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "RF (-log-log(age))", CpGs = length(mlm_test$coefficients)-1)

RF_eval_rel_tuned <-  evaluate.model(RF_test_tuned, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "RF (tuned)", CpGs = length(mlm_test$coefficients)-1)

RF_eval_rel_plot <- RF_eval_rel$plot_train + RF_eval_rel$plot_test + RF_eval_rel_tuned$plot_train + RF_eval_rel_tuned$plot_test +
  plot_layout(nrow=2, guides = "collect")

ggsave(filename = paste0("002_plots/005_m_RF_rel_age_all", extension), RF_eval_rel_plot, width = 10, height = 7)

#### Testing support vector regression models ####
library(e1071)
# normal models
SVM_eps_test <- svm(Y ~ ., data = X, type = "eps-regression", kernel = "linear")
SVM_nu_test <- svm(Y ~ ., data = X, type = "nu-regression", kernel = "linear")
# transformed
SVM_eps_test_t <- svm(Y_log ~ ., data = X, type = "eps-regression", kernel = "linear")
SVM_nu_test_t <- svm(Y_log ~ ., data = X, type = "nu-regression", kernel = "linear")
# SVM_test <- svm(Y ~ ., data = X, type = "nu-regression")
# SVM_test <- svm(Y ~ ., data = X, type = "nu-regression", kernel = "polynomial", degree = 3)
# SVM_test <- svm(Y ~ ., data = X, type = "nu-regression", kernel = "polynomial", degree = 1)

# summary(SVM_test)
# plot(abs(coef(SVM_test)))

## evaluation
# normal
SVM_eps_eval_rel <-  evaluate.model(SVM_eps_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "SVM (eps, linear)", CpGs = length(mlm_test$coefficients)-1)

SVM_nu_eval_rel <-  evaluate.model(SVM_nu_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "SVM (nu, linear)", CpGs = length(mlm_test$coefficients)-1)

SVM_eval_rel_plot <- SVM_eps_eval_rel$plot_train + SVM_eps_eval_rel$plot_test + SVM_nu_eval_rel$plot_train + SVM_nu_eval_rel$plot_test +
  plot_layout(nrow = 2, guides = "collect")
# transformed
SVM_eps_eval_rel_t <-  evaluate.model(SVM_eps_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "SVM (eps, linear) log", CpGs = length(mlm_test$coefficients)-1)

SVM_nu_eval_rel_t <-  evaluate.model(SVM_nu_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "SVM (nu, linear) log", CpGs = length(mlm_test$coefficients)-1)

SVM_eval_rel_t_plot <- SVM_eps_eval_rel_t$plot_train + SVM_eps_eval_rel_t$plot_test + SVM_nu_eval_rel_t$plot_train + SVM_nu_eval_rel_t$plot_test +
  plot_layout(nrow = 2, guides = "collect")

ggsave(filename = paste0("002_plots/005_m_SVM_rel_age_all", extension), SVM_eval_rel_plot, width = 10, height = 7)
ggsave(filename = paste0("002_plots/005_m_SVM_rel_log-age_all", extension), SVM_eval_rel_t_plot, width = 10, height = 7)

#### AE comparison plot ####

## all AE values for plotting REL age
# train
df_AE_training <- cbind(ENR_eval_rel$values_AE_train, mlm_eval_rel$values_AE_train, SVM_eps_eval_rel$values_AE_train, SVM_nu_eval_rel$values_AE_train, RF_eval_rel$values_AE_train, Type = as.factor("Training"))
colnames(df_AE_training) <- c("ENR", "MLM", "SVM_eps", "SVM_nu", "RF", "type")
# test
df_AE_testing <- cbind(ENR_eval_rel$values_AE_test, mlm_eval_rel$values_AE_test, SVM_eps_eval_rel$values_AE_test, SVM_nu_eval_rel$values_AE_test, RF_eval_rel$values_AE_test,Type = as.factor("Testing"))
colnames(df_AE_testing) <- c("ENR", "MLM", "SVM_eps", "SVM_nu", "RF","type")
# both
df_AE <- rbind(df_AE_training, df_AE_testing)

# long transformation
df_AE_long_REL <- as.data.frame(pivot_longer(df_AE, cols = -type, names_to = "model", values_to = "value"))
colnames(df_AE_long_REL) <- c("type", "model", "AE")

## all AE values for plotting LOGLOG REL age
# train
df_AE_training_t <- cbind(ENR_eval_rel_t$values_AE_train, mlm_eval_rel_t$values_AE_train, SVM_eps_eval_rel_t$values_AE_train, SVM_nu_eval_rel_t$values_AE_train, RF_eval_rel_t$values_AE_train, Type = as.factor("Training"))
colnames(df_AE_training_t) <- c("ENR", "MLM", "SVM_eps", "SVM_nu", "RF", "type")
# test
df_AE_testing_t <- cbind(ENR_eval_rel_t$values_AE_test, mlm_eval_rel_t$values_AE_test, SVM_eps_eval_rel_t$values_AE_test, SVM_nu_eval_rel_t$values_AE_test, RF_eval_rel_t$values_AE_test,Type = as.factor("Testing"))
colnames(df_AE_testing_t) <- c("ENR", "MLM", "SVM_eps", "SVM_nu", "RF","type")
# both
df_AE_t <- rbind(df_AE_training_t, df_AE_testing_t)

save(file= paste0("000_data/007_data_comparison/df_AE", save_file_AE), df_AE)
save(file= paste0("000_data/007_data_comparison/df_AE_t", save_file_AE), df_AE_t)

# long transformation
df_AE_long_t <- as.data.frame(pivot_longer(df_AE_t, cols = -type, names_to = "model", values_to = "value"))
colnames(df_AE_long_t) <- c("type", "model", "AE") 

# MLM only
df_AE_MLM <- df_AE[,c("MLM","type")]
df_AE_MLM_t <- df_AE_t[,c("MLM","type")]

# SVM eps only 
df_AE_SVM_t <- df_AE_t[,c("SVM_eps","type")]

df_AE_long_both <- rbind(df_AE_long, df_AE_long_t)

df_AE_LOOCV <- cbind(values_LOOCV_AE$AE, values_LOOCV_t_AE$AE)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colpal <- hcl.colors(7, "SunsetDark") 

plot_AE_comparison <- ggplot(df_AE_long, aes(x = model, y = AE, fill = model, pattern = type)) +
  geom_boxplot_pattern(
    position = position_dodge(width = .9),
    outlier.alpha = 0.7,
    pattern_fill = "transparent",
    pattern_color = "gray10") + 
  scale_fill_manual(values = colpal) +
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Models", fill = "Model", pattern = "Data") +
  theme_bw() + 
  ylim(c(0,9)) +
  theme(axis.text.x = element_text(color = "black"))

plot_AE_comparison

plot_AE_comparison_t <- ggplot(df_AE_long_t, aes(x = model, y = AE, fill = model, pattern = type)) +
  geom_boxplot_pattern(
    position = position_dodge(width = .9),
    outlier.alpha = 0.7,
    pattern_fill = "transparent",
    pattern_color = "gray10") + 
  scale_fill_manual(values = colpal) +
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Models", fill = "Model", pattern = "Data") +
  theme_bw() + 
  ylim(c(0,9)) +
  theme(axis.text.x = element_text(color = "black"))

plot_AE_comparison_t

ggsave(filename = paste0("002_plots/005_comparison_AE_age", extension), plot_AE_comparison, width = 7, height = 5)
ggsave(filename = paste0("002_plots/005_comparison_AE_log_age", extension), plot_AE_comparison_t, width = 7, height = 5)

## no outliers
plot_AE_comparison_no_out <- ggplot(df_AE_long, aes(x = model, y = AE, fill = model, pattern = type)) +
  geom_boxplot_pattern(
    position = position_dodge(width = .9),
    outlier.shape = NA,
    pattern_fill = "transparent",
    pattern_color = "gray10") + 
  scale_fill_manual(values = colpal) +
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Models", fill = "Model", pattern = "Data") +
  theme_bw() + 
  ylim(c(0,3.5)) +
  theme(axis.text.x = element_text(color = "black"))

plot_AE_comparison_no_out

plot_AE_comparison_no_out_t <- ggplot(df_AE_long_t, aes(x = model, y = AE, fill = model, pattern = type)) +
  geom_boxplot_pattern(
    position = position_dodge(width = .9),
    outlier.shape = NA,
    pattern_fill = "transparent",
    pattern_color = "gray10") + 
  scale_fill_manual(values = colpal) +
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Models", fill = "Model", pattern = "Data") +
  theme_bw() + 
  ylim(c(0,3.5)) +
  theme(axis.text.x = element_text(color = "black"))

plot_AE_comparison_no_out_t

ggsave(filename = paste0("002_plots/005_comparison_AE_age_no_out", extension), plot_AE_comparison_no_out, width = 7, height = 5)
ggsave(filename = paste0("002_plots/005_comparison_AE_log_age_no_out", extension), plot_AE_comparison_no_out_t, width = 7, height = 5)

## AE only for SVM

plot_AE_SVM_t <- ggplot(df_AE_SVM_t, aes(x = type, y = SVM_eps, pattern = type)) +
  geom_boxplot_pattern(
    position = position_dodge(width = .9),
    outlier.shape = NA,
    pattern_fill = "transparent",
    pattern_color = "gray10") + 
  scale_fill_manual(values = colpal) +
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Data", pattern = "Data") +
  theme_bw() + 
  theme(legend.position = "") +
  ylim(c(0,3)) 
  # theme(axis.text.x = element_blank())
plot_AE_SVM_t
#### Final plots ####
# This only works when both "005_stati...R" script have been run and the environments of both are still loaded

## SVM comparison
SVM_comparison_plot <- SVM_eps_eval_chron_t$plot_train + SVM_eps_eval_chron_t$plot_test + plot_AE_SVM_t +
  plot_layout(nrow = 1, guides = "collect", width = c(3,3,1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

SVM_comparison_plot
ggsave(filename = paste0("002_plots/005_comparison_SVM_eps_t", extension), SVM_comparison_plot, width = 10, height = 4)

## AE comparison
AE_comparison_plot <- plot_AE_comparison_no_out + plot_AE_comparison_no_out_t +
  plot_layout(nrow = 2, guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"), legend.position = "")

AE_comparison_plot
ggsave(filename = paste0("002_plots/005_comparison_AE", extension), AE_comparison_plot, width = 10, height = 7)


### all models comparison

all_models_plot <- ENR_eval_rel$plot_train + ENR_eval_rel$plot_test + ENR_eval_rel_t$plot_train + ENR_eval_rel_t$plot_test + 
  mlm_eval_rel$plot_train + mlm_eval_rel$plot_test + mlm_eval_rel_t$plot_train + mlm_eval_rel_t$plot_test +
  RF_eval_rel$plot_train + RF_eval_rel$plot_test + RF_eval_rel_t$plot_train + RF_eval_rel_t$plot_test +
  SVM_eps_eval_rel$plot_train + SVM_eps_eval_rel$plot_test +
  SVM_eps_eval_rel_t$plot_train + SVM_eps_eval_rel_t$plot_test +
  SVM_nu_eval_rel$plot_train + SVM_nu_eval_rel$plot_test +
  SVM_nu_eval_rel_t$plot_train + SVM_nu_eval_rel_t$plot_test +
  plot_layout(ncol = 4, guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0("002_plots/005_all_models_plot_REL", extension), all_models_plot, width = 15, height = 18)
