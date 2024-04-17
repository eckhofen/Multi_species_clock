#### Overview ####
# testing different statistical models

#### Settings ####
# change working directory accordingly extension
# setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
extension <- ".pdf"

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
library(patchwork)
library(caret)
library(glmnet)
library(svglite)

#### Loading data ####

load("000_data/006_model_creation/all_meth_values_selected.RData")

#### modifying age (optional) ####
# AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age / 1.4
# AS_meth_values_selected$rel_age
# EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age / 1
# ZF_meth_values_selected$rel_age

AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age * 25
AS_meth_values_selected$rel_age <- AS_meth_values_selected$rel_age * 54
EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age * 20
ZF_meth_values_selected$rel_age <- ZF_meth_values_selected$rel_age * 5

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
meth_train <- rbind(training(AC_split), training(AS_split), training(EH_split), training(ZF_split))
# meth_train <- rbind(training(ZF_split))

meth_test <- rbind(testing(AC_split), testing(AS_split), testing(EH_split), testing(ZF_split))
# meth_test <- rbind(testing(AC_split), testing(AS_split), testing(EH_split), ZF_meth_values_selected)
# meth_test <- rbind(training(AC_split), training(AS_split), training(EH_split))

# checking how many CpGs are present per data set
nrow(meth_train) #272
nrow(meth_test) #99

## plotting age distribution
plot_age_dist <- ggplot(all_meth_values_selected) +
  geom_density(aes(x = rel_age, color = species), linewidth = 1) +
  geom_density(aes(x = rel_age, fill = "all"), alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = color_species) +
  scale_fill_manual(values = "grey") +
  theme_minimal() +
  labs(title = "Relative age distribution for all samples")

# >> the plot shows that our dependent variable is not normally distributed (as expected from age)

# comparing data sets with Kolmogorov-Smirnov test
ks_test_data <- ks.test(meth_train$rel_age, meth_test$rel_age) # D = 0.058304, p-value = 0.966

ks_test_data$statistic

color_compare_tt <- setNames(color_compare, c("Training", "Testing"))
# visually comparing training and testing sets (BOXPLOTS)
plot_sample_age_dist_box <- ggplot() +
  geom_boxplot(data = meth_train, aes(y = rel_age, fill = "Training", x = -.5)) +
  geom_boxplot(data = meth_test, aes(y = rel_age, fill = "Testing", x = .5)) +
  scale_fill_manual(values = color_compare_tt) +
  labs(fill = "Dataset") +
  xlab("Datasets") +
  ylab("Age") +
  theme_classic() +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Age Distribution in Training and Testing Sets", subtitle = paste0("KS-test: D=", round(ks_test_data$statistic, 4), " p-value=", ks_test_data$p.value))

# visually comparing training and testing sets
plot_sample_age_dist <- ggplot() +
  geom_density(data = meth_train, aes(x = rel_age, fill = "Training"), alpha = 0.5, linewidth = NA) +
  geom_density(data = meth_test, aes(x = rel_age, fill = "Testing"), alpha = 0.5, linewidth = NA) +
  scale_fill_manual(values = color_compare_tt, name = "Dataset") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Age Distribution in Training and Testing Sets")

# visually comparing sets (Q-Q plot)
qqplot(meth_train$rel_age, meth_test$rel_age,
       xlab = "Training data",
       ylab = "Testing data",
       main = "Age Distribution in Training and Testing Sets")
abline(0, 1, col = "black")


# plotting both graphs
ggsave(filename = paste0("002_plots/005_age_distribution", extension), plot_age_dist, width = 7, height = 7)
ggsave(filename = paste0("002_plots/005_sample_age_distribution_box", extension), plot_sample_age_dist_box, width = 7, height = 7)

plot_age_dist + plot_sample_age_dist +
  plot_layout(nrow=1)

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
                           y_lim = c(0,12), x_lim = c(0,12), CpGs = "not defined", s = NA) {
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
    predictions_train <- exp(predictions_train)
    predictions_test <- exp(predictions_test)
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
    labs(title = paste(plot_title, "(Training Set)"), y = "Estimated age", x = "Chronological age",
         subtitle = paste0("R=", metrics_train$R, " MSE=", metrics_train$MSE, " MAE=", metrics_train$MAE, " N=", nrow(X_train), " CpGs=", CpGs)) +
    theme_classic()
  
  plot_test <- ggplot(result_df_test, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_manual(values = colpalOI) +
    ylim(y_lim) +
    xlim(x_lim) +
    labs(title = paste(plot_title, "(Testing Set)"), y = "Estimated age", x = "Chronological age",
         subtitle = paste0("R=", metrics_test$R, " MSE=", metrics_test$MSE, " MAE=", metrics_test$MAE, " N=", nrow(X_test), " CpGs=", CpGs)) +
    theme_classic()
  
  # return list containing metrics and plots
  return(list(metrics_train = metrics_train, metrics_test = metrics_test, plot_train = plot_train, plot_test = plot_test))
}


#### CV GLM ####

# change age here!
Y <- meth_train[,"rel_age"]
Y_test <- meth_test[,"rel_age"]

# age transformation
Y_log <- log(Y)
Y_log_test <- log(Y_test)

# define alpha for either lasso, rigid or elastic net regression
glm_alpha <- 0.5

### 1) normal model
# setting seed for reproducibility 
set.seed(123)

## model
GLM_test <- cv.glmnet(as.matrix(X), Y, alpha = glm_alpha, family = "gaussian")
plot(GLM_test)
coef(GLM_test, s=GLM_test$lambda.min)

### 2) log transformed
# setting seed for reproducibility 
set.seed(123)

## model
GLM_test_log <- cv.glmnet(as.matrix(X), Y_log, alpha = glm_alpha)

### running models on testing data
GLM_eval <-  evaluate.model(GLM_test, s = GLM_test$lambda.min, as.matrix(X), Y, as.matrix(X_test), Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "GLM", CpGs = "40")

GLM_eval_t <-  evaluate.model(GLM_test_log, s = GLM_test_log$lambda.min, as.matrix(X), Y, as.matrix(X_test), Y_test, meth_train$species, meth_test$species, transform = TRUE, 
                              colpalOI= color_species, plot_title = "GLM (log(age))", CpGs = "40")

# saving plots
ggsave(filename = paste0("002_plots/005_m_GLM_age_TE", extension), GLM_eval$plot_test, width = 8, height = 7)
ggsave(filename = paste0("002_plots/005_m_GLM_age_TR", extension), GLM_eval$plot_train, width = 8, height = 7)

ggsave(filename = paste0("002_plots/005_m_GLM_log-age_TE", extension), GLM_eval_t$plot_test, width = 8, height = 7)
ggsave(filename = paste0("002_plots/005_m_GLM_log-age_TR", extension), GLM_eval_t$plot_train, width = 8, height = 7)

#### Testing multivariate linear regression models ####
### pre-testing with base R package
mlm_test <- lm(Y ~., data = X)
summary(mlm_test)
coef(mlm_test)
# with transformed age
mlm_test_t <- lm(log(Y) ~., data = X)

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
mlm_eval <-  evaluate.model(mlm_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, 
                            colpalOI= color_species, plot_title = "MLM", CpGs = length(mlm_test$coefficients)-1)
# mlm_alt_eval <-  evaluate.model(mlm_alt, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= colpal_CB, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)

mlm_eval_t <-  evaluate.model(mlm_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, 
                              colpalOI= color_species, plot_title = "MLM (log(age))", CpGs = length(mlm_test_t$coefficients)-1)

# testing optimized mlm
mlm_eval_opt <-
  evaluate.model(mlm_test_opt, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI = color_species, plot_title = "MLM (sign. CpGs only)", CpGs = length(mlm_test_opt$coefficients) - 1)
mlm_eval_opt_t <-
  evaluate.model(mlm_test_opt_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI = color_species, plot_title = "MLM (sig. only; log(age))", CpGs = length(mlm_test_opt$coefficients) - 1)

## plotting 
# normal age
mlm_eval_plot <- mlm_eval$plot_train + mlm_eval$plot_test  + mlm_eval_opt$plot_train + mlm_eval_opt$plot_test +
  plot_layout(nrow = 2)
# transformed age
mlm_eval_plot_t <- mlm_eval_t$plot_train + mlm_eval_t$plot_test  + mlm_eval_opt_t$plot_train + mlm_eval_opt_t$plot_test +
  plot_layout(nrow = 2)

ggsave(filename = paste0("002_plots/005_m_MLM_age_all", extenstion), mlm_eval_plot, width = 10, height = 7)
ggsave(filename = paste0("002_plots/005_m_MLM_log-age_all", extenstion), mlm_eval_plot_t, width = 10, height = 7)


#### GLM (caret) ####
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
trainingData_t$rel_age <- log(trainingData_t$rel_age)
testingData_t$rel_age <- log(testingData_t$rel_age)

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
MLM_eval <-  evaluate.model(MLM_model, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_eval_tuned <-  evaluate.model(MLM_tuned, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "MLM tuned prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_eval$plot_train + MLM_eval$plot_test + MLM_eval_tuned$plot_train + MLM_eval_tuned$plot_test +
  plot_layout(nrow=2)

# transformed models
MLM_t_eval <-  evaluate.model(MLM_model_t, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "MLM (-log-log(age)) prediction", CpGs = length(mlm_test$coefficients)-1)

MLM_t_eval_tuned <-  evaluate.model(MLM_tuned_t, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "MLM (-log-log(age)) tuned prediction", CpGs = length(mlm_test$coefficients)-1)

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

metrics_LOOCV_t <- data.frame(
  R = round(cor(MLM_LOOCV_model_t$pred$pred, MLM_LOOCV_model_t$pred$obs, method = "pearson"), 4),
  MSE = round(mean((MLM_LOOCV_model_t$pred$pred - MLM_LOOCV_model_t$pred$obs)^2), 4),
  MAE = round(mean(abs(MLM_LOOCV_model_t$pred$pred - MLM_LOOCV_model_t$pred$obs)), 4),
  N = nrow(all_data))

# evaluate model
# MLM_LOOCV_eval <-  evaluate.model(MLM_LOOCV_model, all_data, all_age, testingData, Y_test, all_meth_values_selected$species, meth_test$species, transform = FALSE, colpalOI= colpal_CB_a_01, plot_title = "MLM LOOCV prediction", CpGs = length(mlm_test$coefficients)-1)

#### Testing random forest model ####

### Shortcoming of random forest
# Random Forests aren't good at generalizing cases with completely new data. For example, if I tell you that one ice-cream costs $1, 2 ice-creams cost $2, and 3 ice-creams cost $3, how much do 10 ice-creams cost? A linear regression can easily figure this out, while a Random Forest has no way of finding the answer.
# Random forests are biased towards the categorical variable having multiple levels (categories). It is because feature selection based on impurity reduction is biased towards preferring variables with more categories so variable selection (importance) is not accurate for this type of data.

library(randomForest)
set.seed(123)

## actual model
RF_test <- randomForest(Y ~ ., data = X, mtry = 9, ntree = 1500)
plot(RF_test)
varImpPlot(RF_test)
importance(RF_test)

which.min(RF_test$mse)

## tuning model
set.seed(123)
tuneRF(
  x=X, #define predictor variables
  y=Y, #define response variable
  ntreeTry=500,
  mtryStart=4, 
  stepFactor=1.5,
  improve=0.01,
  trace=TRUE #don't show real-time progress
)

# take suggested mtry
set.seed(123)
RF_test_tuned <- randomForest(Y ~ ., data = X, mtry = 4, ntree = 1500)

## evaluation
RF_eval <-  evaluate.model(RF_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "RF", CpGs = length(mlm_test$coefficients)-1)

RF_eval_tuned <-  evaluate.model(RF_test_tuned, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "RF (tuned)", CpGs = length(mlm_test$coefficients)-1)

RF_eval_plot <- RF_eval$plot_train + RF_eval$plot_test + RF_eval_tuned$plot_train + RF_eval_tuned$plot_test +
  plot_layout(nrow=2)

ggsave(filename = paste0("002_plots/005_m_RF_age_all", extension), RF_eval_plot, width = 10, height = 7)

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
SVM_eps_eval <-  evaluate.model(SVM_eps_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "SVM (eps, linear) prediction", CpGs = length(mlm_test$coefficients)-1)

SVM_nu_eval <-  evaluate.model(SVM_nu_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI= color_species, plot_title = "SVM (nu, linear) prediction", CpGs = length(mlm_test$coefficients)-1)

SVM_eval_plot <- SVM_eps_eval$plot_train + SVM_eps_eval$plot_test + SVM_nu_eval$plot_train + SVM_nu_eval$plot_test +
  plot_layout(nrow = 2)
# transformed
SVM_eps_eval_t <-  evaluate.model(SVM_eps_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "SVM (eps, linear) prediction", CpGs = length(mlm_test$coefficients)-1)

SVM_nu_eval_t <-  evaluate.model(SVM_nu_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI= color_species, plot_title = "SVM (nu, linear) prediction", CpGs = length(mlm_test$coefficients)-1)

SVM_eval_t_plot <- SVM_eps_eval_t$plot_train + SVM_eps_eval_t$plot_test + SVM_nu_eval_t$plot_train + SVM_nu_eval_t$plot_test +
  plot_layout(nrow = 2)

ggsave(filename = paste0("002_plots/005_m_SVM_age_all", extension), SVM_eval_plot, width = 10, height = 7)
ggsave(filename = paste0("002_plots/005_m_SVM_log-age_all", extension), SVM_eval_t_plot, width = 10, height = 7)

#### Testing Bayesian models ####
# install.packages("brms")
# install.packages("rstan")  # Required for brms
library(brms)
library(rstan)

# setting up formula for model
BM_formula <- bf(rel_age ~.)
BM_model <- brm(formula = BM_formula, data = trainingData, family = gaussian(), chains = 4, cores = min(10, parallel::detectCores()), iter = 2000)
# summary(BM_model)
# plot(BM_model)

BM_model_t <- brm(formula = BM_formula, data = trainingData_t, family = gaussian(), chains = 4, cores = min(10, parallel::detectCores()), iter = 2000)

sel_cols <- colnames(trainingData[,-ncol(trainingData)])
# plot(conditional_effects(BM_model, effects = sel_cols))

pp_check(BM_model)

BM_pp <- posterior_predict(BM_model) 

ci <- posterior_interval(BM_pp, prob = 0.95)

# predictions
# normal
BM_predict_train <- predict(BM_model, trainingData[-length(trainingData)])
BM_predict_test <- predict(BM_model, testingData[-length(testingData)])

#transformed
BM_predict_train_t <- predict(BM_model_t, trainingData[-length(trainingData)])
BM_predict_test_t <- predict(BM_model_t, testingData[-length(testingData)])

# prediction data frame
# normal
prediction_data_train <- data.frame(Observed = trainingData$rel_age, Predicted = BM_predict_train[,1], Lwr = BM_predict_train[,3], Upr = BM_predict_train[,4])
prediction_data_test <- data.frame(Observed = testingData$rel_age, Predicted = BM_predict_test[,1], Lwr = BM_predict_test[,3], Upr = BM_predict_test[,4])
# transformed
prediction_data_train_t <- data.frame(Observed = trainingData$rel_age, Predicted = BM_predict_train_t[,1], Lwr = BM_predict_train[,3], Upr = BM_predict_train[,4])
prediction_data_test_t <- data.frame(Observed = testingData$rel_age, Predicted = BM_predict_test_t[,1], Lwr = BM_predict_test[,3], Upr = BM_predict_test[,4])

prediction_data_train_t$Predicted <- exp(-exp(-prediction_data_train_t$Predicted))
prediction_data_test_t$Predicted <- exp(-exp(-prediction_data_test_t$Predicted))

# metrics
# normal
metrics_BM_train <- data.frame(
  R = round(cor(prediction_data_train$Predicted, prediction_data_train$Observed, method = "pearson"), 4),
  MSE = round(mean((prediction_data_train$Predicted - prediction_data_train$Observed)^2), 4),
  MAE = round(mean(abs(prediction_data_train$Predicted - prediction_data_train$Observed)), 4),
  N = nrow(prediction_data_train))

metrics_BM_test <- data.frame(
  R = round(cor(prediction_data_test$Predicted, prediction_data_test$Observed, method = "pearson"), 4),
  MSE = round(mean((prediction_data_test$Predicted - prediction_data_test$Observed)^2), 4),
  MAE = round(mean(abs(prediction_data_test$Predicted - prediction_data_test$Observed)), 4),
  N = nrow(prediction_data_test))

# transformed
metrics_BM_train_t <- data.frame(
  R = round(cor(prediction_data_train_t$Predicted, prediction_data_train_t$Observed, method = "pearson"), 4),
  MSE = round(mean((prediction_data_train_t$Predicted - prediction_data_train_t$Observed)^2), 4),
  MAE = round(mean(abs(prediction_data_train_t$Predicted - prediction_data_train_t$Observed)), 4),
  N = nrow(prediction_data_train_t))

metrics_BM_test_t <- data.frame(
  R = round(cor(prediction_data_test_t$Predicted, prediction_data_test_t$Observed, method = "pearson"), 4),
  MSE = round(mean((prediction_data_test_t$Predicted - prediction_data_test_t$Observed)^2), 4),
  MAE = round(mean(abs(prediction_data_test_t$Predicted - prediction_data_test_t$Observed)), 4),
  N = nrow(prediction_data_test_t))

# ggplot(prediction_data, aes(x = Observed, y = Predicted)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = Lwr, ymax = Upr), width = 0.2, alpha = 0.2) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
#   theme_minimal() +
#   ylim(0,0.25) +
#   xlim(0,0.25) +
#   labs(x = "Observed Age", y = "Predicted Age", title = "Bayesian Model Predictions with Uncertainty Intervals")
# 
# BM_eval <- evaluate.model(BM_model, trainingData[-length(trainingData)], trainingData$rel_age, testingData[-length(testingData)], testingData$rel_age, meth_train$species, meth_test$species,  CpGs = ncol(training_data)-1)
# 
# BM_eval_t <- evaluate.model(BM_model_t, trainingData[-length(trainingData)], trainingData$rel_age, testingData[-length(testingData)], testingData$rel_age, meth_train$species, meth_test$species,  CpGs = ncol(training_data)-1, transform = TRUE)

#### Testing deep learning models ####
# install.packages("keras")
# # library(reticulate)
# library(keras)
# DL_test <- keras_model_sequential() %>%
#   layer_dense(units = 64, activation = 'relu', input_shape = c(num_features)) %>%
#   layer_dense(units = 1)
# DL_test %>% compile(
#   optimizer = 'rmsprop',
#   loss = 'mse',
#   metrics = 'mae'
# )
# DL_test %>% fit(X, Y, epochs = 100, batch_size = 10, validation_split = 0.2)

#### Summarizing all models ####
# training data 
df_all_training <- rbind(GLM_eval$metrics_train, mlm_eval$metrics_train, MLM_eval$metrics_train, metrics_BM_train, SVM_eps_eval$metrics_train, SVM_nu_eval$metrics_train, RF_eval$metrics_train, RF_eval_tuned$metrics_train, GLM_eval_t$metrics_train, mlm_eval_t$metrics_train, MLM_t_eval$metrics_train, metrics_BM_train_t, metrics_LOOCV, metrics_LOOCV_t)
rownames(df_all_training) <- c("GLM", "mlm", "MLM (caret)", "BM", "SVM_eps", "SVM_nu", "RF", "RF_tuned", "GLM_t", "mlm_t", "MLM_t", "BM_t", "MLM_LOOCV", "MLM_LOOCV_t")
df_all_training

df_all_testing <- rbind(GLM_eval$metrics_test, mlm_eval$metrics_test, MLM_eval$metrics_test, metrics_BM_test, SVM_eps_eval$metrics_test, SVM_nu_eval$metrics_test, RF_eval$metrics_test, RF_eval_tuned$metrics_test, GLM_eval_t$metrics_test, mlm_eval_t$metrics_test, MLM_t_eval$metrics_test, metrics_BM_test_t, metrics_LOOCV, metrics_LOOCV_t)
rownames(df_all_testing) <- c("GLM", "mlm", "MLM (caret)", "BM", "SVM_eps", "SVM_nu", "RF", "RF_tuned", "GLM_t", "mlm_t", "MLM_t", "BM_t", "MLM_LOOCV", "MLM_LOOCV_t")
df_all_testing

df_eval <- cbind(df_all_training, df_all_testing)

df_eval

#### plotting ####
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>%
  .[c(-1,-9)]
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

