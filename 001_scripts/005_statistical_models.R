#### Overview ####
# testing different statistical models 

#### Settings ####
setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)
library(caret)

#### Loading data ####

load("/workspace/cfngle/results-data/06_model_creation/all_meth_values_selected.RData")

#### Data splitting ####

# option A):
meth_train <- all_meth_values_selected[-seq(1, nrow(all_meth_values_selected), 4),]
meth_test <- all_meth_values_selected[seq(1, nrow(all_meth_values_selected), 4),]

# save(all_meth_values_selected, file = "/workspace/cfngle/results-data/06_model_creation/HS_all_meth_values_selected.RData")

# option B):
# using a stratified splitting technique
set.seed(123)
AC_split <- initial_split(AC_meth_values_selected, strata = "rel_age")

set.seed(123)
AS_split <- initial_split(AS_meth_values_selected, strata = "rel_age")

set.seed(123)
EH_split <- initial_split(EH_meth_values_selected, strata = "rel_age")

set.seed(123)
ZF_split <- initial_split(ZF_meth_values_selected, strata = "rel_age")

meth_train <- rbind(training(AC_split),
                    training(AS_split),
                    training(EH_split),
                    training(ZF_split))

meth_test <- rbind(testing(AC_split),
                   testing(AS_split),
                   testing(EH_split),
                   testing(ZF_split))

# checking how many CpGs are present per data set
nrow(meth_train)
nrow(meth_test)

## plotting age distribution
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>%
  .[c(-1,-9)]

plot_age_dist <- ggplot(all_meth_values_selected) +
  geom_density(aes(x = rel_age, color = species), linewidth = 2) +
  geom_density(aes(x = rel_age, fill = "all"), alpha = 0.2, size = 0) +
  scale_color_manual(values = colpalOI) +
  scale_fill_manual(values = colpalOI[5]) +
  theme_minimal() +
  labs(title = "Relative age distribution for all samples")

# >> the plot shows that our dependent variable is not normally distributed (as expected from age) 

plot_sample_age_dist <- ggplot() +
  geom_density(data = meth_train, aes(x = rel_age, fill = "Training"), alpha = 0.5, size = 0) +
  geom_density(data = meth_test, aes(x = rel_age, fill = "Testing"), alpha = 0.5, size = 0) +
  scale_fill_manual(values = colpalOI) +
  labs(fill = "Dataset") +
  theme_minimal() +
  ggtitle("Age Distribution in Training and Testing Sets")

plot(density(EH_meth_values_selected$rel_age))

# plotting both graphs 
plot_age_dist + plot_sample_age_dist +
  plot_layout(nrow=1)

### defining training data
X <- meth_train %>% 
  select(-rel_age, -species)

Y <- meth_train[,"rel_age"]

# testing data
X_test <- meth_test %>% 
  select(-rel_age, -species)

Y_test <- meth_test[,"rel_age"]

#### model evaluation ####

## function for model testing plus calculating metrics (MSE, MAE, R)
evaluate.model <- function(model, X_train, Y_train, X_test, Y_test, species_train,
                           species_test, transform = FALSE, colpalOI, plot_title = "Model evaluation:", 
                           y_lim = c(0,.3), x_lim = c(0,.3), CpGs = "not defined") {
  # Calculate predictions
  predictions_train <- predict(model, X_train)
  predictions_test <- predict(model, X_test)
  
  # Optional transformation
  if (transform) {
    predictions_train <- exp(-exp(-predictions_train))
    predictions_test <- exp(-exp(-predictions_test))
    Y_train
    Y_test
  }
  
  # Calculate metrics
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
  
  # Prepare data frames for plotting
  result_df_train <- data.frame(age_predicted = predictions_train, age = Y_train, species = species_train)
  result_df_test <- data.frame(age_predicted = predictions_test, age = Y_test, species = species_test)
  
  # Create plots
  plot_train <- ggplot(result_df_train, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_smooth(method = "lm", se = FALSE, color = "black") +
    # geom_line(data = result_df_train, aes(x = age, y = age_predicted), color = "black", alpha = 0.7, linetype = "dashed") +
    scale_color_manual(values = colpalOI) +
    ylim(y_lim) +
    xlim(x_lim) +
    labs(title = paste(plot_title, "(Training Set)"), y = "Estimated age (relative)", x = "Chronological age",
         subtitle = paste0("R=", metrics_train$R, " MSE=", metrics_train$MSE, " MAE=", metrics_train$MAE, " N=", nrow(X_train), " CpGs=", CpGs)) +
    theme_classic()
  
  plot_test <- ggplot(result_df_test, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    # geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_manual(values = colpalOI) +
    ylim(y_lim) +
    xlim(x_lim) +
    labs(title = paste(plot_title, "(Testing Set)"), y = "Estimated age (relative)", x = "Chronological age",
         subtitle = paste0("R=", metrics_test$R, " MSE=", metrics_test$MSE, " MAE=", metrics_test$MAE, " N=", nrow(X_test), " CpGs=", CpGs)) +
    theme_classic()
  
  # Return list containing metrics and plots
  return(list(metrics_train = metrics_train, metrics_test = metrics_test, plot_train = plot_train, plot_test = plot_test))
}

X <- X
#### Testing model CV GLM ####

# change age here!
# Y <- meth_train[,"rel_age"]
# Y_test <- meth_test[,"rel_age"]

#changing value
Y <- -log(-log(Y))
Y_test <- -log(-log(Y_test))

# define alpha
glm_alpha <- 0.5

# setting seed for reproducibility 
set.seed(123)

# model
GLM_test <- cv.glmnet(as.matrix(X), Y, alpha = glm_alpha)
plot(GLM_test)

coef(GLM_test, s=GLM_test$lambda.min)

### running model on testing data
### plotting
GLM_predictions_test <- predict(GLM_test, as.matrix(X_test), s=GLM_test$lambda.min) %>% 
  as.vector()
GLM_predictions_cor_test <- cor(GLM_predictions_test, Y_test, method = "pearson") %>% 
  round(4)
GLM_mse_test <- mean((GLM_predictions_test - Y_test)^2) %>% 
  round(4) 
GLM_mae_test <- mean(abs(GLM_predictions_test - Y_test)) %>% 
  round(4)

# for normal Y
GLM_result_df <- data.frame(predictions = GLM_predictions_test,
                            rel_age = Y_test,
                            species = "ZF")
# for transformed Y!

GLM_result_df <- data.frame(predictions = exp(-exp(-GLM_predictions_test)),
                            rel_age = exp(-exp(-Y_test)),
                            species = "ZF")
# for transformed Y!
GLM_mse_test <- mean((exp(-exp(-GLM_predictions_test)) - exp(-exp(-Y_test)))^2) %>% 
  round(4) 
GLM_mae_test <- mean(abs(exp(-exp(-GLM_predictions_test)) - exp(-exp(-Y_test)))) %>% 
  round(4)


colnames(GLM_result_df) <- c("age_predicted", "age", "species")

GLM_plot_test <- ggplot(GLM_result_df, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = colpalOI[4]) +
  ylab("Estimated age (relative)") +
  xlab("Chronological age") +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction GLM (test data set)",
       subtitle = paste0("R=", GLM_predictions_cor_test,
                         " MSE=", GLM_mse_test,
                         " MAE=", GLM_mae_test,
                         " N=", nrow(X_test), 
                         " alpha=", glm_alpha)) +
  theme_classic()

### running model on training data
GLM_predictions_train <- predict(GLM_test, newx= as.matrix(X), s=GLM_test$lambda.min)
GLM_predictions_cor_train <- cor(GLM_predictions_train, Y, method = "pearson") %>% 
  round(4)
GLM_mse_train <- mean((GLM_predictions_train - Y)^2) %>% 
  round(4)
GLM_mae_train <- mean(abs(GLM_predictions_train - Y)) %>% 
  round(4)

# for normal Y
GLM_result_df_train <- data.frame(predictions = GLM_predictions_train,
                                  rel_age = Y,
                                  species = all_meth_values_selected[all_meth_values_selected$species != "ZF",]$species)

# for transformed Y!
GLM_result_df_train <- data.frame(predictions = exp(-exp(-GLM_predictions_train)),
                                  rel_age = exp(-exp(-Y)),
                                  species = all_meth_values_selected[all_meth_values_selected$species != "ZF",]$species)
# for transformed Y!
GLM_mse_train <- mean((exp(-exp(-GLM_predictions_train)) - exp(-exp(-Y)))^2) %>% 
  round(4) 
GLM_mae_train <- mean(abs(exp(-exp(-GLM_predictions_train)) - exp(-exp(-Y)))) %>% 
  round(4)

colnames(GLM_result_df_train) <- c("age_predicted", "age", "species")

GLM_plot_train <- ggplot(GLM_result_df_train, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted), size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = colpalOI) +
  ylab("Estimated age (relative)") +
  xlab("Chronological age") +
  ylim(0,0.30) +
  xlim(0,0.30) +
  labs(title = "Multispecies age (relative, -log(-log(x))) prediction GLM (train data set)",
       subtitle = paste0("R=", GLM_predictions_cor_train, " MSE=", GLM_mse_train, " MAE=", GLM_mae_train, " N=", nrow(X))) +
  theme_classic()

GLM_plot_train + GLM_plot_test +
  plot_layout(nrow = 1)


#### Testing multivariate linear regression models ####


### pre-testing with base R package
mlm_test <- lm(Y ~., data = X)
summary(mlm_test)

# with transformed age
mlm_test_t <- lm(-log(-log(Y)) ~., data = X)

## selecting only significant values
sign_vec <- as.vector((summary(mlm_test)$coefficients[,4] < 0.05)[-1])
# stat tests 
plot(density(mlm_test$residuals))
shapiro.test(mlm_test$residuals)

## selecting only significant values
mlm_test_opt <- lm(Y ~., data = X[,sign_vec])
summary(mlm_test_significant)
# with transformed age
mlm_test_opt_t <- lm(-log(-log(Y)) ~.,X[sign_vec])

## using centered and scaled dataset 
# did not change anything!
mlm_alt <- lm(Y ~., trainingData)
predict(mlm_alt, testingData)

### testing and plotting model
# color palettes
colpal_CB <- c("#c06d00", "#f9cf6e", "#6a5d00", "#44a02b", "#008649", "#12ebf0", "#65a9ff", "#004588", "#660077", "#ff98f7", "#954674", "#630041")
colpal_CB_01 <- colpal_CB[c(TRUE, FALSE)]
colpal_CB_02 <- colpal_CB[c(FALSE, TRUE)]

colpal_CB_a <- c("#f8cbb1","#006786","#182057","#6b6300","#ff8ab9","#f1aaff","#bb005a","#013aa8","#01ef9a","#fa8200","#ee0028","#26c100")
colpal_CB_a_01 <- colpal_CB_a[1:6]
colpal_CB_a_02 <- colpal_CB_a[7:12]

# testing normal mlm
mlm_eval <-  evaluate.model(mlm_test, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, 
                            colpalOI= colpal_CB, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)
mlm_alt_eval <-  evaluate.model(mlm_alt, trainingData, Y, testingData, Y_test, meth_train$species, meth_test$species, transform = FALSE, 
                                colpalOI= colpal_CB, plot_title = "MLM prediction", CpGs = length(mlm_test$coefficients)-1)
mlm_eval_t <-  evaluate.model(mlm_test_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, 
                              colpalOI= colpal_CB, plot_title = "MLM (age -log-log transformed) prediction", CpGs = length(mlm_test_t$coefficients)-1)

# testing optimized mlm
mlm_eval_opt <-
  evaluate.model(mlm_test_opt, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = FALSE, colpalOI = colpal_CB_02, plot_title = "MLM (sign. CpGs only) prediction", CpGs = length(mlm_test_opt$coefficients) - 1
  )
mlm_eval_opt_t <-
  evaluate.model(mlm_test_opt_t, X, Y, X_test, Y_test, meth_train$species, meth_test$species, transform = TRUE, colpalOI = colpal_CB_02, plot_title = "MLM (sign. CpGs only; age -log-log transformed) prediction", CpGs = length(mlm_test_opt$coefficients) - 1
  )

## plotting 
# normal age
mlm_eval$plot_train + mlm_eval$plot_test  + mlm_eval_opt$plot_train + mlm_eval_opt$plot_test +
  plot_layout(nrow = 2)
# transformed age
mlm_eval_t$plot_train + mlm_eval_t$plot_test  + mlm_eval_opt_t$plot_train + mlm_eval_opt_t$plot_test +
  plot_layout(nrow = 2)

#### MLM (caret) ####

## defining data
# training data
X <- meth_train %>% 
  select(-rel_age, -species)

Y <- meth_train[,"rel_age"]

# testing data
X_test <- meth_test %>% 
  select(-rel_age, -species)

Y_test <- meth_test[,"rel_age"]

## pre processing
# centering means that all the the mean is being substracted from all values and scale divides them by the standard deviation. 

# training data
preProcValues <- preProcess(X, method = c("center", "scale"))
trainingData <- predict(preProcValues, X)
trainingData$rel_age <- Y
# testing data
preProcValues_test <- preProcess(X_test, method = c("center", "scale"))
testingData <- predict(preProcValues_test, X_test)
testingData$rel_age <- Y_test

# training model 
set.seed(123)
trainControl <- trainControl(method = "cv", number = 10) # 10-fold CV
MLM_model <- train(rel_age ~ ., data = trainingData, 
                   method = "glmnet", 
                   trControl = trainControl)
# prediction test
MLM_predict_test <- predict(MLM_model, newdata = X_test)
# caret has functions to interpret prediction results
performanceResults <- postResample(MLM_predict_test, Y_test)
print(performanceResults)

#tuning model 
MLM_tuned <- train(rel_age ~ ., data = testingData, 
                   method = "glmnet",
                   tuneLength = 10, # Number of tuning parameter values
                   trControl = trainControl)

importance <- varImp(MLM_tuned, scale = FALSE)
plot(importance)


#### Testing random forest model ####

### Shortcoming of random forest
# Random Forests aren't good at generalizing cases with completely new data. For example, if I tell you that one ice-cream costs $1, 2 ice-creams cost $2, and 3 ice-creams cost $3, how much do 10 ice-creams cost? A linear regression can easily figure this out, while a Random Forest has no way of finding the answer.
# Random forests are biased towards the categorical variable having multiple levels (categories). It is because feature selection based on impurity reduction is biased towards preferring variables with more categories so variable selection (importance) is not accurate for this type of data.

library(randomForest)
set.seed(123)
# change age here!
Y <- meth_train[,"rel_age"]
Y_test <- meth_test[,"rel_age"]
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
                           species = meth_test$species)
# for transformed Y!

RF_result_df <- data.frame(predictions = exp(-exp(-RF_predictions_test)),
                           rel_age = exp(-exp(-Y_test)),
                           species = meth_test$species)
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
                                 species = meth_train$species)

# for transformed Y!
RF_result_df_train <- data.frame(predictions = exp(-exp(-RF_predictions_train)),
                                 rel_age = exp(-exp(-Y)),
                                 species = meth_train$species)
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
Y <- meth_train[,"rel_age"]
Y_test <- meth_test[,"rel_age"]
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
                            species = meth_test$species)
# for transformed Y!

SVM_result_df <- data.frame(predictions = exp(-exp(-SVM_predictions_test)),
                            rel_age = exp(-exp(-Y_test)),
                            species = meth_test$species)
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
                                  species = meth_train$species)

# for transformed Y!
SVM_result_df_train <- data.frame(predictions = exp(-exp(-SVM_predictions_train)),
                                  rel_age = exp(-exp(-Y)),
                                  species = meth_train$species)
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
install.packages("brms")
install.packages("rstan")  # Required for brms
library(brms)
library(rstan)

# setting up formula for model
BM_formula <- bf(rel_age ~.)
BM_model <- brm(formula = BM_formula, data = trainingData, family = gaussian(), chains = 4, cores = min(10, parallel::detectCores()), iter = 2000)
summary(BM_model)
plot(BM_model)


sel_cols <- colnames(trainingData[,-ncol(trainingData)])
plot(conditional_effects(BM_model, effects = sel_cols))

pp_check(BM_model)

BM_pp <- posterior_predict(BM_model) 

BM_predict_test <- predict(BM_model, testingData)


ci <- posterior_interval(BM_pp, prob = 0.95)

# Merge means and intervals
prediction_data <- data.frame(Observed = testingData$rel_age, Predicted = BM_predict_test[,1], Lwr = BM_predict_test[,3], Upr = BM_predict_test[,4])

# Use ggplot2 for plotting
library(ggplot2)
ggplot(prediction_data, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lwr, ymax = Upr), width = 0.2, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  ylim(0,0.25) +
  xlim(0,0.25) +
  labs(x = "Observed Age", y = "Predicted Age", title = "Bayesian Model Predictions with Uncertainty Intervals")

evaluate.model(BM_model, trainingData[-length(trainingData)], trainingData$rel_age, testingData[-length(testingData)], testingData$rel_age, meth_train$species, meth_test$species,  CpGs = ncol(training_data)-1)

#### Testing deep learning models ####
install.packages("keras")
# library(reticulate)
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