# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: June 2024

#### Overview ####
# Testing LOOCV for ENR, MLM and RF

#### Testing ####

### ENR
# leave one out 

fit_results <- c()
  for(i in 1:nrow(X)){
    # splitting data into train and test
    X_train <- X[-i,]
    Y_train <- Y[-i]
    X_test <- X[i,]
  
    fit <- cv.glmnet(as.matrix(X_train), Y_train, alpha = glm_alpha, family = "gaussian")
    Y_fitted <- predict(fit, as.matrix(X_test), s = fit$lambda.min)
    fit_results[i] <- Y_fitted
}
prediction <- data.frame(age_predicted = fit_results, age = Y)
return(prediction)

ggplot(prediction, aes(x = age, y = age_predicted)) +
  geom_point()


# leave one sample out 
all_values <- all_meth_values_selected[,-c(length(all_meth_values_selected),length(all_meth_values_selected)-1)]
all_age <- all_meth_values_selected$rel_age

# Automatically adding second-degree terms for all numeric predictors in X_train
X_train <- all_values
numeric_predictors <- sapply(X_train, is.numeric)
X_train_poly <- as.data.frame(lapply(X_train[, numeric_predictors], function(x) x^2))
names(X_train_poly) <- paste0(names(X_train_poly), "_2")

# Bind the new polynomial terms back to the original X_train
X_train <- cbind(X_train, X_train_poly)

fit_results <- c()
for(i in 1:nrow(all_values)){ 
  # splitting data into train and test
  X_train <- X_train[-i,]
  Y_train <- all_age[-i]
  X_test <- data.frame(X_train[i,,drop = FALSE])
  
  fit <- lm(Y_train ~., data = X_train)
  Y_fitted <- predict(fit, newdata = X_test)
  fit_results[i] <- Y_fitted
}

prediction <- data.frame(age_predicted = fit_results, age = all_age)
# prediction <- data.frame(age_predicted = Y_fitted, age = Y_train)

metrics <- data.frame(
  R = round(cor(prediction$age_predicted, prediction$age, method = "pearson"), 4),
  MSE = round(mean((prediction$age_predicted - prediction$age)^2), 4),
  MAE = round(mean(abs(prediction$age_predicted - prediction$age)), 4))
metrics


### MLM

## funciton for LOOCV
test.loocv.lm <- function(X, Y, transformation = NULL, inv_transformation) {
  fit_results <- numeric(nrow(X)) # pre-allocate for better performance
  
  for(i in 1:nrow(X)){ 
    # splitting data into train and test
    X_train <- X[-i,]
    Y_train <- if(!is.null(transformation)) {
      sapply(Y[-i], transformation) # apply transformation if not NULL
    } else {
      Y[-i]
    }
    
    X_test <- data.frame(X[i,, drop = FALSE])
    
    fit <- lm(Y_train ~ ., data = X_train)
    Y_fitted <- predict(fit, newdata = X_test)
    fit_results[i] <- Y_fitted
  }
  
  if(!is.null(transformation)) {
    fit_results <- sapply(fit_results, inv_transformation) # apply inverse transformation if needed
  }
  
  prediction <- data.frame(age_predicted = fit_results, age = Y)
  
  metrics <- data.frame(
    R = round(cor(prediction$age_predicted, prediction$age, method = "pearson"), 4),
    MSE = round(mean((prediction$age_predicted - prediction$age)^2), 4),
    MAE = round(mean(abs(prediction$age_predicted - prediction$age)), 4))
  
  plot <- ggplot(prediction, aes(x = age, y = age_predicted)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    ylim(0,0.30) +
    xlim(0,0.30)
  
  return(list(metrics = metrics, plot = plot))
}


test.loocv.lm(all_values, all_age)
test.loocv.lm(all_values, all_age, function(x) -log(x), function(x) exp(-x))
test.loocv.lm(all_values, all_age, function(x) -log(-log(x)), function(x) exp(-exp(-x)))

test.loocv.lm(log(all_values), all_age)
test.loocv.lm(log(all_values), all_age, function(x) -log(x), function(x) exp(-x))
test.loocv.lm(log(all_values), all_age, function(x) -log(-log(x)), function(x) exp(-exp(-x)))

test.loocv.lm(poly(all_values, degree = 2), all_age)


ggplot(prediction, aes(x = age, y = age_predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  ylim(0,0.30) +
  xlim(0,0.30)

### RF

## funciton for LOOCV RF
test.loocv.rf <- function(X, Y, transformation = NULL, inv_transformation) {
  fit_results <- numeric(nrow(X)) # pre-allocate for better performance
  
  for(i in 1:nrow(X)){ 
    # splitting data into train and test
    X_train <- X[-i,]
    Y_train <- if(!is.null(transformation)) {
      sapply(Y[-i], transformation) # apply transformation if not NULL
    } else {
      Y[-i]
    }
    
    X_test <- data.frame(X[i,, drop = FALSE])
    
    # model
    fit <- randomForest(Y_train ~ ., data = X_train, mtry = 4, ntree = 1500)
    Y_fitted <- predict(fit, newdata = X_test)
    fit_results[i] <- Y_fitted
  }
  # apply inverse transformation if needed
  if(!is.null(transformation)) {
    fit_results <- sapply(fit_results, inv_transformation) 
  }
  
  prediction <- data.frame(age_predicted = fit_results, age = Y)
  
  metrics <- data.frame(
    R = round(cor(prediction$age_predicted, prediction$age, method = "pearson"), 4),
    MSE = round(mean((prediction$age_predicted - prediction$age)^2), 4),
    MAE = round(mean(abs(prediction$age_predicted - prediction$age)), 4))
  
  plot <- ggplot(prediction, aes(x = age, y = age_predicted)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    ylim(0,0.30) +
    xlim(0,0.30)
  
  return(list(metrics = metrics, plot = plot))
}

# test.loocv.rf(all_values, all_age)
