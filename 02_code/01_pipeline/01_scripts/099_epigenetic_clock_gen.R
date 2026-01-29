### Showcasing epigenetic age models 
### generating data
set.seed(1999)

age <- runif(267, min = 3, max = 87)

# Function to generate estimated_age with target correlation
generate_estimated_age <- function(age, target_cor = 0.9) {
  age_std <- scale(age)  # standardize to mean 0, sd 1
  noise <- rnorm(length(age), mean = 0, sd = (0.006 * age))
  noise2 <- rnorm(length(age), mean = 0, sd = abs(rnorm(267, mean = 0.5, sd = 1)))
  estimated_age_std <- age_std + noise  # combine age signal + noise
  estimated_age <- as.numeric(scale(estimated_age_std)) * sd(age) + mean(age)  # rescale to match age distribution
  return(estimated_age)
}

estimated_age <- generate_estimated_age(age, target_cor = 0.8)

# Create dataframe
data_gen <- data.frame(age = age, estimated_age = estimated_age)

# Check correlation
data_gen_cor <- round(cor(data_gen$age, data_gen$estimated_age), 4)  # should be ~0.9


ggplot(data_gen, aes(x = age, y = estimated_age)) +
  geom_point(size = 2, alpha = .9, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_smooth(method = "lm", se = TRUE, color = "orange", fill = "grey90", fullrange = TRUE) +
  # geom_line(data = result_df_train, aes(x = age, y = age_predicted), color = "black", alpha = 0.7, linetype = "dashed") +
  labs(y = "Estimated age (years)", x = "Chronological age (years)", color = "Species") +
  # labs(subtitle = paste0("R=", metrics_test$R, "\nMSE=", metrics_test$MSE, "\nMAE=", metrics_test$MAE, "\nN=", nrow(X_test), " CpGs=", CpGs)) +
  theme_bw() +
  xlim(c(0, 100)) +
  ylim(c(0, 100)) +
  # annotate("text", x = 0, y = 10, label = paste0("R=", data_gen_cor), size = 3.5, hjust = 0, vjust = -35) +
  theme(plot.title = element_text(hjust = .5, face = "bold"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_text(face = "bold", size = 14))

ggsave(filename = "002_plots/999_epigenetic_clock.pdf", width = 5.5, height = 5)
 