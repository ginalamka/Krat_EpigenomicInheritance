# Analyze bawmeth output for DMR analysis

setwd("/home/gfl0003/kratepi/work/w_500_250")
setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/krat_epi/krat_epi/krat_epi_analysis/bwameth/Output") 
directory = getwd()
data <- read.table("methylation_imputed_matrix_82.csv", header = TRUE, sep = ",", row.names = 1) #filt_data_82_NA # NOTE this is without the first 5 columns!!
meta <- read.table("kratepi_metadata_full_numeric.csv", header = TRUE, sep=",")

# impute missing values
# https://people.math.ethz.ch/~maechler/adv_topics_compstat/softImpute_MM.html
{
library(softImpute)
anyNA(data)
dat <- as.matrix(data)

# find the best rank.max by finding the "elbow" from the non-NA dataset (must be matrix)
svd_vals <- svd(na.omit(dat)) 
plot(svd_vals$d, type = "b")
# suggests rank.max = 2

# find the best lambda
# first, create a mask for missing values
missing_mask <- is.na(dat)
# simulate missing values for cross-validation
test_mask <- matrix(runif(length(dat)) < 0.1, nrow = nrow(dat))
test_mask[missing_mask] <- FALSE  # Keep original NAs untouched
train_data <- dat
train_data[test_mask] <- NA  # Hide some values for validation

# define range of lambda values
lambda_grid <- seq(0.001, 1, length.out = 10)
rmse_values <- numeric(length(lambda_grid))

# Perform cross-validation
for (i in seq_along(lambda_grid)) {
  lambda <- lambda_grid[i]
  imputed_model <- softImpute(train_data, lambda = lambda, trace = FALSE)
  imputed_data <- complete(train_data, imputed_model)
  
  # Compute RMSE on held-out test data
  rmse_values[i] <- sqrt(mean((imputed_data[test_mask] - dat[test_mask])^2, na.rm = TRUE))
}

# Find the best lambda
best_lambda <- lambda_grid[which.min(rmse_values)]
plot(lambda_grid, rmse_values, type = "b", xlab = "Lambda", ylab = "RMSE")
# The lambda with the lowest RMSE is the best choice.

# Fit softImpute models for lambda = 0.778 and lambda = 0.7
model_778 <- softImpute(train_data, lambda = 0.778, trace = FALSE)
model_700 <- softImpute(train_data, lambda = 0.7, trace = FALSE)

# Get imputed values
imputed_778 <- complete(train_data, model_778)
imputed_700 <- complete(train_data, model_700)

# Compute RMSE for both models
rmse_778 <- sqrt(mean((imputed_778[test_mask] - dat[test_mask])^2, na.rm = TRUE))
rmse_700 <- sqrt(mean((imputed_700[test_mask] - dat[test_mask])^2, na.rm = TRUE))

# Print RMSE results
print(c(RMSE_lambda_778 = rmse_778, RMSE_lambda_700 = rmse_700))
# RMSE_lambda_778 RMSE_lambda_700 
# 11.46917        11.63905 

par(mfrow = c(1, 2))  # Side-by-side plots
hist(imputed_778, breaks = 50, main = "Lambda = 0.778", xlab = "Imputed Values", col = "blue", xlim = c(0,100))
hist(imputed_700, breaks = 50, main = "Lambda = 0.7", xlab = "Imputed Values", col = "red", xlim = c(0,100))

# try the suggested combos, looking for a range between 0 - 100
h <- softImpute(dat, rank.max = 1, lambda = 0.7, trace = TRUE, type = "als") 
j <- softImpute(dat, rank.max = 1, lambda = 0.778, trace = TRUE, type = "als") # best_lambda = 0.778
h$d
# up1 <- complete(dat, i) #i, rank.max = 2, lamba = 0.7
up2 <- complete(dat, h) #h, rank.max = 1, lamba = 0.7
up3 <- complete(dat, j) #j, rank.max = 1, lamba = 0.778
range(up3)
summary(up3)
# lamba = 0.778 slightly overinflates estimates
}
# best model has lambda = 0.778 and rank.max = 1
# up next: use up3 for ridge regressions!
imp_data <- up3

# double check imputed and orig values are similar
hist(dat)
hist(imp_data)

rm(imputed_700, imputed_778, imputed_data, imputed_model, j, missing_mask, 
   model_700, model_778, svd_vals, test_mask, train_data, up2, up3, best_lambda, 
   h, i, lambda, lambda_grid, rmse_700, rmse_778, rmse_values) #clean up

# this is the ChatGPT convo that got me here : https://chatgpt.com/share/67a3a83c-3e20-8012-b543-0def502b8fef

# at this point, 
#.    dat (matrix) = data (df) = orig csv with NAs
#.    meta = prepped metadata
#.    imp_data = data with imputed values with rank.max = 1 and lambda = 0.778

# Transpose the data frame
transposed_data <- as.matrix(t(imp_data))
dataset <- cbind(meta, transposed_data)
rownames(dataset) <- rownames(transposed_data)
meta <- dataset[,1:8]
data <- dataset[,9:ncol(dataset)]
rm(imp_data) #, transposed_data

write.csv(transposed_data, file = paste(directory, "methylation_imputed_matrix_82.csv"), sep = ",", row.names = TRUE, quote = FALSE)
write.csv(dataset, file = paste(directory, "dataset_imputed_df_82.csv"), sep = ",", row.names = TRUE, quote = FALSE)

dataset <- read.table("dataset_imputed_df_82.csv", header = TRUE, sep = ",", row.names = 1)
transposed_data <- as.matrix(read.table("methylation_imputed_matrix_82.csv", header = TRUE, sep = ",", row.names = 1))

data <- dataset[,9:ncol(dataset)]
meta <- dataset[,1:8]
}
# let's look at averages
{
avg_indv = matrix(nrow = 48, ncol = 1)
avg_wind = matrix(nrow = ncol(data), ncol = 1)
for(i in 1:ncol(data)){
  avg_wind[i,] <- mean(data[1:nrow(data),i])
}
for(j in 1:nrow(data)){
  avg_indv[j,] <- sum(data[j,])/ncol(data)   #mean(data[j, 1:ncol(data)])
}
summary(avg_wind)
summary(avg_indv)
plot(avg_indv)
plot(avg_wind)
}

## added 6/3/25 -- compare the imputed and missing datasets
impu = read.table("methylation_imputed_matrix_82.csv", header = TRUE, sep = ",") # , row.names = 1
miss = read.table("filt_data_82_NA.csv", header = TRUE, sep = ",")
trans_miss <- as.data.frame(t(miss))
trans_miss$X <- impu$X
# trans_miss$grp <- "missing"
# impu$grp <- "imputed"

## added 6/3/2026 -- calculate the amount of missing data before imputation
sum(is.na(miss)) # To find the absolute total number of NA values across the entire dataset
#85144 / (48*14785) = 0.1199752

impu_mean = NULL
for(i in 2:ncol(impu)){
  m <- as.numeric(mean(impu[,i]))
  impu_mean = c(impu_mean, m)
}

miss_mean = NULL
for(i in 1:(ncol(trans_miss)-1)){ # this is cuz the sample names are at the last column, not first
  mm <- as.numeric(mean(trans_miss[,i], na.rm = TRUE))
  miss_mean = c(miss_mean, mm)
}

imputeddata <- cbind(impu_mean, 1, 1:nrow(imputeddata))
missingnessdata <- cbind(miss_mean, 0, 1:nrow(imputeddata))
compare <- rbind(imputeddata, missingnessdata)
colnames(compare) <- c("mean", "group", "window")
colnames(imputeddata) <- c("mean", "group", "window")
colnames(missingnessdata) <- c("mean", "group", "window")
compare <- as.data.frame(compare)
compare$group <- as.factor(compare$group)

missingnessdata = as.data.frame(missingnessdata)
lm_miss <- lm(window ~ mean, data = na.omit(missingnessdata))
summary(lm_miss)

imputeddata = as.data.frame(imputeddata)
lm_imput <- lm(window ~ mean, data = imputeddata)
summary(lm_imput)

library(ggplot2)

# Convert group to a factor for better labeling
compare$group <- factor(compare$group, labels = c("Original", "Imputed"))

# Plot the densities
ggplot(compare, aes(x = mean, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 1) +
  labs(
    title = "Density Plot of Original vs. Imputed Methylation Values",
    x = "Mean Methylation (%)",
    y = "Density",
    fill = "Group",
    color = "Group"
  ) +
  theme_classic()

# Check if distributions are sig different
wilcox.test(mean ~ group, 
            data = compare, 
            exact = FALSE)

