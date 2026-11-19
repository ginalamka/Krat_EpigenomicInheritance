# Scaffolding approach for methylation analysis, started 3/4/25
setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/krat_epi/krat_epi/krat_epi_analysis/bwameth/Output") 
directory = getwd()
data = read.table("DATASET.csv", header=T, sep=",")

# load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(gridExtra)
library(tidyverse)
library(VennDiagram)

# Step 1: Calculate methylation variation (SD/Mean) for each window
{
  methylation_vars <- data %>%
    select(year, 35:14818) %>%  # Select year and methylation data
    pivot_longer(cols = -year, names_to = "window", values_to = "methylation") %>%
    group_by(year, window) %>%
    summarise(variation = sd(methylation, na.rm = TRUE) / mean(methylation, na.rm = TRUE), .groups = "drop")
  
  # Scatterplot: Individual variation points per window
  ggplot(methylation_vars, aes(x = year, y = variation)) +
    geom_point(alpha = 0.3) +
    labs(title = "Methylation Variation by Year",
         x = "Year",
         y = "Methylation Variation (SD/Mean)") +
    theme_minimal()
}  

## STEP 2:
## AVERAGE HERITABILITY
{
  # Reshape methylation data to long format
  methylation_long <- data %>%
    select(lineage, 35:14818) %>%
    pivot_longer(cols = -lineage, names_to = "window", values_to = "methylation")
  
  # Fit a mixed-effects model for heritability/repeatability estimation
  model <- lmer(methylation ~ (1 | lineage), data = methylation_long, REML = TRUE)
  
  # Extract variance components
  var_components <- as.data.frame(VarCorr(model))
  var_lineage <- var_components$vcov[1]  # Lineage variance
  var_residual <- attr(VarCorr(model), "sc")^2  # Residual variance
  
  # Calculate repeatability (R)
  repeatability <- var_lineage / (var_lineage + var_residual)
  print(paste("Repeatability (R):", round(repeatability, 3)))
  
  ## HERITABILITY PER WINDOW
  library(lme4)
  library(dplyr)
  library(tidyr)
  
  # Initialize a dataframe to store per-window heritability estimates
  window_heritability <- data.frame(window = character(), heritability = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each methylation window
  for (window_name in names(data)[35:14818]) {
    temp_data <- data %>%
      select(lineage, all_of(window_name)) %>%
      rename(methylation = all_of(window_name)) %>%
      drop_na()
    
    if (n_distinct(temp_data$lineage) > 1) {  # Ensure at least 2 lineages are present
      model <- lmer(methylation ~ (1 | lineage), data = temp_data, REML = TRUE)
      
      var_components <- as.data.frame(VarCorr(model))
      var_lineage <- var_components$vcov[1]  # Lineage variance
      var_residual <- attr(VarCorr(model), "sc")^2  # Residual variance
      
      # Compute heritability (repeatability)
      heritability <- var_lineage / (var_lineage + var_residual)
      
      # Store results
      window_heritability <- rbind(window_heritability, data.frame(window = window_name, heritability = heritability))
    }
  }
  # CONSIDER IF THIS SHOULD BE MANUAL OR USE HERITABILITY FUNCTION / OR rptR FOR REPEATABILITY!!
  # also note, basically all of these are singular -- consider implications
  
  # Print summary of per-window heritability
  print(head(window_heritability))
  # 382 greater than 50%
  
  # Save results
  # write.csv(window_heritability, "window_heritability.csv", row.names = FALSE)
  
  ## HERITABILITY PER LINEAGE
  # Initialize dataframe for lineage heritability estimates
  lineage_heritability <- data.frame(lineage = character(), heritability = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each lineage
  for (lineage_name in unique(data$lineage)) {
    temp_data <- data %>%
      filter(lineage == lineage_name) %>%
      select(35:14818) %>%
      drop_na()
    
    if (nrow(temp_data) > 1) {  # Ensure at least 2 individuals exist
      lineage_mean_var <- apply(temp_data, 2, var, na.rm = TRUE)
      lineage_total_var <- var(as.vector(as.matrix(temp_data)), na.rm = TRUE)
      
      # Heritability estimate as mean variance within lineage / total variance
      heritability <- mean(lineage_mean_var, na.rm = TRUE) / lineage_total_var
      
      # Store results
      lineage_heritability <- rbind(lineage_heritability, data.frame(lineage = lineage_name, heritability = heritability))
    }
  }
  
  # Print summary of per-lineage heritability
  print((lineage_heritability))
  
  # Save results
  # write.csv(lineage_heritability, "lineage_heritability.csv", row.names = FALSE)
  
}

######################
## March 10, 2025
## approach for heritability H^2
{
  winds = data[,35:ncol(data)]
  meta = data[,1:34]
  
  #### Run correlation on gen 0 to gen 3 to see how heritable marks are across generations
  ### only grab Z score > 1.96 because those are significantly positively correlated. 
  ### more info here: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/what-is-a-z-score-what-is-a-p-value.htm
  {
    # Initialize vectors to store results
    cor_values <- c()
    
    # Loop through each methylation window
    for (col in colnames(data)[35:14818]) {
      
      # Subset data for gen 0 and gen 3
      gen0 <- data[data$generation == 0, c("lineage", col)]
      gen3 <- data[data$generation == 3, c("lineage", col)]
      
      # Compute mean methylation for gen 0 per lineage
      gen0_mean <- aggregate(gen0[[col]], by = list(gen0$lineage), mean)
      colnames(gen0_mean) <- c("lineage", "X")
      
      # Merge with gen 3 data to ensure matching lineages
      merged_data <- merge(gen0_mean, gen3, by = "lineage")
      colnames(merged_data)[3] <- "Y"
      
      # Compute Pearson correlation between X (gen 0 mean) and Y (gen 3)
      cor_result <- cor(merged_data$X, merged_data$Y, use = "complete.obs", method = "pearson")
      
      # Store correlation value
      cor_values <- c(cor_values, cor_result)
    }
    
    # Store results in a data frame
    regall <- data.frame(Window = colnames(data)[35:14818], Correlation = cor_values)
    
    # View first few results
    head(regall)
    
    # Compute mean and standard deviation of correlation values
    mean_cor <- mean(regall$Correlation, na.rm = TRUE)
    sd_cor <- sd(regall$Correlation, na.rm = TRUE)
    
    # Compute Z-scores
    regall$Zscore <- (regall$Correlation - mean_cor) / sd_cor
    
    # Identify significant windows (e.g., |Z| > 1.96 for ~5% significance level)
    # regall$Significant <- abs(regall$Zscore) > 1.96
    regall$Significant <- regall$Zscore > 1.96
    
    # Count significant windows
    sum(regall$Significant)  # Number of significant windows == 96 (when abs val of z) // == 48 when only > 1.96
    
    # View significant results
    significant_windows <- regall[regall$Significant == TRUE, ]
    head(significant_windows)
    
    hist(regall$Zscore, breaks = 50, main = "Distribution of Correlation Z-Scores",
         xlab = "Z-Score", col = "lightblue", border = "black")
    
    # Add significance cutoff lines (±1.96)
    abline(v = c(-1.96, 1.96), col = "red", lwd = 2, lty = 2)
  }
  
  sig_herit_windows <- significant_windows
  
  # Optionally, save the table as a CSV file
  # write.csv(sig_herit_windows, "significantly_03_heritable_windows.csv", row.names = FALSE)
  
  ### NOTE, MAY NEED TO ADJUST INDV IN LINEAGE 11 TO BE INCLUDED IN THIS!
  
  ## ok now compare gen3 and gen0, with sex separated
  ## NOTE, ONLY GRABS WHEN BOTH PARENTS ARE SEQUENCED
  {
    # Initialize result vectors
    cor_values <- c()
    sex_values <- c()
    window_names <- c()
    
    # Loop through each methylation window
    for (col in colnames(data)[35:14818]) {
      
      # Subset data for gen 0 and gen 3
      gen0 <- data[data$generation == 0, c("lineage", "sex", col)]
      gen3 <- data[data$generation == 3, c("lineage", col)]
      
      # Merge by lineage
      merged_data <- merge(gen0, gen3, by = "lineage")
      colnames(merged_data)[3:4] <- c("X", "Y")  # Rename columns for clarity
      
      # Compute correlation separately for males and females
      for (s in c(0, 1)) {  # 0 = Female, 1 = Male
        subset_data <- merged_data[merged_data$sex == s, ]
        
        if (nrow(subset_data) > 1) {  # Ensure enough data for correlation
          cor_result <- cor(subset_data$X, subset_data$Y, use = "complete.obs", method = "pearson")
          
          # Store results only if correlation is computed
          cor_values <- c(cor_values, cor_result)
          sex_values <- c(sex_values, s)
          window_names <- c(window_names, col)  # Keep track of window names
          #lineage_names <- c(lineage_names, subset_data$lineage)
        }
      }
    }
    
    # Store results in a data frame
    regall <- data.frame(Window = window_names, Correlation = cor_values, Sex = sex_values)
    
    # View first few results
    head(regall)
    
    # Compute mean and standard deviation for correlations
    mean_cor <- mean(regall$Correlation, na.rm = TRUE)
    sd_cor <- sd(regall$Correlation, na.rm = TRUE)
    
    # Compute Z-scores
    regall$Zscore <- (regall$Correlation - mean_cor) / sd_cor
    
    # Identify significant windows (e.g., |Z| > 1.96 for ~5% significance level)
    regall$Significant <- regall$Zscore > 1.96 #abs(regall$Zscore)
    
    # View significant windows
    significant_windows <- regall[regall$Significant == TRUE, ]
    head(significant_windows)
    nrow(significant_windows)
    
    table(significant_windows$Sex) # WOW!! most of the time, female influences!!
    
    dad = regall[regall$Sex==1,,drop=FALSE]
    mom = regall[regall$Sex==0,,drop=FALSE]
    
    t.test(x = mom$Correlation,
           y = dad$Correlation,
           alternative = "two.sided",
           paired = FALSE, 
           conf.level = 0.95)
    
    wilcox.test(Correlation ~ Sex, 
                data = regall, 
                exact = FALSE,
                conf.int = TRUE)
    
    hist(mom$Correlation, breaks = 100)
    hist(dad$Correlation, breaks = 100)
    
    regall$Sex = as.factor(regall$Sex)
    
    # Plot density with correlation on the y-axis
    ggplot(regall, aes(x = Correlation, fill = Sex)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(x = "Correlation", y = "Density", fill = "Sex") +
      theme(legend.position = "top")
    
    p1 <- ggplot(regall, aes(y = Correlation, fill = Sex)) +
      geom_density(data = filter(regall, Sex == 0), aes(x = ..density..), alpha = 1) +
      geom_density(data = filter(regall, Sex == 1), aes(x = -..density..), alpha = 1) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      ggtitle("3 Generation Jump") +
      theme_minimal() +
      labs(y = "Correlation", x = "Density (mirrored)", fill = "Sex") +
      theme_classic() # theme(legend.position = "top") +
    
    ggplot(regall, aes(y = Zscore, fill = Sex)) +
      geom_density(data = filter(regall, Sex == "0"), aes(x = ..density..), alpha = 0.5) +
      geom_density(data = filter(regall, Sex == "1"), aes(x = -..density..), alpha = 0.5) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      theme_minimal() +
      labs(y = "Zscore", x = "Density (mirrored)", fill = "Sex") +
      theme(legend.position = "top")
  }
  
  # same as above but for gen 1 to gen 3 AND gen 0 to gen 2
  {
    ## ok now compare gen3 and gen1 && gen2 and gen0, with sex separated
    ## NOTE, ONLY GRABS WHEN BOTH PARENTS ARE SEQUENCED
    # Initialize result vectors
    cor_values <- c()
    sex_values <- c()
    window_names <- c()
    
    # Loop through each methylation window
    for (col in colnames(data)[35:14818]) {
      
      # Subset data for gen 0 and gen 3
      gen1 <- data[data$generation == 1, c("lineage", "sex", col)]
      gen3 <- data[data$generation == 3, c("lineage", col)]
      
      gen0 <- data[data$generation == 0, c("lineage", "sex", col)]
      gen2 <- data[data$generation == 2, c("lineage", "indvs", col)]
      
      # Select the row with the smallest 'indvs' for each lineage for the OFFSPRING in that generation
      # Keep the smallest 'indvs' per lineage without changing format
      smallest_indvs <- gen2 %>%
        group_by(lineage) %>%
        filter(indvs == min(indvs)) %>%
        ungroup() %>%
        as.data.frame()  # Convert back to a standard data frame
      
      gen2 <- cbind(smallest_indvs$lineage,smallest_indvs[[3]])
      colnames(gen2) <- c("lineage", col)
      
      # Merge by lineage
      merged_data <- merge(gen1, gen3, by = "lineage")
      merged_data2 <- merge(gen0, gen2, by = "lineage")
      all_data <- rbind(merged_data, merged_data2)
      colnames(all_data)[3:4] <- c("X", "Y")  # Rename columns for clarity
      
      # Compute correlation separately for males and females
      for (s in c(0, 1)) {  # 0 = Female, 1 = Male
        subset_data <- all_data[all_data$sex == s, ]
        
        if (nrow(subset_data) > 1) {  # Ensure enough data for correlation
          cor_result <- cor(subset_data$X, subset_data$Y, use = "complete.obs", method = "pearson")
          
          # Store results only if correlation is computed
          cor_values <- c(cor_values, cor_result)
          sex_values <- c(sex_values, s)
          window_names <- c(window_names, col)  # Keep track of window names
        }
      }
    }
    
    # Store results in a data frame
    correlation13 <- data.frame(Window = window_names, Correlation = cor_values, Sex = sex_values)
    
    # View first few results
    head(correlation13)
    
    # Compute mean and standard deviation for correlations
    mean_cor <- mean(correlation13$Correlation, na.rm = TRUE)
    sd_cor <- sd(correlation13$Correlation, na.rm = TRUE)
    
    # Check if distributions are sig different
    wilcox.test(Correlation ~ Sex, 
                data = correlation13, 
                exact = FALSE, 
                conf.int = TRUE)
    
    # Compute Z-scores
    correlation13$Zscore <- (correlation13$Correlation - mean_cor) / sd_cor
    
    # Identify significant windows (e.g., |Z| > 1.96 for ~5% significance level)
    correlation13$Significant <- correlation13$Zscore > 1.96 #abs(correlation13$Zscore)
    
    # View significant windows
    significant_windows <- correlation13[correlation13$Significant == TRUE, ]
    head(significant_windows)
    nrow(significant_windows)
    
    table(significant_windows$Sex) # 43
    
    dad = correlation13[correlation13$Sex==1,,drop=FALSE]
    mom = correlation13[correlation13$Sex==0,,drop=FALSE]
    
    hist(mom$Correlation, breaks = 100)
    hist(dad$Correlation, breaks = 100)
    
    correlation13$Sex = as.factor(correlation13$Sex)
    
    # Plot density with correlation on the y-axis
    ggplot(correlation13, aes(x = Correlation, fill = Sex)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(x = "Correlation", y = "Density", fill = "Sex") +
      theme(legend.position = "top")
    
    p2 <- ggplot(correlation13, aes(y = Correlation, fill = Sex)) +
      geom_density(data = filter(correlation13, Sex == "0"), aes(x = ..density..), alpha = 1) +
      geom_density(data = filter(correlation13, Sex == "1"), aes(x = -..density..), alpha = 1) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      ggtitle("2 Generation Jump") +
      theme_minimal() +
      labs(y = "Correlation", x = "Density (mirrored)", fill = "Sex") +
      theme_classic()  # + theme(legend.position = "top") 
    
    ggplot(correlation13, aes(y = Zscore, fill = Sex)) +
      geom_density(data = filter(correlation13, Sex == "0"), aes(x = ..density..), alpha = 0.5) +
      geom_density(data = filter(correlation13, Sex == "1"), aes(x = -..density..), alpha = 0.5) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      theme_minimal() +
      labs(y = "Zscore", x = "Density (mirrored)", fill = "Sex") +
      theme(legend.position = "top")
  }
  
  # same as above but for gen 2 to gen 3 AND gen 0 to gen 1
  {
    ## ok now compare gen2 and gen1, with sex separated
    ## NOTE, ONLY GRABS WHEN BOTH PARENTS ARE SEQUENCED
    # Initialize result vectors
    cor_values <- c()
    sex_values <- c()
    window_names <- c()
    
    # Loop through each methylation window
    for (col in colnames(data)[35:14818]) {
      
      # Subset data for gen 0 and gen 3
      gen2 <- data[data$generation == 2, c("lineage", "sex", col)]
      gen3 <- data[data$generation == 3, c("lineage", col)]
      
      gen0 <- data[data$generation == 0, c("lineage", "sex", col)]
      gen1 <- data[data$generation == 1, c("lineage", "indvs",col)]
      
      # Select the row with the smallest 'indvs' for each lineage for the OFFSPRING in that generation
      # Keep the smallest 'indvs' per lineage without changing format
      smallest_indvs <- gen1 %>%
        group_by(lineage) %>%
        filter(indvs == min(indvs)) %>%
        ungroup() %>%
        as.data.frame()  # Convert back to a standard data frame
      
      gen1 <- cbind(smallest_indvs$lineage,smallest_indvs[[3]])
      colnames(gen1) <- c("lineage", col)
      
      # Merge by lineage
      merged_data <- merge(gen2, gen3, by = "lineage")
      merged_data2 <- merge(gen0, gen1, by = "lineage")
      all_data <- rbind(merged_data, merged_data2)
      colnames(all_data)[3:4] <- c("X", "Y")  # Rename columns for clarity
      
      # Compute correlation separately for males and females
      for (s in c(0, 1)) {  # 0 = Female, 1 = Male
        subset_data <- all_data[all_data$sex == s, ]
        
        if (nrow(subset_data) > 1) {  # Ensure enough data for correlation
          cor_result <- cor(subset_data$X, subset_data$Y, use = "complete.obs", method = "pearson")
          
          # Store results only if correlation is computed
          cor_values <- c(cor_values, cor_result)
          sex_values <- c(sex_values, s)
          window_names <- c(window_names, col)  # Keep track of window names
        }
      }
    }
    
    # Store results in a data frame
    correlation23 <- data.frame(Window = window_names, Correlation = cor_values, Sex = sex_values)
    
    # View first few results
    head(correlation23)
    
    # Compute mean and standard deviation for correlations
    mean_cor <- mean(correlation23$Correlation, na.rm = TRUE)
    sd_cor <- sd(correlation23$Correlation, na.rm = TRUE)
    
    # Check if distributions are sig different
    wilcox.test(Correlation ~ Sex, 
                data = correlation23, 
                exact = FALSE, 
                conf.int = TRUE)
    
    # Compute Z-scores
    correlation23$Zscore <- (correlation23$Correlation - mean_cor) / sd_cor
    
    # Identify significant windows (e.g., |Z| > 1.96 for ~5% significance level)
    correlation23$Significant <- correlation23$Zscore > 1.96 #abs(correlation23$Zscore)
    
    # View significant windows
    significant_windows <- correlation23[correlation23$Significant == TRUE, ]
    head(significant_windows)
    nrow(significant_windows)
    
    table(significant_windows$Sex) # 86
    
    dad = correlation23[correlation23$Sex==1,,drop=FALSE]
    mom = correlation23[correlation23$Sex==0,,drop=FALSE]
    
    hist(mom$Correlation, breaks = 100)
    hist(dad$Correlation, breaks = 100)
    
    correlation23$Sex = as.factor(correlation23$Sex)
    
    # Plot density with correlation on the y-axis
    ggplot(correlation23, aes(x = Correlation, fill = Sex)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(x = "Correlation", y = "Density", fill = "Sex") +
      theme(legend.position = "top")
    
    p3 <- ggplot(correlation23, aes(y = Correlation, fill = Sex)) +
      geom_density(data = filter(correlation23, Sex == "0"), aes(x = ..density..), alpha = 1) +
      geom_density(data = filter(correlation23, Sex == "1"), aes(x = -..density..), alpha = 1) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      theme_minimal() +
      ggtitle("1 Generation Jump") +
      labs(y = "Correlation", x = "Density (mirrored)", fill = "Sex") +
      theme(legend.position = "top") +
      theme_classic()
    
    ggplot(correlation23, aes(y = Zscore, fill = Sex)) +
      geom_density(data = filter(correlation23, Sex == "0"), aes(x = ..density..), alpha = 0.5) +
      geom_density(data = filter(correlation23, Sex == "1"), aes(x = -..density..), alpha = 0.5) +
      scale_fill_manual(values = c("0" = "pink", "1" = "lightblue")) +
      theme_minimal() +
      labs(y = "Zscore", x = "Density (mirrored)", fill = "Sex") +
      theme(legend.position = "top")
  }
  
  # first we compared male vs female lines
  # now, compare the generation jumps
  correlation23$jump <- "one"
  correlation13$jump <- "two"
  regall$jump <- "three"
  heritab <- rbind(regall, correlation23, correlation13)
  
  # Check if distributions are sig different
  kruskal.test(Correlation ~ jump, 
              data = heritab)
  kruskal.test(Correlation ~ Sex,
               data = heritab)
  pairwise.wilcox.test(heritab$Correlation, heritab$jump,
                       p.adjust.method = "BH", conf.int=TRUE)
  
  wilcox.test(Correlation ~ Sex, 
              data = heritab, 
              exact = FALSE,
              conf.int = TRUE)
  
  library(gridExtra)
  grid.arrange(p3,p2,p1, nrow = 1)
  
  Maternal <- c(35, 29, 23)
  Paternal <- c(51, 14, 1)
  Jumps <- c("1", "2", "3")
  bar <- data.frame(Maternal, Paternal, Jumps)
  bar <- data.frame(Paternal, Maternal, Jumps)
  
  library(tidyverse)
  bar %>% pivot_longer(cols = !Jumps, names_to = "Parent", values_to = "NumberHeritableWindows") %>% 
    ggplot(aes(fill = Parent, x = Jumps, y = NumberHeritableWindows)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = c("Paternal" = "lightblue", "Maternal" = "pink")) +
    labs(x = "Number of Generations", y = "Number of Heritable Windows", fill = "Parent") +
    theme_classic()
  
  bar %>%
    pivot_longer(cols = !Jumps, names_to = "Parent", values_to = "NumberHeritableWindows") %>%
    mutate(Parent = factor(Parent, levels = c("Paternal", "Maternal"))) %>%
    ggplot(aes(fill = Parent, x = Jumps, y = NumberHeritableWindows)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = c("Paternal" = "lightblue", "Maternal" = "pink")) +
    labs(x = "Number of Generations", y = "Number of Heritable Windows", fill = "Parent") +
    theme_classic()
  
  #### NEXT, Evaluate genetic correlations
  ## this uses the same code as above, however, you will need to change to:
  # for (col in colnames(data)[12:13])
  ## this works, because columns 12 and 13 are the inbreeding and kinship estimates
  
  #### next, look at half and full sib comparisons
  sibs <- data[data$indvs %in% c(3, 4, 8), ]
  
  # Extract methylation data (assuming methylation data is from columns 35 to 14818)
  methylation_matrix <- sibs[, 35:14818]
  
  # Ensure row names are set to individual IDs for easy matching
  rownames(methylation_matrix) <- sibs$indvs
  
  # Initialize results dataframe
  cor_results <- data.frame(Window = colnames(methylation_matrix), 
                            Diff_3_4 = NA, 
                            Diff_3_8 = NA, 
                            Diff_4_8 = NA)
  
  # Loop over windows (columns)
  for (i in 1:ncol(methylation_matrix)) {
    window_data <- methylation_matrix[, i]  # Extract column
    
    # calculate absolute differences, as pearsons correlations need more indv per group
    cor_results$Diff_3_4[i] <- abs(window_data[1] - window_data[2])
    cor_results$Diff_3_8[i] <- abs(window_data[1] - window_data[3])
    cor_results$Diff_4_8[i] <- abs(window_data[2] - window_data[3])
    
  }
  
  # since I'm looking for those that are consistent (i.e., heritable), I THINK I want SD/mean (CV) rather than zscores
  # cor_results$Z_Diff_3_4 <- scale(cor_results$Diff_3_4)
  # cor_results$Z_Diff_3_8 <- scale(cor_results$Diff_3_8)
  # cor_results$Z_Diff_4_8 <- scale(cor_results$Diff_4_8)
  # 
  # CV_Diff_3_4 <- sd(cor_results$Diff_3_4)/mean(cor_results$Diff_3_4)
  # CV_Diff_3_8 <- sd(cor_results$Diff_3_8)/mean(cor_results$Diff_3_8)
  # CV_Diff_4_8 <- sd(cor_results$Diff_4_8)/mean(cor_results$Diff_4_8)
  
  # head(cor_results)
  
  # instead, we want to merge sibs differences and get ONE CV value for all - to compare on same scale
  sibs_Diff_3_4_8 = NULL
  sibs_Diff_3_4_8$Diffs <- rbind(cor_results$Diff_3_4, cor_results$Diff_3_8, cor_results$Diff_4_8)
  CV_Diff_sibs <- sd(sibs_Diff_3_4_8$Diffs)/mean(sibs_Diff_3_4_8$Diffs)
  CV_Diff_sibs # 0.704 -- that means the difference in meth cutoff is < 1% !!
  
  sig34 <- cor_results[cor_results$Diff_3_4 < CV_Diff_sibs,,] # 581 ; full sib, opposite sex ## CV_Diff_3_4 -- 607
  sig38 <- cor_results[cor_results$Diff_3_8 < CV_Diff_sibs,,] # 466 ; 1/2 sib, same sex ## CV_Diff_3_8 -- 476
  sig48 <- cor_results[cor_results$Diff_4_8 < CV_Diff_sibs,,] # 196 ; 1/2 sib, opposite sex ## CV_Diff_4_8 -- 161
  
  # Create a Venn diagram
  # remember, this is for those with LOW variation in the DIFFERENCE among indv
  # consider if this is what we want exactly
  grid.arrange(newpage)
  venn.plot <- draw.triple.venn(
    area1 = nrow(sig34),         # 581
    area2 = nrow(sig38),         # 466
    area3 = nrow(sig48),         # 196
    n12 = nrow(intersect(sig34, sig38)),  # Overlap between 3_4 and 3_8
    n13 = nrow(intersect(sig34, sig48)),  # Overlap between 3_4 and 4_8
    n23 = nrow(intersect(sig38, sig48)),  # Overlap between 3_8 and 4_8
    n123 = nrow(intersect(intersect(sig34, sig38), sig48)),  # Overlap among all
    category = c("Full Sibs", "Half Sibs, Same Sex", "Half Sibs, Diff Sex"),  # c("3 vs 4", "3 vs 8", "4 vs 8")
    fill = c("red", "blue", "purple"),
    alpha = 0.5,
    cex = 2,
    cat.cex = 1.25,
    top = "Low Variation Among Sibs",
    main.just = c("center")
  )
  
  ### add non sibs
  sibs <- data[data$indvs %in% c(3, 4, 8, 9, 14, 32), ] # 9 OR 32 (or 13 or 14)
  # essentailly will compare sibs and non sibs
  
  # Extract methylation data (assuming methylation data is from columns 35 to 14818)
  methylation_matrix <- sibs[, 35:14818] #sibs[, 38:14821]
  
  # Ensure row names are set to individual IDs for easy matching
  rownames(methylation_matrix) <- sibs$indvs
  
  # Initialize results dataframe
  cor_results <- data.frame(Window = colnames(methylation_matrix), 
                            Diff_3_4 = NA, 
                            Diff_3_8 = NA, 
                            Diff_4_8 = NA,
                            Diff_9_14 = NA,
                            Diff_9_32 = NA,
                            Diff_14_32 = NA)
  
  # Loop over windows (columns)
  for (i in 1:ncol(methylation_matrix)) {
    window_data <- methylation_matrix[, i]  # Extract column
    
    # calculate absolute differences, as pearsons correlations need more indv per group
    cor_results$Diff_3_4[i] <- abs(window_data[1] - window_data[2])
    cor_results$Diff_3_8[i] <- abs(window_data[1] - window_data[3])
    cor_results$Diff_4_8[i] <- abs(window_data[2] - window_data[3])
    cor_results$Diff_9_14[i] <- abs(window_data[4] - window_data[5])
    cor_results$Diff_9_32[i] <- abs(window_data[4] - window_data[6])
    cor_results$Diff_14_32[i] <- abs(window_data[5] - window_data[6])
  }
  
  # instead, we want to merge sibs differences and get ONE CV value for all - to compare on same scale
  sibs_Diff_3_4_8_9_14_32 = NULL
  sibs_Diff_3_4_8_9_14_32$Diffs <- rbind(cor_results$Diff_3_4, cor_results$Diff_3_8, 
                                      cor_results$Diff_4_8, cor_results$Diff_9_14, 
                                      cor_results$Diff_9_32, cor_results$Diff_14_32)
  CV_Diff_sibs <- sd(sibs_Diff_3_4_8_9_14_32$Diffs)/mean(sibs_Diff_3_4_8_9_14_32$Diffs)
  CV_Diff_sibs # 0.774 -- that means the difference in meth cutoff is < 1% !!
  
  sig34 <- cor_results[cor_results$Diff_3_4 < CV_Diff_sibs,,] # 638 ; full sib, opposite sex 
  sig38 <- cor_results[cor_results$Diff_3_8 < CV_Diff_sibs,,] # 518 ; 1/2 sib, same sex 
  sig48 <- cor_results[cor_results$Diff_4_8 < CV_Diff_sibs,,] # 220 ; 1/2 sib, opposite sex 
  sig34[,2] <- "full_oppsx"
  sig38[,2] <- "half_smsx"
  sig48[,2] <- "half_oppsx"
  sig_sib_windows_june2025 <- rbind(sig34[,1:2], sig38[,1:2], sig48[,1:2]) # 1342 unique
  # write.table(sig_sib_windows_june2025, "significant_sibling_windows_june2025.csv", col.names = T, row.names = F, sep = ",")
  
  sig914 <- cor_results[cor_results$Diff_9_14 < CV_Diff_sibs,,] # 739 ; NONsib, opposite sex 
  sig932 <- cor_results[cor_results$Diff_9_32 < CV_Diff_sibs,,] # 900 ; NONsib, same sex 
  sig1432 <- cor_results[cor_results$Diff_14_32 < CV_Diff_sibs,,] # 662 ; NONsib, opposite sex 
  
  # by adding the non-sibs for the cuttoff, the cutoff got higher, but that increased the number of sig windows
  # will need to think about why there are more similar windows with the nonsibs
  
  # means for each row: 10.615, 12.616, 18.337, 10.043, 8.713, 10.813
  
  # Create a Venn diagram
  # remember, this is for those with LOW variation in the DIFFERENCE among indv
  # consider if this is what we want exactly
  grid.arrange(newpage)
  venn.plot <- draw.triple.venn(
    area1 = nrow(sig34),         # 
    area2 = nrow(sig38),         # 
    area3 = nrow(sig48),         # 
    n12 = nrow(intersect(sig34, sig38)),  # Overlap between 3_4 and 3_8
    n13 = nrow(intersect(sig34, sig48)),  # Overlap between 3_4 and 4_8
    n23 = nrow(intersect(sig38, sig48)),  # Overlap between 3_8 and 4_8
    n123 = nrow(intersect(intersect(sig34, sig38), sig48)),  # Overlap among all
    category = c("Full Sibs", "Half Sibs, Same Sex", "Half Sibs, Diff Sex"),  # c("3 vs 4", "3 vs 8", "4 vs 8")
    fill = c("red", "blue", "purple"),
    alpha = 0.5,
    cex = 2,
    cat.cex = 1.25,
    top = "Low Variation Among Sibs",
    main.just = c("center")
  )
  
  # run a T-test comparing the two groups
  # NOTE I am combining these differences, NOT averaging them
  t.test(x = cor_results[,2:4], 
         y = cor_results[,5:7], 
         alternative = "two.sided",
         paired = FALSE, 
         conf.level = 0.95)
  # significant p < 0.001

# NOTE consider adding indv 28 here -- have 1, 2, 3 gen offspring for that indv

######################
## March 31, 2025
## approach for comparing sig windows and locations across the genome

## updated june 6
# for heritable windows, grab genes of interest
a <- regall[regall$Significant == "TRUE",,drop=FALSE]
b <- correlation13[correlation13$Significant == "TRUE",, drop=FALSE]
c <- correlation23[correlation23$Significant == "TRUE",, drop=FALSE]
a$Significant <- "three"
b$Significant <- "two"
c$Significant <- "one"

abc <- rbind(a, b[,1:5], c)

setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/krat_epi/krat_epi/krat_epi_analysis/bwameth/Output") 
filt <- read.table("filt_data_82.csv", header=T, sep=",")
data = read.table("DATASET.csv", header=T, sep=",")
v_columns = grep("^V", names(data), value = TRUE) 

short <- matrix(v_columns, nrow = nrow(filt), ncol = 1)

winds <- cbind(short, filt)
window_names_82 = cbind(winds$short, winds$contig, winds$start, winds$end)
colnames(window_names_82) <- c("short", "contig", "start", "end")

abc_windows <- abc$Window
abc_windows <- window_names_82[window_names_82[,1] %in% abc_windows,,drop=FALSE]
write.csv(abc_windows, "signif_heritable_windows_contig_start_end_june2025.csv", row.names = FALSE)

sib_windows <- sig_sib_windows_june2025$Window
sig_sib_windows_june2025 <- window_names_82[window_names_82[,1] %in% sib_windows,,drop=FALSE]
write.csv(sig_sib_windows_june2025, "signif_sibling_windows_contig_start_end_june2025.csv", row.names = FALSE)

# sig_regall = regall[regall$Significant==TRUE,,drop=FALSE]
# inherited_windows = sig_regall$Window
# inherited_windows = window_names_82[window_names_82[,1] == sig_regall$Window,,]
# 
# # Subset window_names_82 where the first column matches any value in sig_regall$Window
# inherited_windows <- window_names_82[window_names_82[,1] %in% sig_regall$Window, ]

### DO THIS IN EASLEY
setwd("/home/gfl0003/kratepi/work/w_500_250")

# Load required packages
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)

# Load significant methylation windows
sig_windows <- read_csv("windows_of_note.csv")  # Ensure it has columns: "chr", "start", "end" # OLD
sig_windows <- read.csv("signif_heritable_windows_contig_start_end_june2025.csv") # NEW
sib_windows <- read.csv("signif_sibling_windows_contig_start_end_june2025.csv")

for(x in 1:nrow(sig_windows)){
  sig_windows$start[x] <- sig_windows$start[x] - 1000
  sig_windows$end[x] <- sig_windows$end[x] + 1000
}

# Load GFF file
gff_file <- "GCF_019054845.1_ASM1905484v1_genomic.gff"
gff_data <- import(gff_file, format = "gff")

# Convert significant windows to GRanges
sig_gr <- GRanges(seqnames = sig_windows$contig,
                  ranges = IRanges(start = sig_windows$start, end = sig_windows$end))

# Find overlaps
overlaps <- findOverlaps(sig_gr, gff_data)

# Extract overlapping annotation details
overlapping_annotations <- data.frame(
  chr = seqnames(gff_data)[subjectHits(overlaps)],
  start = start(gff_data)[subjectHits(overlaps)],
  end = end(gff_data)[subjectHits(overlaps)],
  feature = gff_data$type[subjectHits(overlaps)],  # Feature type (gene, exon, etc.)
  gene_id = mcols(gff_data)$ID[subjectHits(overlaps)]  # Adjust based on your GFF annotation
)

# Check if 'Dbxref' exists and extract Entrez ID safely
if ("Dbxref" %in% colnames(mcols(gff_data))) {
  entrez_ids <- mcols(gff_data)$Dbxref[subjectHits(overlaps)]
  
  # Ensure consistency by handling cases where Dbxref is a list
  entrez_ids <- sapply(entrez_ids, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)  # Handle missing values
    } 
    if (any(grepl("GeneID:", x))) {
      return(sub(".*GeneID:([0-9]+).*", "\\1", paste(x, collapse = ";")))  # Extract numeric Entrez ID
    } 
    return(NA)  # If no GeneID found, return NA
  })
  
  # Add the Entrez ID column
  overlapping_annotations$entrez <- entrez_ids
} else {
  overlapping_annotations$entrez <- NA  # If 'Dbxref' doesn't exist, fill with NA
}

# Merge with significant windows
annotated_windows <- cbind(sig_windows[queryHits(overlaps), ], overlapping_annotations)

# Save results
write_csv(annotated_windows, "annotated_significant_windows.csv") # old
write_csv(annotated_windows, "annotated_significant_heritable_windows_promoters_june2025.csv") # new
write_csv(annotated_windows, "annotated_significant_sibling_windows_promoters_june2025.csv")

#### NOW MOVE BACK TO LOCAL COMPUTER

anno = read.table("annotated_significant_heritable_windows_promoters_june2025.csv", header=T, sep=",")
anno = read.table("annotated_significant_sibling_windows_promoters_june2025.csv", header=T, sep=",")
table(anno$feature)

features = anno[anno$feature != "region",,drop=FALSE] # all 149 mapped to a region, which is good!
length(unique(features$short)) # 11 for herit, 58 for sibs

# write.table(features, "heritable_annotated_features_june2025.csv", col.names = T, row.names = F, sep=",")
# write.table(features, "sibling_annotated_features_june2025.csv", col.names = T, row.names = F, sep=",")

ann_wind_herit <- unique(features$short)
sig_ann_wind_herit <- abc[abc$Window %in% ann_wind_herit,,drop=FALSE]

# # Create pie chart
# pie(table(features$feature), 
#     main = "Distribution of Significant Windows Across the Genome", 
#     col = rainbow(length(unique(features$feature))))

# Convert table to dataframe
feature_counts <- as.data.frame(table(features$feature))
colnames(feature_counts) <- c("Feature", "Count")

# Create pie chart
ggplot(feature_counts, aes(x = "", y = Count, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Distribution of Significant Windows Across the Genome") +
  scale_fill_brewer(palette = "Set3")

which(features$short %in% inherited_windows) # V14601

### Cross check windows of interest with the metadata
sig_windows <- read_csv("windows_of_note.csv") 
windows_of_note <- as.data.frame(sig_windows)

meta = data[,1:37]
methyl = data[,38:ncol(data)]

# windows_of_note$short <- as.character(windows_of_note[,1]) # 878

# Ensure all values in grab$short exist in methyl column names
matching_columns <- colnames(methyl) %in% as.character(windows_of_note$short)

# Subset methyl by selecting only matching columns
methyl_subset <- methyl[, matching_columns, drop = FALSE]

# View results
head(methyl_subset)
dim(methyl_subset)
# this will now be used for the significance / AOV / Mr.PP / ridgeregression



####

# compare if there is a difference in methylation between oldest and newer samples
methvalues <- data[,35:14818]
avgmeth = NULL
methvalues$avgmeth <- rowMeans(methvalues)
data$avgmeth <- rowMeans(methvalues)

wind_avg <- colMeans(methvalues)
t.test(wind_avg)

summary(lm(avgmeth ~ year, data)

# Check if distributions are sig different
kruskal.test(Correlation ~ jump, 
             data = heritab)
kruskal.test(Correlation ~ Sex,
             data = heritab)
pairwise.wilcox.test(heritab$Correlation, heritab$jump,
                     p.adjust.method = "BH", conf.int=TRUE)

wilcox.test(Correlation ~ Sex, 
            data = heritab, 
            exact = FALSE,
            conf.int = TRUE)

kruskal.test(avgmeth ~ indvs, data)

t.test(data$avgmeth)

summary(lm(avgmeth ~ indvs, data)

library(rstatix)
kruskal_effsize(avgmeth ~ indvs, data = data, ci = TRUE, conf.level = 0.95)

library(effsize)
data$indvs <- as.factor(data$indvs)
cliff.delta(avgmeth ~ indvs, data = data, conf.level = 0.95, use.unbiased = TRUE, use.normal = TRUE)
