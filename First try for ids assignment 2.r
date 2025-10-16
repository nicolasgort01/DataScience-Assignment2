library(readr)
library(tidyverse)
library(janitor)
library(ggplot2)
library(dplyr)
library(stringr)




# no scientific notations
options(scipen = 999)

# Read data
stopifnot(file.exists("a2_data_group_9.csv"))
data <- readr::read_csv("a2_data_group_9.csv", show_col_types = FALSE)

names(data)


# Preprocessing
data <- janitor::clean_names(data) # cleans all the cols names automatically 
names(data)


data <- data %>%
  select(-location_id)


# --- 1. Create X and y ---

# Build the X matrix: all variables, except the last row
X <- as.matrix(data[-nrow(data), ]) # this includes the date col


# Build the y vector: the max temperature column, except its first observation
y <- data$max_temperature_c[-1]

# Check alignment and dimensions
dim(X)  # should be (n-1) x p
length(y)  # should be (n-1)



# --- 2. Compute correlations ---

# 0) Ensure all-numeric predictors and aligned X/y
df_corr <- data %>% dplyr::select(where(is.numeric))

X_corr <- df_corr[-nrow(df_corr), ]            # data.frame (predictors at time t)
y_corr <- df_corr$max_temperature_c[-1]        # vector (target at time t+1)

# 1) Correlations
correlations <- sapply(X_corr, function(x) cor(x, y_corr, use = "complete.obs"))
correlations <- sort(correlations, decreasing = TRUE)
print(correlations)

# 2) Plot ALL variables
num_vars <- ncol(X_corr)
par(mfrow = c(4, 5))   # 4 rows × 5 columns grid (enough for 17 plots)

for (col in names(correlations)) {
  x <- X_corr[[col]]
  plot(x, y_corr,
       main = paste("Next-day Max Temp vs", col),
       xlab = col,
       ylab = "Next-day Max Temp (°C)",
       pch = 19, col = "darkblue")
  abline(lm(y_corr ~ x), col = "red", lwd = 2)
}

par(mfrow = c(1, 1))  # reset plotting layout


# 3) Correlation summary bar chart for the report

# correlations should already be a named numeric vector
corr_df <- tibble(
  variable = names(correlations),
  corr = as.numeric(correlations)
) %>%
  arrange(desc(corr)) %>%  # positives first, then negatives
  mutate(
    label = str_replace_all(variable, "_", " "),
    label = str_to_sentence(label),
    label = str_wrap(label, width = 25),
    label = factor(label, levels = rev(label))  # keep sorted order
  )

# Plot
ggplot(corr_df, aes(x = label, y = corr, fill = corr > 0)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  coord_flip() +
  labs(
    title = "Correlation with Next-day Max Temperature",
    x = NULL,
    y = "Pearson correlation"
  ) +
  scale_fill_manual(values = c("TRUE" = "#3B82F6", "FALSE" = "#EF4444")) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.2)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",      # remove legend
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )





# --- 3. Train Test Split (Time Series) ---

# 0) Ensure time order
data <- data %>% arrange(date)

# 1) Keep only numeric predictors (this automatically drops Date)
X_df <- data %>% select(where(is.numeric))  # includes max_temperature_c

# 2) Build aligned X (t) and y (t+1)
X <- as.matrix(X_df[-nrow(X_df), , drop = FALSE])
y <- data$max_temperature_c[-1]

# 3) Time indices (keep separately; do NOT include in X)
time_X <- data$date[-nrow(data)]
time_y <- data$date[-1]

# --- train/test split by proportion, no leakage ---
n <- nrow(X); k <- floor(0.8 * n)
X_train <- X[1:k, , drop = FALSE]; y_train <- y[1:k]
X_test  <- X[(k+1):n, , drop = FALSE]; y_test <- y[(k+1):n]


# --- 4. Generate PC on X_train ---

# 1. Perform PCA *only* on training data
pca_model <- prcomp(X_train, center = TRUE, scale. = TRUE)

# 2. Check how much variance each component explains
summary(pca_model)

# 3. Look at the contribution of each variable to each PC
round(pca_model$rotation,2)  




# --- 6. 2 methods to find the the optimal number of PC ---

# 1. Visualize variance explained 
plot(pca_model, type = "b", main = "Scree Plot") 



# 2. Eigenvalues (variances of PCs) and explained variance
eig  <- pca_model$sdev^2

# Show only eigenvalues > 1 (Kaiser's rule)
data.frame(
  PC = paste0("PC", seq_along(eig)),
  Eigenvalue = eig
) %>% 
  dplyr::filter(Eigenvalue > 1)


# 3. Using Cumulative VAF
prop_var <- eig / sum(eig)
cum_var  <- cumsum(prop_var)
tau   <- 0.85
k_vaf <- which(cum_var >= tau)[1]
k_vaf

cum_var

# 4
source("permtestPCA.R")          # Load the permtestPCA() function
perm_range <- permtestPCA(X_train)

# --- Permutation (Parallel Analysis) for choosing number of PCs ---

# X_train: numeric matrix/data.frame (observations x variables)
# B: number of permutations (e.g., 1000)
# alpha: significance level for the cutoff (0.05 uses 95th percentile of permuted eigs)
# Returns a list with suggested k and helpful objects, and draws a scree-style plot.
perm_parallel_pca <- function(X_train, B = 1000, alpha = 0.05, seed = 42,
                              center = TRUE, scale. = TRUE, make_plot = TRUE) {
  stopifnot(is.data.frame(X_train) || is.matrix(X_train))
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X_train)
  n <- nrow(X)
  p <- ncol(X)
  
  # 1) PCA on the real data
  pca_real <- prcomp(X, center = center, scale. = scale.)
  eig_real <- pca_real$sdev^2  # length p
  
  # 2) Permutation: shuffle rows within each column to break correlation, keep marginals
  perm_eigs <- replicate(B, {
    Xb <- apply(X, 2, function(col) sample(col, size = n, replace = FALSE))
    # PCA on permuted data
    pb <- prcomp(Xb, center = center, scale. = scale.)
    pb$sdev^2
  })
  
  perm_eigs <- t(perm_eigs)  # B x p (rows=replicates)
  
  # 3) For each component j, get the (1-alpha) quantile across permutations
  q_perm <- apply(perm_eigs, 2, quantile, probs = 1 - alpha, names = FALSE)
  
  # 4) Suggested k: how many components have real eigenvalue > permutation cutoff
  keep_vec <- eig_real > q_perm
  k_suggested <- if (any(keep_vec)) max(which(keep_vec)) else 0
  
  # 5) Helpful summaries
  vaf_real <- eig_real / sum(eig_real)
  cum_vaf_real <- cumsum(vaf_real)
  vaf_q_perm <- q_perm / sum(eig_real)        # put cutoffs on same scale for plotting
  cum_vaf_q_perm <- cumsum(vaf_q_perm)
  
  # 6) Plot
  if (make_plot) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    
    # Scree with cutoff
    plot(1:p, eig_real, type = "b", pch = 16,
         xlab = "Principal component", ylab = "Eigenvalue",
         main = sprintf("Parallel Analysis (Permutation): suggested k = %d", k_suggested))
    lines(1:p, q_perm, type = "b", pch = 1)
    abline(v = k_suggested + 0.5, lty = 3)  # visual separator after chosen k
    legend("topright",
           legend = c("Real data eigenvalues", sprintf("Permutation cutoff (%.0f%%)", (1 - alpha) * 100)),
           pch = c(16, 1), lty = 1, bty = "n")
  }
  
  list(
    k_suggested = k_suggested,
    eigen_real = eig_real,
    eigen_perm = perm_eigs,   # B x p matrix
    cutoff_quantile = q_perm, # length p (eigenvalue thresholds)
    keep = keep_vec,          # logical vector of length p
    vaf_real = vaf_real,
    cum_vaf_real = cum_vaf_real
  )
}


perm_range <- perm_parallel_pca(X_train)



# --- 7. 2 Construct the biplots ---

# Biplot of first 2 PCs
biplot(pca_model, scale = 0)


library(factoextra)
library(ggplot2)
library(patchwork)

# --- If your dataset is huge, optionally downsample *for plotting only* ---
set.seed(1)
plot_idx <- seq_len(min(4000, nrow(X_train)))  # or sample.int(nrow(X_train), 4000)
pca_plot <- prcomp(X_train[plot_idx, , drop = FALSE], center = TRUE, scale. = TRUE)

# A. PC1 vs PC2
p12 <- fviz_pca_biplot(
  pca_plot,
  axes = c(1, 2),
  label = "var",           # only label variables
  repel = TRUE,            # avoid overlap
  col.var = "#C43E3E",     # red loadings
  arrowsize = 0.7,
  col.ind = "skyblue3",    # points
  alpha.ind = 0.15,        # transparency
  pointsize = 0.6
) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Biplot Components 1 and 2", x = "Comp.1", y = "Comp.2")

# B. PC2 vs PC3
p23 <- fviz_pca_biplot(
  pca_plot,
  axes = c(2, 3),
  label = "var",
  repel = TRUE,
  col.var = "#C43E3E",
  arrowsize = 0.7,
  col.ind = "skyblue3",
  alpha.ind = 0.15,
  pointsize = 0.6
) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Biplot Components 2 and 3", x = "Comp.1", y = "Comp.2")

# Side-by-side layout
p12 + p23



# Loadings
loadings <- pca_model$rotation
round(loadings, 2)

# Optionally visualize
library(factoextra)
fviz_pca_biplot(pca_model)




# --- 8. Bootstrap Confidence Interval for choosen number of PCs ---


set.seed(42)

# Assumes you already have X_train (numeric matrix/data.frame of predictors only)
k <- 3

# Fit PCA once on the original training set (useful for reporting point estimate)
pca_train <- prcomp(X_train, center = TRUE, scale. = TRUE)
eig_train <- pca_train$sdev^2
vaf_k_hat <- sum(eig_train[1:k]) / sum(eig_train)   # point estimate

# Bootstrap (simple percentile CI)
B <- 1000
vaf_boot <- replicate(B, {
  idx <- sample.int(nrow(X_train), size = nrow(X_train), replace = TRUE)
  Xb  <- X_train[idx, , drop = FALSE]
  pca_b <- prcomp(Xb, center = TRUE, scale. = TRUE)
  eig_b <- pca_b$sdev^2
  sum(eig_b[1:k]) / sum(eig_b)
})

ci <- quantile(vaf_boot, c(0.025, 0.975))
list(
  point_estimate = vaf_k_hat,
  ci_95 = ci,
  mean_boot = mean(vaf_boot),
  sd_boot = sd(vaf_boot)
)


# Convert to data frame for ggplot
df_boot <- data.frame(vaf = vaf_boot)

# Plot density + CI + point estimate
ggplot(df_boot, aes(x = vaf)) +
  geom_histogram(aes(y = ..density..),
                 bins = 40, fill = "skyblue3", color = "white", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = vaf_k_hat, color = "red", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = ci[1], color = "darkred", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = ci[2], color = "darkred", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Bootstrap distribution of total variance explained (first 4 PCs)",
    x = "Total variance explained",
    y = "Density",
    subtitle = sprintf("Red line = point estimate (%.2f%%)\nDashed = 95%% CI [%.2f%%, %.2f%%]",
                       100*vaf_k_hat, 100*ci[1], 100*ci[2])
  ) +
  theme_minimal(base_size = 13)


#####
# Install/load boot if needed
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
library(boot)

set.seed(42)

# Inputs you already have:
# X_train: numeric matrix/data.frame of predictors
k <- 3
B <- 1000

# 1) Define a statistic function that returns ALL eigenvalues
#    (just like your professor's example)
eig_stat <- function(data, indices) {
  Xb <- data[indices, , drop = FALSE]
  pc <- prcomp(Xb, center = TRUE, scale. = TRUE)        # analogous to center=TRUE, scale.=TRUE
  pc$sdev^2                              # return eigenvalues
}

# 2) Run the bootstrap
fit.boot  <- boot(data = X_train, statistic = eig_stat, R = B)

# 3) Extract the R x p matrix of eigenvalues
eigs.boot <- fit.boot$t  # rows = bootstrap replicates, cols = eigenvalues

# 4) Convert those eigenvalues to VAF_k for each bootstrap replicate
vaf_boot <- apply(eigs.boot, 1, function(ev) sum(ev[1:k]) / sum(ev))

# 5) Point estimate via the same pipeline on the original sample
pc0 <- princomp(X_train, cor = TRUE)
eig0 <- pc0$sdev^2
vaf_k_hat <- sum(eig0[1:k]) / sum(eig0)

# 6) CI (percentile, to match your current approach)
ci <- quantile(vaf_boot, c(0.025, 0.975))

# 7) Summary output
summary_list <- list(
  point_estimate = vaf_k_hat,
  ci_95 = ci,
  mean_boot = mean(vaf_boot),
  sd_boot = sd(vaf_boot)
)
print(summary_list)

# 8) Plot the bootstrap distribution
hist(vaf_boot,
     breaks = 30,
     main = sprintf("Bootstrap Distribution of VAF for First %d PCs", k),
     xlab = "VAF (Variance Accounted For)",
     col = "skyblue", border = "white")
abline(v = vaf_k_hat, col = "red", lwd = 2, lty = 2)
abline(v = ci, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", legend = c("Point Estimate", "95% CI"),
       col = c("red", "darkgreen"), lty = c(2, 3), bty = "n")

####









# --- 9. Fit the choosen PCs (k) -----


if (!requireNamespace("pls", quietly = TRUE)) install.packages("pls")
library(pls)

set.seed(42)

k <- k

# Fit PCR with chosen number of components
pcr_fit <- pcr(
  y_train ~ .,
  data = data.frame(y_train = y_train, X_train),
  scale = TRUE,
  center = TRUE,
  validation = "CV",
  ncomp = k   # <-- specify number of components to use
)

# Check model summary
summary(pcr_fit)


# --- 10. Fit a benchmark MLR -----

# Fit a standard multiple linear regression (no dimensionality reduction)
lm_fit <- lm(y_train ~ ., data = data.frame(y_train = y_train, X_train))

# Check model summary
summary(lm_fit)


# PCR predictions (specify number of components = 3)
pcr_pred_test <- predict(pcr_fit, newdata = data.frame(X_test), ncomp = 3)

# MLR predictions
lm_pred_test <- predict(lm_fit, newdata = data.frame(X_test))




# --- 11. Performance metrics -----

# Define RMSE function
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# Calculate RMSE for test set
lm_rmse_test  <- rmse(y_test, lm_pred_test)
pcr_rmse_test <- rmse(y_test, pcr_pred_test)

# Display results
lm_rmse_test
pcr_rmse_test

# Compute R² for PCR train set
pcr_pred_train <- predict(pcr_fit, newdata = data.frame(X_train), ncomp = 3)
ss_res_train <- sum((y_train - pcr_pred_train)^2)
ss_tot_train <- sum((y_train - mean(y_train))^2)
pcr_r2_train <- 1 - (ss_res_train / ss_tot_train)

pcr_r2_train


# Compute R² for PCR test set
ss_res <- sum((y_test - pcr_pred_test)^2)
ss_tot <- sum((y_test - mean(y_test))^2)
pcr_r2_test <- 1 - (ss_res / ss_tot)

pcr_r2_test


# --- 12. k + 1 and k - 1 -----

# Fit PCR with chosen number of components + 1 (k + 1)
pcr_fit1 <- pcr(
  y_train ~ .,
  data = data.frame(y_train = y_train, X_train),
  scale = TRUE,
  center = TRUE,
  validation = "CV",
  ncomp = k + 1   # <-- specify number of components to use
)


summary(pcr_fit1)

# PCR predictions
pcr_pred_test <- predict(pcr_fit1, newdata = data.frame(X_test), ncomp = k + 1)


# Fit PCR with chosen number of components - 1 (k - 1)
pcr_fit2 <- pcr(
  y_train ~ .,
  data = data.frame(y_train = y_train, X_train),
  scale = TRUE,
  center = TRUE,
  validation = "CV",
  ncomp = k - 1   # <-- specify number of components to use
)

summary(pcr_fit2)


# PCR predictions
pcr_pred_test <- predict(pcr_fit2, newdata = data.frame(X_test), ncomp = k - 1)





# --- 13 -----



# Combine observed and predicted values
results <- data.frame(
  y_test = y_test,
  lm_pred = lm_pred_test,
  pcr_pred = as.numeric(pcr_pred_test)
)

# Create decile groups based on observed y_test values
results$group <- cut(
  results$y_test,
  breaks = quantile(results$y_test, probs = seq(0, 1, 0.1)),
  include.lowest = TRUE,
  labels = paste0("G", 1:10)
)

rmse_by_group <- results %>%
  group_by(group) %>%
  summarise(
    RMSE_LM  = rmse(y_test, lm_pred),
    RMSE_PCR = rmse(y_test, pcr_pred)
  )

# Convert to long format for easier plotting
rmse_long <- rmse_by_group %>%
  tidyr::pivot_longer(cols = c(RMSE_LM, RMSE_PCR),
                      names_to = "Model",
                      values_to = "RMSE")

# Plot
ggplot(rmse_long, aes(x = group, y = RMSE, color = Model, group = Model)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "RMSE by temperature decile (test data)",
    x = "Decile group (based on observed max temperature)",
    y = "RMSE",
    color = "Model"
  )

