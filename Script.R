# Preliminaries.
rm(list = ls())

library(rio)
library(dplyr)
library(ggplot2)
library(boot)
library(plotrix)

data <- import('a2_data_group_9.csv')
getwd()

independent_matrix <- data %>% 
  select(!`Max. temperature (°C)`)

X <- independent_matrix[- nrow(independent_matrix), ]

Y <- data$`Max. temperature (°C)`[-1]

# Exploratory Analysis. 
par(mfrow = c(2, 2))

for(colname in colnames(X)){
  plot(X[, colname], Y,
  main = paste(colname, 'vs Max Temperature'),
  xlab = colname,
  ylab = 'Max Temperature (°C)',
  pch = 19, 
  col = 'lightgreen')
} # We have 19 different variables. 

colSums(is.na(X))
summary(X) # We need to standardize, PCA is very sensitive towards scales and outliers. 

# Let's do the data split. 
n <- nrow(X)
train_size <- floor(0.8 * n)

X_train <- X[1:train_size, ]
Y_train <- Y[1:train_size]

X_test  <- X[(train_size + 1):n, ]
Y_test  <- Y[(train_size + 1):n]

# Apply PCA excluding non numerical data such as Location and Date. 
X_train_PCA <- X_train %>% 
  select(!c(Date, `Location ID`))

apply(X_train_PCA, 2 , sd)
apply(X_train_PCA, 2 , mean)

fitpca <- princomp(X_train_PCA, cor = T, scores = T) # Argument SCORES makes the standardization. 

summary(fitpca)
# By Kaisers rule we go with 4 PC because the eigenvalues are bigger than 1. 
# By total variation 3 or 4 components already explain around 85% of the variance.
screeplot(fitpca, main = '') # The elbow is seen in the first component. 

# Let's take a look at the biplots. 

par(mfrow = c(1,2))

biplot(fitpca, pc.biplot = T, scale = 1, choices = 1:2, col = c('lightblue', 'red'),
       asp = 1, cex = c(0.5, 1), main = 'Biplot Componentes 1 and 2')

biplot(fitpca, pc.biplot = T, scale = 1, choices = 1:2, col = c('lightblue', 'red'),
       asp = 1, cex = c(0.5, 1), main = 'Biplot Componentes 2 and 3')

# This analysis is easy im just lazy to write it down, let's wait until we're writing the report.

print(summary(fitpca, loadings = T, cutoff = 0.3), ndigits = 2)

corpca <- cor(X_train_PCA, fitpca$scores)[,1:3] # 3 componentes. 
rowSums(corpca * corpca) # We can say that 96% of the variance of the variable Mean Temperature is explained by
# the first three principal components.
temp <- cbind(corpca, rowSums(corpca * corpca))
colnames(temp) <- c('Dim1', 'Dim2', 'Dim3','Fit')
temp <- round(temp, 2)
temp # All the variables are explained pretty well by the first 3 principal components except for Snowfall sum

# Let's find the confidence interval. 
par(mfrow = c(1,1))


plot(fitpca$sdev ** 2 , 
     type = 'b', 
     pch = 18, 
     cex = 1, 
     lwd = 1, 
     xlab = '', 
     ylab = 'Eigen Value')

nTests <- 1000
k = 0
eigs.Xperm = matrix(NA, 
                    nTests, 
                    ncol(X_train_PCA))

Xperm = as.data.frame(
  matrix(NA, 
         nrow(X_train_PCA), 
         ncol(X_train_PCA))
)

Xperm[, 1] <- X_train_PCA[ , 1];

for (i in 1:nTests){
  for (j in 2:ncol(X_train_PCA)){
    set.seed(k)
    ind <- sample(1:nrow(X_train_PCA), replace = F)
    Xperm[, j ] <- X_train_PCA[ind, j]
    k = k + 1
  }
  res.Xperm <- princomp(Xperm, cor = T, scores = T)
  eigs.Xperm[i, ] <- res.Xperm$sdev ** 2
}

plot(fitpca$sdev ** 2, type = 'b', pch = 18, cex = 1, lwd = 1, xlab = '', ylab = 'Variance Explained' )
lines(apply(eigs.Xperm, 2, quantile, c(0.025)), type = 'b', col = 'blue', pch = 18)
lines(apply(eigs.Xperm, 2, quantile, c(0.925)), type = 'b', col = 'blue', pch = 18)
# The first three components actually say something about the data. Im not entirely sure about component 4 we need to check. 

my_boot_pca <- function(x, ind){
  res <- princomp(x[ind, ], cor = TRUE)
  return(res$sdev^2)
}

fit.boot  <- boot(data = X_train_PCA, statistic = my_boot_pca, R = 1000)
eigs.boot <- fit.boot$t         
