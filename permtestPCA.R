#' Permutation Test for PCA
#' 
#' @description Does the permutation test for principal components analysis. The output is a scree 
#' plot with alpha/2 and 1 - alpha/2 confidence intervals for each of the eigenvalues.
#' @export
#' @usage permtestPCA(X, nTests = 100, alpha = 0.05, center.data = TRUE, scale.data = TRUE, ...) 
#' @param X Numeric n by no of variables matrix on which PCA should be done.
#' @param nTests Number of permutations to be done. Default is 100. 
#' @param alpha Confidence level.
#' @param center.data Indicates whether the column means need to be removed. Default is TRUE.
#' @param scale.data Indicates whether the column needs to have sum of squares n. Default is TRUE.
#' @return A matrix with the lower and upper alpha/2 quantiles for each of the eigenvalues.
#' @examples
#' permtestPCA(gsp_share)
permtestPCA <- function (X, nTests = 100, alpha = 0.05, center.data = TRUE, scale.data = TRUE, ...){
  # Perform a PCA permutation test
  n <- nrow(X)
  m <- ncol(X)
  X <- scale(X, center = center.data, scale = scale.data)
  if (scale.data) {a <- 1/(n - 1)} else {a <- 1}
  
  res.X      <- prcomp(X)
  eigs.X     <- res.X$sdev^2 
  eigs.Xperm <- matrix(0, m, nTests)
  Xperm      <- matrix(0, n, m)
  Xperm[, 1] <- X[, 1];
  for (i in 1:nTests){
    for (j in 2:m) {
      ind <- sort(runif(n), index.return = TRUE)$ix   # Find random permutation of values 1:n
      Xperm[, j] <- X[ind, j]
    }
    res.Xperm  <- prcomp(Xperm)
    eigs.Xperm[, i] <- res.Xperm$sdev^2 
  }
  
  perc.alpha <- matrix(0, m, 2)
  for (s in 1:m){
    perc.alpha[s,] <- quantile(eigs.Xperm[s,], c(alpha/2, 1 - alpha/2) )
  }
  plot(1:m, eigs.X, type = "b", col = "red", main = "Permutation test PCA", xlab = "Component", ylab = "Eigenvalue", ...)
  lines(1:m, perc.alpha[, 1], type = "b", col="blue")
  lines(1:m, perc.alpha[, 2], type = "b", col="blue")
  
  #String1 <- sprintf('%4.1f%% Confidence',Alpha/2)
  #String2 <- sprintf('%4.1f%% Confidence',100-Alpha/2)
  #legend('Data',String1,String2)
  string1 <- paste("Confidence: ",formatC(alpha/2, digits=3, width=5, format="f"))
  string2 <- paste("Confidence: ",formatC(1-alpha/2, digits=3, width=5, format="f"))
  
  legend("topright", inset=.05, c("Observed", string1, string2), 
         lty = c(1, 1, 1), col = c("red", "blue", "blue"), pch = c("o", "o", "o"))
  return(perc.alpha)
}
