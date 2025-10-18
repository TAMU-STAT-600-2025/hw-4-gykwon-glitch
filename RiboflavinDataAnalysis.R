# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
fit <- fitLASSO(X, Y, n_lambda = 60)
# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
plot(fit$lambda_seq, colSums(fit$beta_mat != 0))
# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark)
microbenchmark(
  fitLASSO(X, Y, NULL, n_lambda = 60),
  times = 10
)
# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
#Unit: milliseconds
#expr      min       lq     mean   median       uq      max neval
#fitLASSO(X, Y, NULL, n_lambda = 60) 688.3569 697.4986 719.6364 707.0673 725.1173 791.5496    10
# -> 0.71 seconds
# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
cvfit <- cvLASSO(X, Y, n_lambda = 30)
# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
plot(cvfit$lambda_seq, cvfit$cvm, type="b", log="x",
     xlab=expression(lambda), ylab="CV(lambda)")
abline(v=cvfit$lambda_min, col="blue", lty=2)
abline(v=cvfit$lambda_1se, col="red",  lty=2)


