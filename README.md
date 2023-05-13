# Boosting Dirichlet regression models 
Component-wise gradient boosting algorithm for modeling Dirichlet regression models within the framework of generalized additive models for location, scale and shape, which enables the simultaneous modeling of all distribution parameters of the Dirichlet distribution of a multivariate response conditional on
explanatory variables.

# Example 

require("gamboostLSS")
require("DirichletReg")

source("families/trivariateDirichlet.R")

set.seed(1)

n = 100
p = 10

x = matrix(runif(p * n, 0,1), n)

x = data.frame(x)

a1 = exp(2.5*x[,1] - x[,2] + 3*x[,3]) 
a2 = exp(2*x[,4] + 2*x[,5] - x[,6])
a3 = exp(1.5*x[,7] -  1.5*x[,8] + x[,9])
A = cbind(a1,a2,a3)

y = rdirichlet(nrow(A),A)

colnames(y) = c("y1","y2","y3")


# --- model

mod = glmboostLSS(y ~ ., data = x, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')

coef(mod[200], off2int = TRUE)

par(mfrow = c(1,3))
plot(mod)
