# Simple simulation study for trivariate Dirichlet boosting.

require(DirichletReg)
require(gamboostLSS)

source("families/septivariateDirichlet.R")
set.seed(100)

n.train = 500
n.mstop = 500

n = n.train + n.mstop
#weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
p = 2

x.train = matrix(runif(p * n, 0,1), n)

a1 = exp(-.7+2.5*x.train[,1] + 0.5*x.train[,2])
a2 = exp(0.3 + 3*x.train[,1] + 2*x.train[,2])
a3 = exp(-.5 + 2*x.train[,1] + 1.5*x.train[,2])
a4 = exp(-.7+5*x.train[,1] + 0.5*x.train[,2])
a5 = exp(0.3 + x.train[,1] + 4*x.train[,2])
a6 = exp(-.5 + 2.4*x.train[,1] + x.train[,2])
a7 = exp(-.5 + 1.5*x.train[,1] + 1.5*x.train[,2])


A = cbind(a1,a2,a3,a4,a5,a6,a7)

y.train = rdirichlet(nrow(A),A)

# Fit Dirichlet for comparison
YD = DR_data(y.train)
test = DirichReg(YD ~ x.train)
coef(test)

# Boosting test

model1 = glmboostLSS(y.train ~ x.train, families = DirichletSV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL), control = boost_control(trace = TRUE, mstop = 1000, nu  = 0.1), method = 'noncyclic')
coef(model1, off2int = TRUE)
summary(model1)
selected(model1)
print(model1)