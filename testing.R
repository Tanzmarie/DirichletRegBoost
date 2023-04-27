# Simulation study for trivariate Dirichlet boosting in the low dimensional setting (p = 10, n.train = 500, n = 1000)

require(DirichletReg)
require(gamboostLSS)
require(tidyverse)
require(scoringRules)

source("families/trivariateDirichlet.R")

set.seed(100)

n.train = 1000
n.test = 1000
n.mstop = 1000

n = n.train + n.mstop

# weights needed for oobag risk estimation for noncyclical optimal mstop criterion
weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 

p = 10

a1 = c(2.5,-1,3, rep(0,p-3))
a2 = c(2,2,-1 ,rep(0,p-3))
a3 = c(1.5,-1.5,1 ,rep(0,p-3))

TrueBeta =  vector('list')
TrueBeta$alpha1 = a1
TrueBeta$alpha2 = a2
TrueBeta$alpha3 = a3

# Create train data

x.train = matrix(runif(p * n, 0,1), n)
x.train = data.frame(x.train)

a1.train = exp(sqrt(x.train[,1])) 
a2.train = exp(sin(x.train[,2]))
a3.train = exp(cos(x.train[,3]))


A = cbind(a1.train,a2.train,a3.train)

y.train = rdirichlet(nrow(A),A)

x.train_ncyc = x.train[weight.mstop == 1,]
y.train_ncyc = y.train[weight.mstop == 1,]

# Create test data

x.test = matrix(runif(p * n, 0,1), n)
x.test = data.frame(x.test)

a1.test = exp(sqrt(x.train[,1])) 
a2.test = exp(sin(x.train[,2]))
a3.test = exp(cos(x.train[,3]))


A = cbind(a1.test,a2.test,a3.test)

y.test = rdirichlet(nrow(A),A)

x.test_ncyc = x.test[weight.mstop == 1,]
y.test_ncyc = y.test[weight.mstop == 1,]
colnames(y.test_ncyc) = c("y1","y2","y3")

# Model based boosting

model2 = gamboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')

cv25 = cv(model.weights(model2), type = "kfold")
cvr = cvrisk(model2, folds = cv25, grid = 1:1000)
plot(cvr)
StopIT = mstop(cvr)

model2 = gamboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
plot(coef(model2), which = 1)

# DIRIG

ydr = DR_data(y.train_ncyc)
model4 = DirichReg(ydr ~ ., data = x.train_ncyc)
coef(model4)

