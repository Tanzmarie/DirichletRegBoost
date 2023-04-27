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

a1.train = exp(2.5*x.train[,1] - x.train[,2] + 3*x.train[,3]) 
a2.train = exp(2*x.train[,1] + 2*x.train[,2] - x.train[,3])
a3.train = exp(1.5*x.train[,1] -  1.5*x.train[,2] + x.train[,3])


A = cbind(a1.train,a2.train,a3.train)

y.train = rdirichlet(nrow(A),A)

x.train_ncyc = x.train[weight.mstop == 1,]
y.train_ncyc = y.train[weight.mstop == 1,]

# Create test data

x.test = matrix(runif(p * n, 0,1), n)
x.test = data.frame(x.test)

a1.test = exp(2.5*x.test[,1] - x.test[,2] + 3*x.test[,3]) 
a2.test = exp(2*x.test[,1] + 2*x.test[,2] - x.test[,3])
a3.test = exp(1.5*x.test[,1] -  1.5*x.test[,2] + x.test[,3])


A = cbind(a1.test,a2.test,a3.test)

y.test = rdirichlet(nrow(A),A)

x.test_ncyc = x.test[weight.mstop == 1,]
y.test_ncyc = y.test[weight.mstop == 1,]
colnames(y.test_ncyc) = c("y1","y2","y3")

# Boosting via glmboostLSS and finding mstop for cyclical algorithm


model1 = glmboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = c(alpha1 = 500, alpha2 = 500, alpha3 = 500), nu = 0.1), method = 'cyclic')
coef(model1, off2int = TRUE)

grid1 = make.grid(max = c(alpha1 = 500, alpha2 = 500, alpha3 = 500), min = 20, length.out = 10, dense_mu_grid = FALSE)
cv10f = cv(model.weights(model1), type = "kfold")
cvr = cvrisk(model1, folds = cv10f, grid = grid1)
plot(cvr)
mstop(cvr)

model1[mstop(cvr)]
coef(model1[mstop(cvr)],off2int = TRUE)

par(mfrow = c(1,3))
plot(model1[mstop(cvr)])

# Boosting via glmboost and finding mstop via cvrisk for noncyclical algorithm

model2 = glmboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')

cv25 = cv(model.weights(model2), type = "kfold")
cvr = cvrisk(model2, folds = cv25, grid = 1:1000)
plot(cvr)
StopIT = mstop(cvr)

model2 = glmboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
coef(model2, off2int = TRUE)

# Maybe over oobag?

model3 = glmboostLSS(y.train ~ . ,data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 2000, risk = "oobag", nu = 0.1), method = "noncyclic", weights = weight.mstop)

oobag.risk = risk(model3, merge = TRUE)
MSTOP = which.min(risk(model3, merge = TRUE))

model3 = glmboostLSS(y.train_ncyc ~ ., data = x.train_ncyc, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = MSTOP, nu = 0.1), method = "noncyclic")
coef(model3, off2int = TRUE)
par(mfrow = c(1,3))
plot(model3)

# Dirichlet Regression

ydr = DR_data(y.train_ncyc)
model4 = DirichReg(ydr ~ ., data = x.train_ncyc)

#(1) Variable Selection
# View True Positive and False Discovery Rate defined as:
#
# TPR = Number of correctly selected informative variables / Total number of informative variables 
# FDR = Number of falsely selected uninformative variables / Total number of selected variables
# Idea: Probing for Sparse and Fast Variable Selection with Model-Based Boosting, Janek and Hepp (2018)

selectedVar.a1.risk = vector("list")
selectedVar.a2.risk = vector("list")
selectedVar.a3.risk = vector("list")

true.positive.a1.risk = vector("list")
true.positive.a2.risk = vector("list")
true.positive.a3.risk = vector("list")

false.positive.a1.risk = vector('list')
false.positive.a2.risk = vector('list')
false.positive.a3.risk = vector('list')

nameVar = names(x.train)[1:p]
trueVar = nameVar[1:3]
falseVar = nameVar[4:p]

selectedVar.a1.risk = names(coef(model2$alpha1))[-1]
selectedVar.a2.risk = names(coef(model2$alpha2))[-1]
selectedVar.a3.risk = names(coef(model2$alpha3))[-1]

true.positive.a1.risk = length(which(trueVar %in% names(coef(model2$alpha1))))
true.positive.a2.risk = length(which(trueVar %in% names(coef(model2$alpha2))))
true.positive.a3.risk = length(which(trueVar %in% names(coef(model2$alpha3))))

false.positive.a1.risk = length(which(falseVar %in% names(coef(model2$alpha1))))
false.positive.a2.risk = length(which(falseVar %in% names(coef(model2$alpha2))))
false.positive.a3.risk = length(which(falseVar %in% names(coef(model2$alpha3))))


TPR = vector("list")
FDR = vector ("list")

TPR["Overall"] = sum(true.positive.a1.risk + true.positive.a2.risk + true.positive.a3.risk) / (length(trueVar) * 3)
TPR["alpha1"] = true.positive.a1.risk / length(trueVar)
TPR["alpha2"] = true.positive.a2.risk / length(trueVar)
TPR["alpha3"] = true.positive.a3.risk / length(trueVar)

FDR["Overall"] = sum(false.positive.a1.risk + false.positive.a2.risk + false.positive.a3.risk) / (length(selectedVar.a1.risk) + length(selectedVar.a2.risk) + length(selectedVar.a3.risk))
FDR["alpha1"] = false.positive.a1.risk / length(selectedVar.a1.risk)
FDR["alpha2"] = false.positive.a2.risk / length(selectedVar.a2.risk)
FDR["alpha3"] = false.positive.a3.risk / length(selectedVar.a3.risk)

TPR = t(data.frame(TPR))
colnames(TPR) = "TPR"
FDR = t(data.frame(FDR))
colnames(FDR) = "FDR"

cbind(TPR,FDR)

#(2) Estimation
# Plotting coefficients via facet_grid()
#
#coefficients = coef(model1, off2int = TRUE)
#coeff_df = data.frame(Coefficient = names(coefficients$alpha1),
#                      Value = as.vector(coefficients$alpha1))

#ggplot(coeff_df, aes(x = Coefficient, y = Value)) +
#  geom_point(color = "steelblue") +
#  labs(x = "Coefficient", y = "Value") +
#  ggtitle("OLS Coefficients") +
#  theme_minimal() +
#  facet_grid(. ~ Coefficient, scales = "free_x", space = "free_x")
#
# -> Needs further evaluation

# multiple solutions in one plot via facet_grid

test = data.frame(coef(model4))
test = test[2:11,]
test$Var = row.names(test)
test$Var = factor(test$Var, level = nameVar)
colnames(test)[1:3] = c("alpha1", "alpha2", "alpha3")

cof_df = data.frame(coef(model2, which = c(2:11)))
cof_df$Var = row.names(cof_df)
cof_df$Var = factor(cof_df$Var, level = nameVar)

TBet = data.frame(TrueBeta)
TBet$Var = row.names(cof_df)

df = gather(cof_df, alpha, Val, -Var)
true.df = gather(TBet, alpha, Val, -Var)
DR = gather(test, alpha, Val, -Var)

ggplot() +
  geom_point(data = df, aes(x = Var, y = Val), color = "black") +
  geom_point(data = true.df, aes(x = Var, y = Val), color = "red", size = 2.5) +
 # geom_point(data = DR, aes(x = Var, y = Val), color = "green") +
  facet_grid(rows = vars(alpha)) +
  theme_light()



#(3) Predictive Performance
# -> Negative Log-Likelihood (NLL)
# -> MSEP
# -> Energy score

# Start by doing predictions

pred.a1 = predict(model2$alpha1, newdata = x.test_ncyc, type = "response")
pred.a2 = predict(model2$alpha2, newdata = x.test_ncyc, type = "response")
pred.a3 = predict(model2$alpha3, newdata = x.test_ncyc, type = "response")
pred.A = cbind(pred.a1,pred.a2,pred.a3)
pred.mu = pred.A / rowSums(pred.A)

pred.DR = predict(model4, newdata = x.test_ncyc, mu = TRUE)

# MSEP 

MSEPB = vector("list")
MSEPDR = vector("list")

MSEPB$alpha1 = mean((pred.mu[,1] - y.test_ncyc[,1])**2)
MSEPB$alpha2 = mean((pred.mu[,2] - y.test_ncyc[,2])**2)
MSEPB$alpha3 = mean((pred.mu[,3] - y.test_ncyc[,3])**2)

MSEPDR$alpha1 = mean((pred.DR[,1] - y.test_ncyc[,1])**2)
MSEPDR$alpha2 = mean((pred.DR[,2] - y.test_ncyc[,2])**2)
MSEPDR$alpha3 = mean((pred.DR[,3] - y.test_ncyc[,3])**2)

MSEPB = t(data.frame(MSEPB))
MSEPDR = t(data.frame(MSEPDR))

MSEP = vector("list")
MSEP$MSEP = cbind(MSEPB,MSEPDR)
colnames(MSEP$MSEP) = c("Boosting", "DirichReg")
MSEP

# NLL loss

pred.DR = predict(model4, newdata = x.test_ncyc, mu = FALSE, alpha = TRUE)

NLL = vector("list")

loss = function(alpha1, alpha2, alpha3, y) {
  y3 = y[,3]
  y2 = y[,2]
  y1 = y[,1]
  
  - (lgamma(alpha1 + alpha2 + alpha3) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3))
     + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3)))
  
}

NLL$Boosting = sum(loss(alpha1 = pred.a1, alpha2 = pred.a2, alpha3 = pred.a3, y = y.test_ncyc))
NLL$DirichReg = sum(loss(alpha1 = pred.DR[,1], alpha2 = pred.DR[,2], alpha3 = pred.DR[,3], y = y.test_ncyc))
NLL

# Energy Score

es_boost = vector()
es_DR = vector()

for (i in 1:length(pred.a1)) {
  
  pred_sample_boost = matrix(NA, nrow = 3, ncol = 10000)
  pred_sample_DR = matrix(NA, nrow = 3, ncol = 10000)
  
  # Boosting approach
  
  sample_boost = rdirichlet(10000, pred.A[i,])
  
  pred_sample_boost[1,] = sample_boost[,1]
  pred_sample_boost[2,] = sample_boost[,2]
  pred_sample_boost[3,] = sample_boost[,3]

  es_boost[i] <- es_sample(y = c(y.test_ncyc[i,1], y.test_ncyc[i,2], y.test_ncyc[i,3]), dat = pred_sample_boost) 
  
  # DirichReg 
  
  sample_DR = rdirichlet(10000, pred.DR[i,])
  
  pred_sample_DR[1,] = sample_DR[,1]
  pred_sample_DR[2,] = sample_DR[,2]
  pred_sample_DR[3,] = sample_DR[,3]
  
  es_DR[i] <- es_sample(y = c(y.test_ncyc[i,1], y.test_ncyc[i,2], y.test_ncyc[i,3]), dat = pred_sample_DR)
  
  
}

energy_score = list()
energy_score$Boosting = mean(es_boost)
energy_score$DirichReg = mean(es_DR)
energy_score


