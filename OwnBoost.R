# A try for an own boosting algorithm

require(DirichletReg)
require(gamboostLSS)

source("trivariateDirichlet.R")

set.seed(100)

n = 1000

x1 = runif(n)
x2 = runif(n)
x3 = runif(n)
x = data.frame(x1,x2,x3)
x = as.matrix(x)

a1 = exp(-.7+2.5*x[,1] + 0.5*x[,2] + 1.5*x[,3])
a2 = exp(0.3 + 3*x[,1] + 0.5*x[,2] + 2*x[,3])
a3 = exp(-.5 + 2*x[,1] + 0.5*x[,2] + x[,3])

A = cbind(a1,a2,a3)

Y = rdirichlet(nrow(A),A)

par(mfrow=c(1,1),mar=c(4,4,1,1))
matplot(x,Y,pch=20)

# Initialize values (m = 0)
rate = 0.1

alp1 = exp(min(Y[,1]))
alp2 = exp(min(Y[,1]))
alp3 = exp(min(Y[,1]))

# lets try an alternative initialization (method of moments)

#alp1 = exp(((mean(Y[,1]) - mean(Y[,1]**2)) * mean(Y[,1])) / (mean(Y[,1]**2) - Y[,1]**2))
#alp2 = exp(((mean(Y[,1]) - mean(Y[,1]**2)) * mean(Y[,2])) / (mean(Y[,1]**2) - Y[,1]**2))
#alp3 = exp(((mean(Y[,1]) - mean(Y[,1]**2)) * mean(Y[,3])) / (mean(Y[,1]**2) - Y[,1]**2))


y1 = Y[,1]
y2 = Y[,2]
y3 = Y[,3]

coefficients = list()

intercept = list()
intercept[1] = min(Y[,1])
intercept[2] = min(Y[,2])
intercept[3] = min(Y[,3])
intercept

coefficients = list()
coefficients[[1]] = 0
coefficients[[2]] = 0
coefficients[[3]] = 0

loss = list()

for (i in 1:100000){
  
#Increase m by one and calculate gradient vector

u1 = (digamma(alp1 + alp2 + alp3) - digamma(alp1) + log(y1))
u2 = (digamma(alp1 + alp2 + alp3) - digamma(alp2) + log(y2))
u3 = (digamma(alp1 + alp2 + alp3) - digamma(alp3) + log(y3))

#fit seperate linear models for every covariate

u1d = lm(u1 ~ x)

alp1 = exp(log(alp1) + rate * fitted(u1d))

intercept[[1]] = intercept[[1]] + rate * coef(u1d)[1]
coefficients[[1]] = coefficients[[1]] + rate * coef(u1d)[2]



u2d = lm(u2 ~ x)

alp2 = exp(log(alp2) + rate * fitted(u2d))
intercept[[2]] = intercept[[2]] + rate * coef(u2d)[1]
coefficients[[2]] = coefficients[[2]] + rate * coef(u2d)[2]



u3d = lm(u3 ~ x)



alp3 = exp(log(alp3) + rate * fitted(u3d))

intercept[[3]] = intercept[[3]] + rate * coef(u3d)[1]
coefficients[[3]] = coefficients[[3]] + rate * coef(u3d)[2]



loss[[i]] = sum(-(lgamma(alp1 + alp2 + alp3) - (lgamma(alp1) + lgamma(alp2) + lgamma(alp3)) + ((alp1 - 1) * log(y1)
         + (alp2 - 1) * log(y2) + (alp3 - 1) * log(y3))))     

}

intercept
coefficients
View(loss)

# Fit Dirichlet for comparison
YD = DR_data(Y)
test = DirichReg(YD ~ x)
coef(test)

# Boosting via gamboostLSS
# optimal for non is mstop = 2883
model1 = glmboostLSS(Y ~ x, families = DirichletTV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL), control = boost_control(trace = TRUE, mstop = 5000, nu = 0.01), method = 'cyclic')
coef(model1)
plot(model1, parameter = "alpha2")


cv10f = cv(model.weights(model1), type = "kfold")
cvr = cvrisk(model1, fold = cv10f)
plot(cvr)

