# Simple simulation study for trivariate Dirichlet boosting.

require(DirichletReg)
require(gamboostLSS)

source("families/trivariateDirichlet.R")

n = 1000

x1 = runif(n)
x2 = runif(n)
x = data.frame(x1,x2)
x = as.matrix(x)

a1 = exp(-.7+2.5*x[,1] + 0.5*x[,2])
a2 = exp(0.3 + 3*x[,1] + 2*x[,2])
a3 = exp(-.5 + 2*x[,1] + 1.5*x[,2])

A = cbind(a1,a2,a3)

Y = rdirichlet(nrow(A),A)



par(mfrow=c(1,1),mar=c(4,4,1,1))
matplot(x,Y,pch=20)

# Perform a 80-20 train-test split
set.seed(123) # Set a seed for reproducibility

n <- nrow(x)

train_indices <- sample(1:n, floor(0.7*n)) # 80% for training
x_train <- x[train_indices, ]
x_test <- x[-train_indices, ]
y_train = Y[train_indices, ]
y_test = Y[-train_indices, ]

# Fit Dirichlet for comparison
YD = DR_data(Y)
test = DirichReg(YD ~ x)
coef(test)

# Boosting via gamboostLSS
# optimal for non is mstop = 2883
model1 = glmboostLSS(Y ~ x, families = DirichletTV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL), control = boost_control(trace = TRUE, mstop = 5000, nu = 0.1), method = 'cyclic')

#cv10f = cv(model.weights(model1), type = "kfold")
#cvr = cvrisk(model1, fold = cv10f)
#plot(cvr)

MSTOP <- which.min(risk(model1,merge = T))

