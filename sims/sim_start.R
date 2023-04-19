# Simple simulation study for trivariate Dirichlet boosting.

require(DirichletReg)
require(gamboostLSS)

source("families/trivariateDirichlet.R")
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


A = cbind(a1,a2,a3)

y.train = rdirichlet(nrow(A),A)

# Fit Dirichlet for comparison
YD = DR_data(y.train)
test = DirichReg(YD ~ x.train)
coef(test)

# Boosting via gamboostLSS


#model1 = glmboostLSS(y.train ~ x.train, families = DirichletTV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL), control = boost_control(trace = TRUE, risk = 'oobag', mstop = 500, nu  = 0.1),  weights = weight.mstop, method = 'cyclic')
#coef(model1, off2int = TRUE)

#oobag.risk = risk(model1, merge = T)
#MSTOP = which.min(risk(model1, merge = T))

#rm(model1)

#x.train_biv <- x.train[weight.mstop == 1, ]

model1 = glmboostLSS(y.train ~ x.train, families = DirichletTV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL), control = boost_control(trace = TRUE, mstop = 1000, nu  = 0.1), method = 'cyclic')
coef(model1, off2int = TRUE)
summary(model1)
selected(model1)
print(model1)

#coefficients = coef(model1, off2int = TRUE)
#coeff_df = data.frame(Coefficient = names(coefficients$alpha1),
#                      Value = as.vector(coefficients$alpha1))

#ggplot(coeff_df, aes(x = Coefficient, y = Value)) +
#  geom_point(color = "steelblue") +
#  labs(x = "Coefficient", y = "Value") +
#  ggtitle("OLS Coefficients") +
#  theme_minimal() +
#  facet_grid(. ~ Coefficient, scales = "free_x", space = "free_x")




cv10f = cv(model.weights(model1), type = "kfold")
cvr = cvrisk(model1)
plot(cvr)
mstop(cvr)
model1[mstop(cvr)]


#(1) Variable Selection
# View True Positive and False Discovery Rate defined as:
#
# TPR = Number of correctly selected informative variables / Total number of informative variables 
# FDR = Number of falsely selected uninformative variables / Total number of selected variables
# Idea: Probing for Sparse and Fast Variable Selection with Model-Based Boosting, Janek and Hepp (2018)


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

#(3) Predictive Performance
# -> Negative Log-Likelihood (NLL)


