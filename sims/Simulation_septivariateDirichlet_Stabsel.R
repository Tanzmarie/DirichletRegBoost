#########################################
library("batchtools")
library("gamboostLSS")
library("DirichletReg")
library("tidyverse")
library("forcats")
#########################################

source("families/septivariateDirichlet.R")

seed = 10 
set.seed(seed)

p = 50
n = 150

###### data generation
# - training data

x.train = matrix(runif(p * n, 0,1), n)
x.train = data.frame(x.train)

a1.train = exp(2.5*x.train[,1] - x.train[,2]) 
a2.train = exp(2*x.train[,3] + 2*x.train[,4])
a3.train = exp(1.5*x.train[,5] - 0.5*x.train[,6])
a4.train = exp(3*x.train[,7] - 2*x.train[,8]) 
a5.train = exp(0.5*x.train[,9] + 2.5*x.train[,10])
a6.train = exp(1.5*x.train[,11] - 3*x.train[,12])
a7.train = exp(2*x.train[,13] + 2*x.train[,14])

A = cbind(a1.train, a2.train, a3.train, a4.train, a5.train, a6.train, a7.train)

y.train = rdirichlet(nrow(A),A)

colnames(y.train) = c("y1","y2","y3", "y4", "y5", "y6", "y7")


############################################

# For simulation study: cutoff = (0.55, 0.99), q = [15,25,50], p = [50,100,150] #

septDR = glmboostLSS(y.train ~ ., data = x.train, families = DirichletSV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')

s = stabsel(septDR, cutoff = 0.9, PFER = 5)

trueVar = c("X1.alpha1" ,"X3.alpha2" ,"X5.alpha3" ,"X7.alpha4", "X9.alpha5", "X11.alpha6", "X13.alpha7",
            "X2.alpha1" ,"X4.alpha2" ,"X6.alpha3" ,"X8.alpha4", "X10.alpha5", "X12.alpha6", "X14.alpha7")

falseVar = names(s$max)[!names(s$max) %in% trueVar]
reInter = c("(Intercept).alpha1","(Intercept).alpha2","(Intercept).alpha3","(Intercept).alpha4","(Intercept).alpha5","(Intercept).alpha6","(Intercept).alpha7")
falseVar = falseVar[!falseVar %in% reInter]

CT = seq(0.55,0.99, 0.01) 
df = data.frame(TP = numeric(), FP = numeric(), PI = numeric(), p = numeric(), q = numeric(), PFER = numeric())

for (i in CT) {
  
s = stabsel(septDR, q = 15, cutoff = i)

tp.a1 = length(which(trueVar %in% names(selected(s)$alpha1)))
tp.a2 = length(which(trueVar %in% names(selected(s)$alpha2)))
tp.a3 = length(which(trueVar %in% names(selected(s)$alpha3)))
tp.a4 = length(which(trueVar %in% names(selected(s)$alpha4)))
tp.a5 = length(which(trueVar %in% names(selected(s)$alpha5)))
tp.a6 = length(which(trueVar %in% names(selected(s)$alpha6)))
tp.a7 = length(which(trueVar %in% names(selected(s)$alpha7)))

fp.a1 = length(which(falseVar %in% names(selected(s)$alpha1)))
fp.a2 = length(which(falseVar %in% names(selected(s)$alpha2)))
fp.a3 = length(which(falseVar %in% names(selected(s)$alpha3)))
fp.a4 = length(which(falseVar %in% names(selected(s)$alpha4)))
fp.a5 = length(which(falseVar %in% names(selected(s)$alpha5)))
fp.a6 = length(which(falseVar %in% names(selected(s)$alpha6)))
fp.a7 = length(which(falseVar %in% names(selected(s)$alpha7)))


TP = sum(tp.a1 + tp.a2 + tp.a3 + tp.a4 + tp.a5 + tp.a6 + tp.a7)
FP = sum(fp.a1 + fp.a2 + fp.a3 + fp.a4 + fp.a5 + fp.a6 + fp.a7)
PI = s$cutoff
p = 50
q = s$q
PFER = s$PFER

df = rbind(df, cbind(TP,FP,PI,p,q,PFER))

rm(s)

}


results$p150q50 = df

save(results, file = "BalancedStabsRes")
load("UnbalancedStabsRes")

############# PLOTTING ##################

res = do.call(rbind, results)


ggplot(res, aes(x = PI, y = TP)) +
  geom_point(color = "lightgrey") +
  geom_smooth(aes(y = FP), se = FALSE, size = 1, linetype = "dashed", color = "blue") +
  geom_smooth(se = FALSE, linewidth = 1, color = "black") +
  labs(x = "Threshold", y = "Number of True/False positives", title = "") +
  facet_grid(p ~ q, labeller = labeller(q = c("15" = "q = 15", "25" = "q = 25", "50" = "q = 50"),
                                        p = c("50" = "P = 50", "100" = "P = 100", "150" = "P = 150"))) +
  theme_light() +
  theme(strip.text = element_text(color = "black"))

ggplot(res, aes(x = PI, y = PFER)) +
  geom_line(color = "black") +
  #geom_smooth(aes(y = FP), se = FALSE, size = 1, linetype = "dashed", color = "blue") +
  #geom_smooth(se = FALSE, linewidth = 1, color = "black") +
  labs(x = "Threshold", y = "PFER", title = "") +
  facet_grid(p ~ q, labeller = labeller(q = c("15" = "q = 15", "25" = "q = 25", "50" = "q = 50"),
                                        p = c("50" = "P = 50", "100" = "P = 100", "150" = "P = 150"))) +
  theme_light() +
  theme(strip.text = element_text(color = "black"))



vec = c()
j = 1
for (i in CT) {
vec[j] = stabsel_parameters(q = 50, p = 1057, cutoff = i)$PFER
j = j + 1
}



results$p50q50$PFER = 0 

