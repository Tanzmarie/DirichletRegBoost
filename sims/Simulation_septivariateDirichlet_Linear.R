#########################################
library("batchtools")
library("gamboostLSS")
library("scoringRules")
library("DirichletReg")
library("tidyverse")
library("reshape2")
#########################################

source("families/septivariateDirichlet.R")


sim = function(seed,n,p) {
  
  set.seed(seed)
  
  a1 = c(2.5,-1,5, rep(0,p-3))
  a2 = c(2,2,-4 ,rep(0,p-3))
  a3 = c(1.5,-1.5,1 ,rep(0,p-3))
  a4 = c(3,-2,-1, rep(0,p-3))
  a5 = c(0.5,2.5,-1.5 ,rep(0,p-3))
  a6 = c(1.5,-3,-1 ,rep(0,p-3))
  a7 = c(2,2,2 ,rep(0,p-3))
  
  
  TrueBeta =  vector('list')
  
  TrueBeta$alpha1 = a1
  TrueBeta$alpha2 = a2
  TrueBeta$alpha3 = a3
  TrueBeta$alpha4 = a4
  TrueBeta$alpha5 = a5
  TrueBeta$alpha6 = a6
  TrueBeta$alpha7 = a7
  
  ###### data generation
  # - training data
  
  x.train = matrix(runif(p * n, 0,1), n)
  x.train = data.frame(x.train)
  
  a1.train = exp(2.5*x.train[,1] - x.train[,2] + 5*x.train[,3]) 
  a2.train = exp(2*x.train[,1] + 2*x.train[,2] - 4*x.train[,3])
  a3.train = exp(1.5*x.train[,1] -  1.5*x.train[,2] + x.train[,3])
  a4.train = exp(3*x.train[,1] - 2*x.train[,2] + -1*x.train[,3]) 
  a5.train = exp(0.5*x.train[,1] + 2.5*x.train[,2] - 1.5*x.train[,3])
  a6.train = exp(1.5*x.train[,1] - 3*x.train[,2] + -1*x.train[,3])
  a7.train = exp(2*x.train[,1] + 2*x.train[,2] + 2*x.train[,3])
  
  A = cbind(a1.train, a2.train, a3.train, a4.train, a5.train, a6.train, a7.train)
  
  y.train = rdirichlet(nrow(A),A)
  
  colnames(y.train) = c("y1","y2","y3", "y4", "y5", "y6", "y7")
  
  # - test data
  
  x.test = matrix(runif(p * n, 0,1), n)
  x.test = data.frame(x.test)
  
  a1.test = exp(2.5*x.test[,1] - x.test[,2] + 5*x.test[,3]) 
  a2.test = exp(2*x.test[,1] + 2*x.test[,2] - 4*x.test[,3])
  a3.test = exp(1.5*x.test[,1] -  1.5*x.test[,2] + x.test[,3])
  a4.test = exp(3*x.test[,1] - 2*x.test[,2] + -1*x.test[,3]) 
  a5.test = exp(0.5*x.test[,1] + 2.5*x.test[,2] - 1.5*x.test[,3])
  a6.test = exp(1.5*x.test[,1] - 3*x.test[,2] + -1*x.test[,3])
  a7.test = exp(2*x.test[,1] + 2*x.test[,2] + 2*x.test[,3])
  
  A = cbind(a1.test, a2.test, a3.test, a4.test, a5.test, a6.test, a7.test)
  
  y.test = rdirichlet(nrow(A),A)
  
  colnames(y.test) = c("y1","y2","y3", "y4", "y5", "y6", "y7")
  
  # - Model
  
  septDR = glmboostLSS(y.train ~ ., data = x.train, families = DirichletSV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')
  
  cv25 = cv(model.weights(septDR), type = "subsampling")
  cvr = cvrisk(septDR, folds = cv25, grid = 1:500)
  
  StopIT = mstop(cvr)
  
  rm(septDR)
  
  septDR = glmboostLSS(y.train ~ ., data = x.train, families = DirichletSV(), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
  
  mstop.septDR <-  vector('list')
  mstop.septDR$mstop <- StopIT
  mstop.septDR$a1 <- septDR$alpha1$mstop()
  mstop.septDR$a2 <- septDR$alpha2$mstop()
  mstop.septDR$a3 <- septDR$alpha3$mstop()
  mstop.septDR$a4 <- septDR$alpha4$mstop()
  mstop.septDR$a5 <- septDR$alpha5$mstop()
  mstop.septDR$a6 <- septDR$alpha6$mstop()
  mstop.septDR$a7 <- septDR$alpha7$mstop()
  
  coef.septDR = coef(septDR, which = "")
  coef.septDR_ow = coef(septDR)
  
  # # - Dirichlet Regression
  # 
  # ydr = DR_data(y.train)
  # DR = DirichReg(ydr ~ ., data = x.train)
  # 
  
  # - Variable Selection
  
  selectedVar.a1 = vector("list")
  selectedVar.a2 = vector("list")
  selectedVar.a3 = vector("list")
  
  true.positive.a1 = vector("list")
  true.positive.a2 = vector("list")
  true.positive.a3 = vector("list")
  
  false.positive.a1 = vector('list')
  false.positive.a2 = vector('list')
  false.positive.a3 = vector('list')
  
  
  # nameVar = names(x.train)[1:p]
  # trueVar = nameVar[1:3]
  # falseVar = nameVar[4:p]
  # 
  # selectedVar.a1 = names(coef(septDR$alpha1))[-1]
  # selectedVar.a2 = names(coef(septDR$alpha2))[-1]
  # selectedVar.a3 = names(coef(septDR$alpha3))[-1]
  # 
  # true.positive.a1 = length(which(trueVar %in% names(coef(septDR$alpha1))))
  # true.positive.a2 = length(which(trueVar %in% names(coef(septDR$alpha2))))
  # true.positive.a3 = length(which(trueVar %in% names(coef(septDR$alpha3))))
  # 
  # false.positive.a1 = length(which(falseVar %in% names(coef(septDR$alpha1))))
  # false.positive.a2 = length(which(falseVar %in% names(coef(septDR$alpha2))))
  # false.positive.a3 = length(which(falseVar %in% names(coef(septDR$alpha3))))
  
  
  nameVar = names(x.train)[1:p]
  trueVar = nameVar[1:3]
  falseVar = nameVar[4:p]
  
  selectedVar.a1 = names(coef(septDR$alpha1))[-1]
  selectedVar.a2 = names(coef(septDR$alpha2))[-1]
  selectedVar.a3 = names(coef(septDR$alpha3))[-1]
  selectedVar.a4 = names(coef(septDR$alpha4))[-1]
  selectedVar.a5 = names(coef(septDR$alpha5))[-1]
  selectedVar.a6 = names(coef(septDR$alpha6))[-1]
  selectedVar.a7 = names(coef(septDR$alpha7))[-1]
  
  true.positive.a1 = length(which(trueVar %in% names(coef(septDR$alpha1))))
  true.positive.a2 = length(which(trueVar %in% names(coef(septDR$alpha2))))
  true.positive.a3 = length(which(trueVar %in% names(coef(septDR$alpha3))))
  true.positive.a4 = length(which(trueVar %in% names(coef(septDR$alpha4))))
  true.positive.a5 = length(which(trueVar %in% names(coef(septDR$alpha5))))
  true.positive.a6 = length(which(trueVar %in% names(coef(septDR$alpha6))))
  true.positive.a7 = length(which(trueVar %in% names(coef(septDR$alpha7))))
  
  false.positive.a1 = length(which(falseVar %in% names(coef(septDR$alpha1))))
  false.positive.a2 = length(which(falseVar %in% names(coef(septDR$alpha2))))
  false.positive.a3 = length(which(falseVar %in% names(coef(septDR$alpha3))))
  false.positive.a4 = length(which(falseVar %in% names(coef(septDR$alpha4))))
  false.positive.a5 = length(which(falseVar %in% names(coef(septDR$alpha5))))
  false.positive.a6 = length(which(falseVar %in% names(coef(septDR$alpha6))))
  false.positive.a7 = length(which(falseVar %in% names(coef(septDR$alpha7))))
  
  
  TPR = vector("list")
  FDR = vector ("list")
  
  TPR$alpha1 = true.positive.a1 / length(trueVar)
  TPR$alpha2 = true.positive.a2 / length(trueVar)
  TPR$alpha3 = true.positive.a3 / length(trueVar)
  TPR$alpha4 = true.positive.a4 / length(trueVar)
  TPR$alpha5 = true.positive.a5 / length(trueVar)
  TPR$alpha6 = true.positive.a6 / length(trueVar)
  TPR$alpha7 = true.positive.a7 / length(trueVar)
  

  FDR$alpha1 = false.positive.a1 / length(selectedVar.a1)
  FDR$alpha2 = false.positive.a2 / length(selectedVar.a2)
  FDR$alpha3 = false.positive.a3 / length(selectedVar.a3)
  FDR$alpha4 = false.positive.a4 / length(selectedVar.a4)
  FDR$alpha5 = false.positive.a5 / length(selectedVar.a5)
  FDR$alpha6 = false.positive.a6 / length(selectedVar.a6)
  FDR$alpha7 = false.positive.a7 / length(selectedVar.a7)
  
  # # - Predictive Performance
  # 
  # pred.a1 = predict(septDR$alpha1, newdata = x.test, type = "response")
  # pred.a2 = predict(septDR$alpha2, newdata = x.test, type = "response")
  # pred.a3 = predict(septDR$alpha3, newdata = x.test, type = "response")
  # pred.a4 = predict(septDR$alpha4, newdata = x.test, type = "response")
  # pred.a5 = predict(septDR$alpha5, newdata = x.test, type = "response")
  # pred.a6 = predict(septDR$alpha6, newdata = x.test, type = "response")
  # pred.a7 = predict(septDR$alpha7, newdata = x.test, type = "response")
  # pred.A = cbind(pred.a1, pred.a2, pred.a3, pred.a4, pred.a5, pred.a6, pred.a7)
  # pred.mu = pred.A / rowSums(pred.A)
  # 
  # pred.DR = predict(DR, newdata = x.test, mu = TRUE)
  # 
  # # MSEP
  # 
  # MSEPB = vector("list")
  # MSEPDR = vector("list")
  # 
  # MSEPB$alpha1 = mean((pred.mu[,1] - y.test[,1])**2)
  # MSEPB$alpha2 = mean((pred.mu[,2] - y.test[,2])**2)
  # MSEPB$alpha3 = mean((pred.mu[,3] - y.test[,3])**2)
  # MSEPB$alpha4 = mean((pred.mu[,4] - y.test[,4])**2)
  # MSEPB$alpha5 = mean((pred.mu[,5] - y.test[,5])**2)
  # MSEPB$alpha6 = mean((pred.mu[,6] - y.test[,6])**2)
  # MSEPB$alpha7 = mean((pred.mu[,7] - y.test[,7])**2)
  # 
  # MSEPDR$alpha1 = mean((pred.DR[,1] - y.test[,1])**2)
  # MSEPDR$alpha2 = mean((pred.DR[,2] - y.test[,2])**2)
  # MSEPDR$alpha3 = mean((pred.DR[,3] - y.test[,3])**2)
  # MSEPDR$alpha4 = mean((pred.DR[,4] - y.test[,4])**2)
  # MSEPDR$alpha5 = mean((pred.DR[,5] - y.test[,5])**2)
  # MSEPDR$alpha6 = mean((pred.DR[,6] - y.test[,6])**2)
  # MSEPDR$alpha7 = mean((pred.DR[,7] - y.test[,7])**2)
  # 
  # # NLL
  # 
  # pred.DR = predict(DR, newdata = x.test, mu = FALSE, alpha = TRUE)
  # 
  # NLL = vector("list")
  # 
  # loss = function(alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, y) {
  #   y7 = y[,7]
  #   y6 = y[,6]
  #   y5 = y[,5]
  #   y4 = y[,4]
  #   y3 = y[,3]
  #   y2 = y[,2]
  #   y1 = y[,1]
  #   
  #   - (lgamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3) + lgamma(alpha4) + lgamma(alpha5) + lgamma(alpha6) + lgamma(alpha7))
  #      + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
  #   
  # }
  # 
  # NLL$Boosting = sum(loss(alpha1 = pred.a1, alpha2 = pred.a2, alpha3 = pred.a3, alpha4 = pred.a4, alpha5 = pred.a5, alpha6 = pred.a6, alpha7 = pred.a7, y = y.test))
  # NLL$DirichReg = sum(loss(alpha1 = pred.DR[,1], alpha2 = pred.DR[,2], alpha3 = pred.DR[,3], alpha4 = pred.DR[,4], alpha5 = pred.DR[,5], alpha6 = pred.DR[,6], alpha7 = pred.DR[,7], y = y.test))
  # 
  
  return(list(TrueBeta = TrueBeta, MStop = mstop.septDR, Coefficients = coef.septDR_ow, Coefficients_plt = coef.septDR,
              TPR = TPR, FDR = FDR))
}


n = 150
p = 300
cur = 100

results = mclapply(1:cur, sim, n = n, p = p)

save(results, file = "res5")

load("res4")

# Calculating mean TPR and FDR

FP = data.frame(alpha1 = numeric(),
                alpha2 = numeric(),
                alpha3 = numeric(),
                alpha4 = numeric(),
                alpha5 = numeric(),
                alpha6 = numeric(),
                alpha7 = numeric())

TP = data.frame(alpha1 = numeric(),
                alpha2 = numeric(),
                alpha3 = numeric(),
                alpha4 = numeric(),
                alpha5 = numeric(),
                alpha6 = numeric(),
                alpha7 = numeric())

for(i in 1:cur){ 
  
  FP[i,]  = unlist(results[[i]]$FDR)
  TP[i,]  = unlist(results[[i]]$TPR)
  
}

colMeans(TP)
colMeans(FP)

# Plotting FDR and TPR

TP = TP[,-1]
FP = FP[,-1]

tp <- gather(TP, key = "alpha", value = "TP_value")
fp <- gather(FP, key = "alpha", value = "FP_value")

df <- merge(tp, fp, by = c("alpha"))

ggplot(df, aes(x = FP_value, y = TP_value)) +
  geom_point() +
  facet_grid(cols = vars(alpha)) +
  ylim(0,1) +
  xlim(0,1) +
  labs(x = "FDR", y = "TPR") +
  theme_light()


# Plotting estimates vs true values

a1 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a2 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a3 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a4 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a5 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a6 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())

a7 = data.frame(x1 = numeric(),
                x2 = numeric(),
                x3 = numeric(),
                x4 = numeric(),
                x5 = numeric(),
                x6 = numeric(),
                x7 = numeric(),
                x8 = numeric(),
                x9 = numeric(),
                x10 = numeric())



for (i in 1:cur){
  a1[i,] = results[[i]]$Coefficients_plt$alpha1[-1]
  a2[i,] = results[[i]]$Coefficients_plt$alpha2[-1]
  a3[i,] = results[[i]]$Coefficients_plt$alpha3[-1]
  a4[i,] = results[[i]]$Coefficients_plt$alpha4[-1]
  a5[i,] = results[[i]]$Coefficients_plt$alpha5[-1]
  a6[i,] = results[[i]]$Coefficients_plt$alpha6[-1]
  a7[i,] = results[[i]]$Coefficients_plt$alpha7[-1]
}

coeflist = list()
coeflist$alpha1 = a1
coeflist$alpha2 = a2
coeflist$alpha3 = a3
coeflist$alpha4 = a4
coeflist$alpha5 = a5
coeflist$alpha6 = a6
coeflist$alpha7 = a7

coef_df <- do.call(rbind, lapply(coeflist, as.data.frame))

coef_df$model = rep(names(coeflist), each = nrow(coef_df)/length(coeflist))
coef_melted = melt(coef_df, id.vars = "model")

TBet = data.frame(results[[1]]$TrueBeta)
TBet = TBet[1:10,]
TBet$variable = unique(coef_melted$variable)
true.df = gather(TBet, model, value, -variable)

ggplot(coef_melted, aes(x = variable, y = value)) + 
  geom_boxplot() +
  #geom_point(data = true.df, aes(x = variable, y = value), color = "red", size = 2.5) +
  geom_boxplot(data = true.df, aes(x = variable, y = value), color = "red") +
  facet_grid(rows = vars(model), scales = "free_y") +
  ylab("") +
  theme_light() +
  theme(strip.text = element_text(color = "black"))


# Predictive Performance

NLL = data.frame(Boosting = numeric(), DirichReg = numeric())

for (i in 1:cur) {
  NLL[i,] = unlist(results[[i]]$Likelihood, use.names = TRUE)
}

colMeans(NLL)

MSEPBoost = data.frame(alpha1 = numeric(),alpha2 = numeric(),alpha3 = numeric(), alpha4 = numeric(), alpha5 = numeric(), alpha6 = numeric(), alpha7 = numeric())
MSEPDirig = data.frame(alpha1 = numeric(),alpha2 = numeric(),alpha3 = numeric(), alpha4 = numeric(), alpha5 = numeric(), alpha6 = numeric(), alpha7 = numeric())

for (i in 1:cur) {
  MSEPBoost[i,] = unlist(results[[i]]$MSEPB)
  MSEPDirig[i,] = unlist(results[[i]]$MSEPDR)
  
}

colMeans(MSEPBoost)
colMeans(MSEPDirig)


