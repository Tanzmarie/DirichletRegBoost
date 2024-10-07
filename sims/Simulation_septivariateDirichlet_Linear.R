#########################################
library("batchtools")
library("gamboostLSS")
library("scoringRules")
library("DirichletReg")
library("tidyverse")
library("reshape2")
library("vroom")
library("Cairo")
library("latex2exp")
#########################################

source("families/Dirichlet.R")


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
  
  TrueBeta$a1 = a1
  TrueBeta$a2 = a2
  TrueBeta$a3 = a3
  TrueBeta$a4 = a4
  TrueBeta$a5 = a5
  TrueBeta$a6 = a6
  TrueBeta$a7 = a7
  
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
  
  septDR = glmboostLSS(y.train ~ ., data = x.train, families = Dirichlet(K = 7), control = boost_control(trace = TRUE, mstop = 750, nu = 0.1), method = 'noncyclic')
  
  cv25 = cv(model.weights(septDR), type = "subsampling")
  cvr = cvrisk(septDR, folds = cv25, grid = 1:750)
  
  StopIT = mstop(cvr)
  
  rm(septDR)
  
  septDR = glmboostLSS(y.train ~ ., data = x.train, families = Dirichlet(K = 7), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
  
  mstop.septDR =  vector('list')
  mstop.septDR$mstop = StopIT
  mstop.septDR$a1 = septDR$a1$mstop()
  mstop.septDR$a2 = septDR$a2$mstop()
  mstop.septDR$a3 = septDR$a3$mstop()
  mstop.septDR$a4 = septDR$a4$mstop()
  mstop.septDR$a5 = septDR$a5$mstop()
  mstop.septDR$a6 = septDR$a6$mstop()
  mstop.septDR$a7 = septDR$a7$mstop()
  
  coef.septDR = coef(septDR, which = "")
  coef.septDR_ow = coef(septDR)
  
  # # - Dirichlet Regression

  ydr = DR_data(y.train, trafo = FALSE)
  DR = DirichReg(ydr ~ ., data = x.train)


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
  
  
  nameVar = names(x.train)[1:p]
  trueVar = nameVar[1:3]
  falseVar = nameVar[4:p]
  
  selectedVar.a1 = names(coef(septDR$a1))[-1]
  selectedVar.a2 = names(coef(septDR$a2))[-1]
  selectedVar.a3 = names(coef(septDR$a3))[-1]
  selectedVar.a4 = names(coef(septDR$a4))[-1]
  selectedVar.a5 = names(coef(septDR$a5))[-1]
  selectedVar.a6 = names(coef(septDR$a6))[-1]
  selectedVar.a7 = names(coef(septDR$a7))[-1]
  
  true.positive.a1 = length(which(trueVar %in% names(coef(septDR$a1))))
  true.positive.a2 = length(which(trueVar %in% names(coef(septDR$a2))))
  true.positive.a3 = length(which(trueVar %in% names(coef(septDR$a3))))
  true.positive.a4 = length(which(trueVar %in% names(coef(septDR$a4))))
  true.positive.a5 = length(which(trueVar %in% names(coef(septDR$a5))))
  true.positive.a6 = length(which(trueVar %in% names(coef(septDR$a6))))
  true.positive.a7 = length(which(trueVar %in% names(coef(septDR$a7))))
  
  false.positive.a1 = length(which(falseVar %in% names(coef(septDR$a1))))
  false.positive.a2 = length(which(falseVar %in% names(coef(septDR$a2))))
  false.positive.a3 = length(which(falseVar %in% names(coef(septDR$a3))))
  false.positive.a4 = length(which(falseVar %in% names(coef(septDR$a4))))
  false.positive.a5 = length(which(falseVar %in% names(coef(septDR$a5))))
  false.positive.a6 = length(which(falseVar %in% names(coef(septDR$a6))))
  false.positive.a7 = length(which(falseVar %in% names(coef(septDR$a7))))
  
  
  TPR = vector("list")
  TNR = vector("list")
  FDR = vector("list")
  PPV = vector("list")
  NPV = vector("list")
  
  
  TPR$a1 = true.positive.a1 / length(trueVar)
  TPR$a2 = true.positive.a2 / length(trueVar)
  TPR$a3 = true.positive.a3 / length(trueVar)
  TPR$a4 = true.positive.a4 / length(trueVar)
  TPR$a5 = true.positive.a5 / length(trueVar)
  TPR$a6 = true.positive.a6 / length(trueVar)
  TPR$a7 = true.positive.a7 / length(trueVar)
  
  TNR$a1 = 1 - (false.positive.a1 / length(falseVar))
  TNR$a2 = 1 - (false.positive.a2 / length(falseVar))
  TNR$a3 = 1 - (false.positive.a3 / length(falseVar))
  TNR$a4 = 1 - (false.positive.a4 / length(falseVar))
  TNR$a5 = 1 - (false.positive.a5 / length(falseVar))
  TNR$a6 = 1 - (false.positive.a6 / length(falseVar))
  TNR$a7 = 1 - (false.positive.a7 / length(falseVar))

  FDR$a1 = false.positive.a1 / length(selectedVar.a1)
  FDR$a2 = false.positive.a2 / length(selectedVar.a2)
  FDR$a3 = false.positive.a3 / length(selectedVar.a3)
  FDR$a4 = false.positive.a4 / length(selectedVar.a4)
  FDR$a5 = false.positive.a5 / length(selectedVar.a5)
  FDR$a6 = false.positive.a6 / length(selectedVar.a6)
  FDR$a7 = false.positive.a7 / length(selectedVar.a7)
  
  PPV$a1 = 1 - (false.positive.a1 / length(selectedVar.a1))
  PPV$a2 = 1 - (false.positive.a2 / length(selectedVar.a2))
  PPV$a3 = 1 - (false.positive.a3 / length(selectedVar.a3))
  PPV$a4 = 1 - (false.positive.a4 / length(selectedVar.a4))
  PPV$a5 = 1 - (false.positive.a5 / length(selectedVar.a5))
  PPV$a6 = 1 - (false.positive.a6 / length(selectedVar.a6))
  PPV$a7 = 1 - (false.positive.a7 / length(selectedVar.a7))
  
  NPV$a1 = ((1 - (false.positive.a1 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a1))
  NPV$a2 = ((1 - (false.positive.a2 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a2))
  NPV$a3 = ((1 - (false.positive.a3 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a3))
  NPV$a4 = ((1 - (false.positive.a4 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a4))
  NPV$a5 = ((1 - (false.positive.a5 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a5))
  NPV$a6 = ((1 - (false.positive.a6 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a6))
  NPV$a7 = ((1 - (false.positive.a7 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a7))
  
  # - Predictive Performance

  pred.a1 = predict(septDR$a1, newdata = x.test, type = "response")
  pred.a2 = predict(septDR$a2, newdata = x.test, type = "response")
  pred.a3 = predict(septDR$a3, newdata = x.test, type = "response")
  pred.a4 = predict(septDR$a4, newdata = x.test, type = "response")
  pred.a5 = predict(septDR$a5, newdata = x.test, type = "response")
  pred.a6 = predict(septDR$a6, newdata = x.test, type = "response")
  pred.a7 = predict(septDR$a7, newdata = x.test, type = "response")
  pred.A = cbind(pred.a1, pred.a2, pred.a3, pred.a4, pred.a5, pred.a6, pred.a7)
  pred.mu = pred.A / rowSums(pred.A)

  pred.DR = predict(DR, newdata = x.test, mu = TRUE)

  # MSEP

  MSEPB = vector("list")
  MSEPDR = vector("list")

  MSEPB$a1 = mean((pred.mu[,1] - y.test[,1])**2)
  MSEPB$a2 = mean((pred.mu[,2] - y.test[,2])**2)
  MSEPB$a3 = mean((pred.mu[,3] - y.test[,3])**2)
  MSEPB$a4 = mean((pred.mu[,4] - y.test[,4])**2)
  MSEPB$a5 = mean((pred.mu[,5] - y.test[,5])**2)
  MSEPB$a6 = mean((pred.mu[,6] - y.test[,6])**2)
  MSEPB$a7 = mean((pred.mu[,7] - y.test[,7])**2)

  MSEPDR$a1 = mean((pred.DR[,1] - y.test[,1])**2)
  MSEPDR$a2 = mean((pred.DR[,2] - y.test[,2])**2)
  MSEPDR$a3 = mean((pred.DR[,3] - y.test[,3])**2)
  MSEPDR$a4 = mean((pred.DR[,4] - y.test[,4])**2)
  MSEPDR$a5 = mean((pred.DR[,5] - y.test[,5])**2)
  MSEPDR$a6 = mean((pred.DR[,6] - y.test[,6])**2)
  MSEPDR$a7 = mean((pred.DR[,7] - y.test[,7])**2)

  # NLL

  pred.DR = predict(DR, newdata = x.test, mu = FALSE, a = TRUE)

  NLL = vector("list")

  loss = function(a1, a2, a3, a4, a5, a6, a7, y) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]

    - (lgamma(a1 + a2 + a3 + a4 + a5 + a6 + a7) - (lgamma(a1) + lgamma(a2) + lgamma(a3) + lgamma(a4) + lgamma(a5) + lgamma(a6) + lgamma(a7))
       + ((a1 - 1) * log(y1) + (a2 - 1) * log(y2) + (a3 - 1) * log(y3) + (a4 - 1) * log(y4) + (a5 - 1) * log(y5) + (a6 - 1) * log(y6) + (a7 - 1) * log(y7)))

  }

  NLL$Boosting = sum(loss(a1 = pred.a1, a2 = pred.a2, a3 = pred.a3, a4 = pred.a4, a5 = pred.a5, a6 = pred.a6, a7 = pred.a7, y = y.test))
  NLL$DirichReg = sum(loss(a1 = pred.DR[,1], a2 = pred.DR[,2], a3 = pred.DR[,3], a4 = pred.DR[,4], a5 = pred.DR[,5], a6 = pred.DR[,6], a7 = pred.DR[,7], y = y.test))

  
  return(list(TrueBeta = TrueBeta,
              NLL = NLL,
              MSEPB = MSEPB,
              MSEPDR = MSEPDR,
              MStop = mstop.septDR,
              Coefficients = coef.septDR_ow,
              Coefficients_plt = coef.septDR,
              TPR = TPR, FDR = FDR, TNR = TNR, PPV = PPV, NPV = NPV))
}


n = 150
p = 10
cur = 50

set.seed(123)

# Generate a list of reproducible seeds
seeds = sample.int(1e6, cur)

# Run simulations using different seeds
results = lapply(1:cur, function(i) sim(seed = seeds[i], n = n, p = p))

# saveRDS(results, file = "49Dir50.RData")

# results = readRDS("49Dir50.RData")

# Calculating mean TPR and FDR
FalDis = data.frame(a1 = numeric(),
                    a2 = numeric(),
                    a3 = numeric(),
                    a4 = numeric(),
                    a5 = numeric(),
                    a6 = numeric(),
                    a7 = numeric())

TruePos = data.frame(a1 = numeric(),
                     a2 = numeric(),
                     a3 = numeric(),
                     a4 = numeric(),
                     a5 = numeric(),
                     a6 = numeric(),
                     a7 = numeric())

TrueNeg = data.frame(a1 = numeric(),
                     a2 = numeric(),
                     a3 = numeric(),
                     a4 = numeric(),
                     a5 = numeric(),
                     a6 = numeric(),
                     a7 = numeric())

PosPred = data.frame(a1 = numeric(),
                     a2 = numeric(),
                     a3 = numeric(),
                     a4 = numeric(),
                     a5 = numeric(),
                     a6 = numeric(),
                     a7 = numeric())

NegPred = data.frame(a1 = numeric(),
                     a2 = numeric(),
                     a3 = numeric(),
                     a4 = numeric(),
                     a5 = numeric(),
                     a6 = numeric(),
                     a7 = numeric())



for(i in 1:cur){ 
  
  FalDis[i,]  = unlist(results[[i]]$FDR)
  TruePos[i,]  = unlist(results[[i]]$TPR)
  TrueNeg[i,]  = unlist(results[[i]]$TNR)
  PosPred [i,]  = unlist(results[[i]]$PPV)
  NegPred[i,]  = unlist(results[[i]]$NPV)
  
}

colMeans(TruePos) * 100
colMeans(TrueNeg) * 100
colMeans(FalDis) * 100
colMeans(PosPred) * 100
colMeans(NegPred) * 100

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
  a1[i,] = results[[i]]$Coefficients_plt$a1[-1]
  a2[i,] = results[[i]]$Coefficients_plt$a2[-1]
  a3[i,] = results[[i]]$Coefficients_plt$a3[-1]
  a4[i,] = results[[i]]$Coefficients_plt$a4[-1]
  a5[i,] = results[[i]]$Coefficients_plt$a5[-1]
  a6[i,] = results[[i]]$Coefficients_plt$a6[-1]
  a7[i,] = results[[i]]$Coefficients_plt$a7[-1]
}

coeflist = list()
coeflist$a1 = a1
coeflist$a2 = a2
coeflist$a3 = a3
coeflist$a4 = a4
coeflist$a5 = a5
coeflist$a6 = a6
coeflist$a7 = a7

coef_df = do.call(rbind, lapply(coeflist, as.data.frame))

coef_df$model = rep(names(coeflist), each = nrow(coef_df)/length(coeflist))
coef_df = coef_df[,c(1:4,11)]
coef_melted = melt(coef_df, id.vars = "model")

TBet = data.frame(results[[1]]$TrueBeta)
TBet = TBet[1:4,]
TBet$variable = unique(coef_melted$variable)[1:4]
true.df = gather(TBet, model, value, -variable)

custom_labels = c(a1 = "\u03b11",  
                  a2 = "\u03b12",  
                  a3 = "\u03b13",
                  a4 = "\u03b14",
                  a5 = "\u03b15",
                  a6 = "\u03b16",
                  a7 = "\u03b17")

CairoPDF("test.pdf", width = 9, height = 5)
ggplot(coef_melted, aes(x = variable, y = value)) + 
  geom_boxplot() +
  geom_boxplot(data = true.df, aes(x = variable, y = value), color = "red") +
  facet_wrap(~ model, scales = "free_y", labeller = labeller(model = custom_labels), nrow = 4) +
  labs(title = "", y = "Coefficient", x = "Variable") +
  scale_x_discrete(labels = function(x) ifelse(x == "x4", TeX("x\\{4,...,P\\}"), x)) +
  theme_bw()
dev.off()


# Predictive Performance

NLL = data.frame(Boosting = numeric(), DirichReg = numeric())

for (i in 1:cur) {
  NLL[i,] = unlist(results[[i]]$NLL, use.names = TRUE)
}

colMeans(NLL)

MSEPBoost = data.frame(a1 = numeric(),a2 = numeric(),a3 = numeric(), a4 = numeric(), a5 = numeric(), a6 = numeric(), a7 = numeric())
MSEPDirig = data.frame(a1 = numeric(),a2 = numeric(),a3 = numeric(), a4 = numeric(), a5 = numeric(), a6 = numeric(), a7 = numeric())

for (i in 1:cur) {
  MSEPBoost[i,] = unlist(results[[i]]$MSEPB)
  MSEPDirig[i,] = unlist(results[[i]]$MSEPDR)
  
}

sqrt(colMeans(MSEPBoost))
sqrt(colMeans(MSEPDirig))



