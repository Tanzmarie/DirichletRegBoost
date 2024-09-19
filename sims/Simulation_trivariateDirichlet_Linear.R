#########################################
library("batchtools")
library("gamboostLSS")
library("scoringRules")
library("DirichletReg")
library("tidyverse")
library("reshape2")
library("Cairo")
#########################################

source("families/Dirichlet.R")


sim = function(seed,n,p) {
  
  set.seed(seed)
  
  a1 = c(2.5,-1,3, rep(0,p-3))
  a2 = c(0,0,0,2,2,-1 ,rep(0,p-6))
  a3 = c(0,0,0,0,0,0,1.5,-1.5,1 ,rep(0,p-9))
  
  TrueBeta =  vector('list')
  
  TrueBeta$a1 = a1
  TrueBeta$a2 = a2
  TrueBeta$a3 = a3
  
  
  ###### data generation
  # - training data
  
  x.train = matrix(runif(p * n, 0,1), n)
  x.train = data.frame(x.train)
  
  a1.train = exp(2.5*x.train[,1] - x.train[,2] + 3*x.train[,3]) 
  a2.train = exp(2*x.train[,4] + 2*x.train[,5] - x.train[,6])
  a3.train = exp(1.5*x.train[,7] -  1.5*x.train[,8] + x.train[,9])
  
  A = cbind(a1.train,a2.train,a3.train)
  
  y.train = rdirichlet(nrow(A),A)
  
  colnames(y.train) = c("y1","y2","y3")
  
  # - test data
  
  x.test = matrix(runif(p * n, 0,1), n)
  x.test = data.frame(x.test)
  
  a1.test = exp(2.5*x.test[,1] - x.test[,2] + 3*x.test[,3]) 
  a2.test = exp(2*x.test[,4] + 2*x.test[,5] - x.test[,6])
  a3.test = exp(1.5*x.test[,7] -  1.5*x.test[,8] + x.test[,9])
  
  
  A = cbind(a1.test,a2.test,a3.test)
  
  y.test = rdirichlet(nrow(A),A)
  
  colnames(y.test) = c("y1","y2","y3")
  
  # - Model
  
  trivDR = glmboostLSS(y.train ~ ., data = x.train, families = Dirichlet(K = 3), control = boost_control(trace = TRUE, mstop = 500, nu = 0.1), method = 'noncyclic')
  
  cvr = cvrisk(trivDR, folds = cv(model.weights(trivDR), type = "kfold"), grid = 1:500)
  
  StopIT = mstop(cvr)
  
  rm(trivDR)
  
  trivDR = glmboostLSS(y.train ~ ., data = x.train, families = Dirichlet(K = 3), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')

  mstop.trivDR =  vector('list')
  mstop.trivDR$mstop = StopIT
  mstop.trivDR$a1 = trivDR$a1$mstop()
  mstop.trivDR$a2 = trivDR$a2$mstop()
  mstop.trivDR$a3 = trivDR$a3$mstop()
  
  coef.trivDR = coef(trivDR, which = "")
  coef.trivDR_ow = coef(trivDR)
  
  # - Dirichlet Regression
  # 
  # ydr = DR_data(y.train)
  # DR = DirichReg(ydr ~ ., data = x.train)
  # 
  # 
  # - Variable Selection
  
  nameVar = names(x.train)[1:p]
  trueVar = nameVar[1:3]
  falseVar = nameVar[4:p]

  selectedVar.a1 = names(coef(trivDR$a1))[-1]
  selectedVar.a2 = names(coef(trivDR$a2))[-1]
  selectedVar.a3 = names(coef(trivDR$a3))[-1]

  true.positive.a1 = length(which(nameVar[1:3] %in% names(coef(trivDR$a1))))
  true.positive.a2 = length(which(nameVar[4:6] %in% names(coef(trivDR$a2))))
  true.positive.a3 = length(which(nameVar[7:9] %in% names(coef(trivDR$a3))))

  false.positive.a1 = length(which(nameVar[4:p] %in% names(coef(trivDR$a1))))
  false.positive.a2 = length(which(nameVar[c(1:3,7:p)] %in% names(coef(trivDR$a2))))
  false.positive.a3 = length(which(nameVar[c(1:6,10:p)] %in% names(coef(trivDR$a3))))

  
  TPR = vector("list")
  TNR = vector("list")
  FDR = vector("list")
  PPV = vector("list")
  NPV = vector("list")
  
 
  TPR$a1 = true.positive.a1 / length(nameVar[1:3])
  TPR$a2 = true.positive.a2 / length(nameVar[4:6])
  TPR$a3 = true.positive.a3 / length(nameVar[7:9])
  
  TNR$a1 = 1 - (false.positive.a1 / length(falseVar))
  TNR$a2 = 1 - (false.positive.a2 / length(falseVar))
  TNR$a3 = 1 - (false.positive.a3 / length(falseVar))
  
  FDR$a1 = false.positive.a1 / length(selectedVar.a1)
  FDR$a2 = false.positive.a2 / length(selectedVar.a2)
  FDR$a3 = false.positive.a3 / length(selectedVar.a3)
  
  PPV$a1 = 1 - (false.positive.a1 / length(selectedVar.a1))
  PPV$a2 = 1 - (false.positive.a2 / length(selectedVar.a2))
  PPV$a3 = 1 - (false.positive.a3 / length(selectedVar.a3))
  
  NPV$a1 = ((1 - (false.positive.a1 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a1))
  NPV$a2 = ((1 - (false.positive.a2 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a2))
  NPV$a3 = ((1 - (false.positive.a3 / length(falseVar))) * length(falseVar)) / (length(nameVar) - length(selectedVar.a3))
  
  # - Predictive Performance

  pred.a1 = predict(trivDR$a1, newdata = x.test, type = "response")
  pred.a2 = predict(trivDR$a2, newdata = x.test, type = "response")
  pred.a3 = predict(trivDR$a3, newdata = x.test, type = "response")
  pred.A = cbind(pred.a1,pred.a2,pred.a3)
  pred.mu = pred.A / rowSums(pred.A)

  # pred.DR = predict(DR, newdata = x.test, mu = TRUE)

  # MSEP

  MSEPB = vector("list")
  MSEPDR = vector("list")

  MSEPB$a1 = mean((pred.mu[,1] - y.test[,1])**2)
  MSEPB$a2 = mean((pred.mu[,2] - y.test[,2])**2)
  MSEPB$a3 = mean((pred.mu[,3] - y.test[,3])**2)

  # MSEPDR$a1 = mean((pred.DR[,1] - y.test[,1])**2)
  # MSEPDR$a2 = mean((pred.DR[,2] - y.test[,2])**2)
  # MSEPDR$a3 = mean((pred.DR[,3] - y.test[,3])**2)

  # NLL

  # pred.DR = predict(DR, newdata = x.test, mu = FALSE, alpha = TRUE)

  NLL = vector("list")

  loss = function(a1, a2, a3, y) {
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]

    - (lgamma(a1 + a2 + a3) - (lgamma(a1) + lgamma(a2) + lgamma(a3))
       + ((a1 - 1) * log(y1) + (a2 - 1) * log(y2) + (a3 - 1) * log(y3)))

  }

  NLL$Boosting = sum(loss(a1 = pred.a1, a2 = pred.a2, a3 = pred.a3, y = y.test))
  # NLL$DirichReg = sum(loss(a1 = pred.DR[,1], a2 = pred.DR[,2], a3 = pred.DR[,3], y = y.test))

  return(list(TrueBeta = TrueBeta, MStop = mstop.trivDR, Coefficients = coef.trivDR_ow, Coefficients_plt = coef.trivDR,
              Likelihood = NLL, MSEPB = MSEPB, # MSEPDR = MSEPDR,
              TPR = TPR, FDR = FDR, TNR = TNR, PPV = PPV, NPV = NPV))
}

n = 150
p = 300
cur = 50
  
set.seed(123)

# Generate a list of reproducible seeds
seeds = sample.int(1e6, cur)

# Run simulations using different seeds
results = lapply(1:cur, function(i) sim(seed = seeds[i], n = n, p = p))

# Performance Criteria

FalDis = data.frame(a1 = numeric(),
                a2 = numeric(),
                a3 = numeric())

TruePos = data.frame(a1 = numeric(),
                a2 = numeric(),
                a3 = numeric())

TrueNeg = data.frame(a1 = numeric(),
                a2 = numeric(),
                a3 = numeric())

PosPred = data.frame(a1 = numeric(),
                a2 = numeric(),
                a3 = numeric())

NegPred = data.frame(a1 = numeric(),
                a2 = numeric(),
                a3 = numeric())



for(i in 1:cur){ 
  FalDis[i,]  = unlist(results[[i]]$FDR)
  TruePos[i,]  = unlist(results[[i]]$TPR)
  TrueNeg[i,]  = unlist(results[[i]]$TNR)
  PosPred [i,]  = unlist(results[[i]]$PPV)
  NegPred[i,]  = unlist(results[[i]]$NPV)

}

colMeans(FalDis)
colMeans(TruePos)
colMeans(TrueNeg)
colMeans(PosPred)
colMeans(NegPred)


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


for (i in 1:cur){
  a1[i,] = results[[i]]$Coefficients_plt$a1[-1]
  a2[i,] = results[[i]]$Coefficients_plt$a2[-1]
  a3[i,] = results[[i]]$Coefficients_plt$a3[-1]
}

coeflist = list()
coeflist$a1 = a1
coeflist$a2 = a2
coeflist$a3 = a3

coef_df = do.call(rbind, lapply(coeflist, as.data.frame))

coef_df$model = rep(names(coeflist), each = nrow(coef_df)/length(coeflist))
coef_melted = melt(coef_df, id.vars = "model")

TBet = data.frame(results[[1]]$TrueBeta)
TBet = TBet[1:10,]
TBet$variable = unique(coef_melted$variable)
true.df = gather(TBet, model, value, -variable)

custom_labels = c(a1 = "\u03b11",  
                   a2 = "\u03b12",  
                   a3 = "\u03b13")

CairoPDF("test.pdf", width = 7, height = 5)
ggplot(coef_melted, aes(x = variable, y = value)) + 
  geom_boxplot() +
  geom_boxplot(data = true.df, aes(x = variable, y = value), color = "red") +
  facet_grid(rows = vars(model), scales = "free_y", labeller = labeller(model = custom_labels)) +
  labs(
    title = "",
    y = "", x = "Variable"
  ) +
  theme_bw()
dev.off()



# Predictive Performance

NLL = data.frame(Boosting = numeric(), DirichReg = numeric())

for (i in 1:cur) {
  NLL[i,] = unlist(results[[i]]$Likelihood, use.names = TRUE)
}

colMeans(NLL)

MSEPBoost = data.frame(a1 = numeric(),a2 = numeric(),a3 = numeric())
MSEPDirig = data.frame(a1 = numeric(),a2 = numeric(),a3 = numeric())

for (i in 1:cur) {
  MSEPBoost[i,] = unlist(results[[i]]$MSEPB)
  MSEPDirig[i,] = unlist(results[[i]]$MSEPDR)

}

colMeans(MSEPBoost)
colMeans(MSEPDirig)


