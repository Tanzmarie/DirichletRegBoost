#########################################
library("batchtools")
library("gamboostLSS")
library("scoringRules")
library("DirichletReg")
library("tidyverse")
library("reshape2")
#########################################

source("families/trivariateDirichlet.R")


sim = function(seed,n,p) {
  
  set.seed(seed)
  
  a1 = c(2.5,-1,3, rep(0,p-3))
  a2 = c(0,0,0,2,2,-1 ,rep(0,p-6))
  a3 = c(0,0,0,0,0,0,1.5,-1.5,1 ,rep(0,p-9))
  
  TrueBeta =  vector('list')
  
  TrueBeta$alpha1 = a1
  TrueBeta$alpha2 = a2
  TrueBeta$alpha3 = a3
  
  
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
  
  trivDR = glmboostLSS(y.train ~ ., data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')
  
  cv25 = cv(model.weights(trivDR), type = "kfold")
  cvr = cvrisk(trivDR, folds = cv25, grid = 1:1000)
  
  StopIT = mstop(cvr)
  
  rm(trivDR)
  
  trivDR = glmboostLSS(y.train ~ ., data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')

  mstop.trivDR <-  vector('list')
  mstop.trivDR$mstop <- StopIT
  mstop.trivDR$a1 <- trivDR$alpha1$mstop()
  mstop.trivDR$a2 <- trivDR$alpha2$mstop()
  mstop.trivDR$a3 <- trivDR$alpha3$mstop()
  
  coef.trivDR = coef(trivDR, which = "")
  coef.trivDR_ow = coef(trivDR)
  
  # - Dirichlet Regression
  
  ydr = DR_data(y.train)
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
  
  
  # nameVar = names(x.train)[1:p]
  # trueVar = nameVar[1:3]
  # falseVar = nameVar[4:p]
  # 
  # selectedVar.a1 = names(coef(trivDR$alpha1))[-1]
  # selectedVar.a2 = names(coef(trivDR$alpha2))[-1]
  # selectedVar.a3 = names(coef(trivDR$alpha3))[-1]
  # 
  # true.positive.a1 = length(which(trueVar %in% names(coef(trivDR$alpha1))))
  # true.positive.a2 = length(which(trueVar %in% names(coef(trivDR$alpha2))))
  # true.positive.a3 = length(which(trueVar %in% names(coef(trivDR$alpha3))))
  # 
  # false.positive.a1 = length(which(falseVar %in% names(coef(trivDR$alpha1))))
  # false.positive.a2 = length(which(falseVar %in% names(coef(trivDR$alpha2))))
  # false.positive.a3 = length(which(falseVar %in% names(coef(trivDR$alpha3))))
   
  
  nameVar = names(x.train)[1:p]
  trueVar = nameVar[1:3]
  falseVar = nameVar[4:p]

  selectedVar.a1 = names(coef(trivDR$alpha1))[-1]
  selectedVar.a2 = names(coef(trivDR$alpha2))[-1]
  selectedVar.a3 = names(coef(trivDR$alpha3))[-1]

  true.positive.a1 = length(which(nameVar[1:3] %in% names(coef(trivDR$alpha1))))
  true.positive.a2 = length(which(nameVar[4:6] %in% names(coef(trivDR$alpha2))))
  true.positive.a3 = length(which(nameVar[7:9] %in% names(coef(trivDR$alpha3))))

  false.positive.a1 = length(which(nameVar[4:p] %in% names(coef(trivDR$alpha1))))
  false.positive.a2 = length(which(nameVar[c(1:3,7:p)] %in% names(coef(trivDR$alpha2))))
  false.positive.a3 = length(which(nameVar[c(1:6,10:p)] %in% names(coef(trivDR$alpha3))))

  
  
  TPR = vector("list")
  FDR = vector ("list")
  
  TPR$Overall = sum(true.positive.a1 + true.positive.a2 + true.positive.a3) / (length(nameVar[1:3]) * 3)
  TPR$alpha1 = true.positive.a1 / length(nameVar[1:3])
  TPR$alpha2 = true.positive.a2 / length(nameVar[4:6])
  TPR$alpha3 = true.positive.a3 / length(nameVar[7:9])
  
  FDR$Overall = sum(false.positive.a1 + false.positive.a2 + false.positive.a3) / (length(selectedVar.a1) + length(selectedVar.a2) + length(selectedVar.a3))
  FDR$alpha1 = false.positive.a1 / length(selectedVar.a1)
  FDR$alpha2 = false.positive.a2 / length(selectedVar.a2)
  FDR$alpha3 = false.positive.a3 / length(selectedVar.a3)
  
  
  # - Predictive Performance

  pred.a1 = predict(trivDR$alpha1, newdata = x.test, type = "response")
  pred.a2 = predict(trivDR$alpha2, newdata = x.test, type = "response")
  pred.a3 = predict(trivDR$alpha3, newdata = x.test, type = "response")
  pred.A = cbind(pred.a1,pred.a2,pred.a3)
  pred.mu = pred.A / rowSums(pred.A)

  pred.DR = predict(DR, newdata = x.test, mu = TRUE)

  # MSEP

  MSEPB = vector("list")
  MSEPDR = vector("list")

  MSEPB$alpha1 = mean((pred.mu[,1] - y.test[,1])**2)
  MSEPB$alpha2 = mean((pred.mu[,2] - y.test[,2])**2)
  MSEPB$alpha3 = mean((pred.mu[,3] - y.test[,3])**2)

  MSEPDR$alpha1 = mean((pred.DR[,1] - y.test[,1])**2)
  MSEPDR$alpha2 = mean((pred.DR[,2] - y.test[,2])**2)
  MSEPDR$alpha3 = mean((pred.DR[,3] - y.test[,3])**2)

  # NLL

  pred.DR = predict(DR, newdata = x.test, mu = FALSE, alpha = TRUE)

  NLL = vector("list")

  loss = function(alpha1, alpha2, alpha3, y) {
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]

    - (lgamma(alpha1 + alpha2 + alpha3) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3)))

  }

  NLL$Boosting = sum(loss(alpha1 = pred.a1, alpha2 = pred.a2, alpha3 = pred.a3, y = y.test))
  NLL$DirichReg = sum(loss(alpha1 = pred.DR[,1], alpha2 = pred.DR[,2], alpha3 = pred.DR[,3], y = y.test))

  # # Energy Score
  # 
  # es_boost = vector()
  # es_DR = vector()
  # 
  # for (i in 1:length(pred.a1)) {
  # 
  #   pred_sample_boost = matrix(NA, nrow = 3, ncol = 10000)
  #   pred_sample_DR = matrix(NA, nrow = 3, ncol = 10000)
  # 
  #   # Boosting approach
  # 
  #   sample_boost = rdirichlet(10000, pred.A[i,])
  # 
  #   pred_sample_boost[1,] = sample_boost[,1]
  #   pred_sample_boost[2,] = sample_boost[,2]
  #   pred_sample_boost[3,] = sample_boost[,3]
  # 
  #   es_boost[i] <- es_sample(y = c(y.test[i,1], y.test[i,2], y.test[i,3]), dat = pred_sample_boost)
  # 
  #   # DirichReg
  # 
  #   sample_DR = rdirichlet(10000, pred.DR[i,])
  # 
  #   pred_sample_DR[1,] = sample_DR[,1]
  #   pred_sample_DR[2,] = sample_DR[,2]
  #   pred_sample_DR[3,] = sample_DR[,3]
  # 
  #   es_DR[i] <- es_sample(y = c(y.test[i,1], y.test[i,2], y.test[i,3]), dat = pred_sample_DR)
  # 
  # 
  # }
  # 
  # energy_score = list()
  # energy_score$Boosting = mean(es_boost)
  # energy_score$DirichReg = mean(es_DR)
  # energy_score
  
  n.ges = n
  
  return(list(TrueBeta = TrueBeta, n = n.ges, MStop = mstop.trivDR, Coefficients = coef.trivDR_ow, Coefficients_plt = coef.trivDR,
              Likelihood = NLL, MSEPB = MSEPB, MSEPDR = MSEPDR, TPR = TPR, FDR = FDR))
}


n = 150
p = 49
cur = 100

results = mclapply(1:cur, sim, n = n, p = p)

save(results, file = "results2")

# Calculating mean TPR and FDR

FP = data.frame(Overall = numeric(),
                  alpha1 = numeric(),
                  alpha2 = numeric(),
                  alpha3 = numeric())

TP = data.frame(Overall = numeric(),
                 alpha1 = numeric(),
                 alpha2 = numeric(),
                 alpha3 = numeric())

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


for (i in 1:cur){
  a1[i,] = results[[i]]$Coefficients_plt$alpha1[-1]
  a2[i,] = results[[i]]$Coefficients_plt$alpha2[-1]
  a3[i,] = results[[i]]$Coefficients_plt$alpha3[-1]
}

coeflist = list()
coeflist$alpha1 = a1
coeflist$alpha2 = a2
coeflist$alpha3 = a3

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
  facet_grid(rows = vars(model)) +
  theme_light()

# Predictive Performance

NLL = data.frame(Boosting = numeric(), DirichReg = numeric())

for (i in 1:cur) {
  NLL[i,] = unlist(results[[i]]$Likelihood, use.names = TRUE)
}

colMeans(NLL)

MSEPBoost = data.frame(alpha1 = numeric(),alpha2 = numeric(),alpha3 = numeric())
MSEPDirig = data.frame(alpha1 = numeric(),alpha2 = numeric(),alpha3 = numeric())

for (i in 1:cur) {
  MSEPBoost[i,] = unlist(results[[i]]$MSEPB)
  MSEPDirig[i,] = unlist(results[[i]]$MSEPDR)

}

colMeans(MSEPBoost)
colMeans(MSEPDirig)


