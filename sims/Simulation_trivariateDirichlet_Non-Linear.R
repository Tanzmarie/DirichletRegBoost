#########################################
library("batchtools")
library("gamboostLSS")
library("scoringRules")
library("DirichletReg")
library("tidyverse")
library("reshape2")
#########################################

source("families/trivariateDirichlet.R")

sim = function(seed, n.train, p){
  
  
  ###### data generation
  # - training data
  
  n = n.train + 1000
  weight.mstop = c(rep(1, times = n.train),rep(0, times = 1000))
  
  x.train = matrix(rnorm(p * n, 0,1), n)
  x.train = data.frame(x.train)
  
  a1.train = exp(x.train[,1]**2) 
  a2.train = exp(2*tanh(3*x.train[,2]))
  a3.train = exp(cos(x.train[,3]))
  
  
  A = cbind(a1.train,a2.train,a3.train)
  
  y.train = rdirichlet(nrow(A),A)
  
  
  # - Model
  
  x = paste(c(paste("bbs(X", 1:p, ")", sep = "")), collapse = "+")
  form = as.formula((paste("y.train ~",  x)))
  
  trivDRNON = gamboostLSS(form, data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, risk = "oobag", nu = 0.1), weights = weight.mstop, method = 'noncyclic')
  #trivDRNON = gamboostLSS(form, data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 1000, nu = 0.1), method = 'noncyclic')
  
  
  # cv25 = cv(model.weights(trivDRNON), type = "kfold")
  # cvr = cvrisk(trivDRNON, folds = cv25, grid = 1:1000)
  # StopIT = mstop(cvr)
  
  StopIT = which.min(risk(trivDRNON,merge = T))
  
  x.train = x.train[weight.mstop == 1, ]
  y.train = y.train[weight.mstop == 1, ]
  
  trivDRNON = gamboostLSS(form, data = x.train, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
  
  mstop.trivDRNON =  vector('list')
  mstop.trivDRNON$mstop = StopIT
  mstop.trivDRNON$a1 = trivDRNON$alpha1$mstop()
  mstop.trivDRNON$a2 = trivDRNON$alpha2$mstop()
  mstop.trivDRNON$a3 = trivDRNON$alpha3$mstop()
  
  coef.trivDRNON = coef(trivDRNON, which = "")
  coef.trivDRNON_ow = coef(trivDRNON)
  
  # - Variable Selection 
  
  nameVar = paste(c(paste("bbs(X",1:p,")", sep = "")))
  
  selectedVar.a1 = names(coef(trivDRNON$alpha1))
  selectedVar.a2 = names(coef(trivDRNON$alpha2))
  selectedVar.a3 = names(coef(trivDRNON$alpha3))
  
  true.positive.a1 = length(which(nameVar[1] %in% names(coef(trivDRNON$alpha1))))
  true.positive.a2 = length(which(nameVar[2] %in% names(coef(trivDRNON$alpha2))))
  true.positive.a3 = length(which(nameVar[3] %in% names(coef(trivDRNON$alpha3))))
  
  false.positive.a1 = length(which(nameVar[-1] %in% names(coef(trivDRNON$alpha1))))
  false.positive.a2 = length(which(nameVar[-2] %in% names(coef(trivDRNON$alpha2))))
  false.positive.a3 = length(which(nameVar[-3] %in% names(coef(trivDRNON$alpha3))))
  
  
  TPR = vector("list")
  FDR = vector ("list")
  
  TPR$Overall = sum(true.positive.a1 + true.positive.a2 + true.positive.a3) / 3
  TPR$alpha1 = true.positive.a1 / length(nameVar[1])
  TPR$alpha2 = true.positive.a2 / length(nameVar[2])
  TPR$alpha3 = true.positive.a3 / length(nameVar[3])
  
  FDR$Overall = sum(false.positive.a1 + false.positive.a2 + false.positive.a3) / (length(selectedVar.a1) + length(selectedVar.a2) + length(selectedVar.a3))
  FDR$alpha1 = false.positive.a1 / length(selectedVar.a1)
  FDR$alpha2 = false.positive.a2 / length(selectedVar.a2)
  FDR$alpha3 = false.positive.a3 / length(selectedVar.a3)
  
  # - Effect Plotting
  
  plotDR = vector("list")
  
  newX = matrix(rep(seq(from = -3, to = 3,length.out = 500),times = p),nrow =500)
  
  plotDR$newX_alpha1 = predict(trivDRNON, data.frame(newX), parameter = 'alpha1',which = 1, type = 'link')
  plotDR$newX_alpha2 = predict(trivDRNON, data.frame(newX), parameter = 'alpha2',which = 2, type = 'link')
  plotDR$newX_alpha3 = predict(trivDRNON, data.frame(newX), parameter = 'alpha3',which = 3, type = 'link')
  
  n.ges = n
  
  return(list(n = n.ges, MStop = mstop.trivDRNON, Coefficients = coef.trivDRNON_ow, Coefficients_plt = coef.trivDRNON,
              TPR = TPR, FDR = FDR, plt = plotDR, newdata = newX ))
}


n.train = 1000
p = 10
cur = 100

results = mclapply(1:cur, sim, n.train = n.train, p = p)

save(results, file = "NONLIN2")

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

# Effect Plotting

df.a1 = data.frame(x = results[[1]]$newdata[,1],trueVal = results[[1]]$newdata[,1]**2)
df.a2 = data.frame(x = results[[1]]$newdata[,2],trueVal = 2*tanh(3*results[[1]]$newdata[,2]))
df.a3 = data.frame(x = results[[1]]$newdata[,3],trueVal = cos(results[[1]]$newdata[,3]))

for (i in 1:cur) {
df.a1 = cbind(df.a1,results[[i]]$plt$newX_alpha1)
df.a2 = cbind(df.a2,results[[i]]$plt$newX_alpha2)
df.a3 = cbind(df.a3,results[[i]]$plt$newX_alpha3)
}

names(df.a1)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))
names(df.a2)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))
names(df.a3)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))

df.combined = bind_rows(list(alpha1 = df.a1, alpha2 = df.a2, alpha3 = df.a3),.id = "alpha")

df.long = df.combined %>% pivot_longer(cols = starts_with("y"),names_to = "var", values_to = "value")


ggplot(df.long, aes(x, value)) +
  geom_line(aes(group = var), color = "black", linewidth = 1, alpha = 0.1) +
  geom_line(data = unique(df.long[,1:3]),
            aes(x = x, y = trueVal), color = "red", linewidth = 2, linetype = "dotted") +
  facet_grid(~ alpha, scales = "free_x") +
  labs(x = "x", y = "f(x)", title = "") +
  theme_light() +
  theme(strip.text = element_text(color = "black"))

