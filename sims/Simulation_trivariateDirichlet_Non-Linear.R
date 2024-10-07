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

sim = function(seed, n.train, p){
  
  set.seed(seed)
  
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
  
  trivDRNON = gamboostLSS(form, data = x.train, families = Dirichlet(K = 3), control = boost_control(trace = TRUE, mstop = 750, risk = "oobag", nu = 0.1), weights = weight.mstop, method = 'noncyclic')
  
  StopIT = which.min(risk(trivDRNON,merge = T))
  
  x.train = x.train[weight.mstop == 1, ]
  y.train = y.train[weight.mstop == 1, ]
  
  trivDRNON = gamboostLSS(form, data = x.train, families = Dirichlet(K = 3), control = boost_control(trace = TRUE, mstop = StopIT, nu = 0.1), method = 'noncyclic')
  
  mstop.trivDRNON =  vector('list')
  mstop.trivDRNON$mstop = StopIT
  mstop.trivDRNON$a1 = trivDRNON$a1$mstop()
  mstop.trivDRNON$a2 = trivDRNON$a2$mstop()
  mstop.trivDRNON$a3 = trivDRNON$a3$mstop()
  
  coef.trivDRNON = coef(trivDRNON, which = "")
  coef.trivDRNON_ow = coef(trivDRNON)
  
  # - Variable Selection 
  
  nameVar = paste(c(paste("bbs(X",1:p,")", sep = "")))
  falseVar = paste(c(paste("bbs(X",2:p,")", sep = "")))
  
  selectedVar.a1 = names(coef(trivDRNON$a1))
  selectedVar.a2 = names(coef(trivDRNON$a2))
  selectedVar.a3 = names(coef(trivDRNON$a3))
  
  true.positive.a1 = length(which(nameVar[1] %in% names(coef(trivDRNON$a1))))
  true.positive.a2 = length(which(nameVar[2] %in% names(coef(trivDRNON$a2))))
  true.positive.a3 = length(which(nameVar[3] %in% names(coef(trivDRNON$a3))))
  
  false.positive.a1 = length(which(nameVar[-1] %in% names(coef(trivDRNON$a1))))
  false.positive.a2 = length(which(nameVar[-2] %in% names(coef(trivDRNON$a2))))
  false.positive.a3 = length(which(nameVar[-3] %in% names(coef(trivDRNON$a3))))
  
  
  TPR = vector("list")
  TNR = vector("list")
  FDR = vector("list")
  PPV = vector("list")
  NPV = vector("list")
  
  TPR$a1 = true.positive.a1 / length(nameVar[1])
  TPR$a2 = true.positive.a2 / length(nameVar[2])
  TPR$a3 = true.positive.a3 / length(nameVar[3])
  
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
  
  
  # - Effect Plotting
  
  plotDR = vector("list")
  
  newX = matrix(rep(seq(from = -3, to = 3,length.out = 500),times = p),nrow =500)
  
  plotDR$newX_a1 = predict(trivDRNON, data.frame(newX), parameter = 'a1',which = 1, type = 'link')
  plotDR$newX_a2 = predict(trivDRNON, data.frame(newX), parameter = 'a2',which = 2, type = 'link')
  plotDR$newX_a3 = predict(trivDRNON, data.frame(newX), parameter = 'a3',which = 3, type = 'link')
  
  n.ges = n
  
  return(list(n = n.ges,
              MStop = mstop.trivDRNON,
              Coefficients = coef.trivDRNON_ow,
              Coefficients_plt = coef.trivDRNON,
              TPR = TPR, FDR = FDR, TNR = TNR, PPV = PPV, NPV = NPV,
              plt = plotDR,
              newdata = newX ))
}


n.train = 1000
p = 10
cur = 100

set.seed(123)

# Generate a list of reproducible seeds
seeds = sample.int(1e6, cur)

# Run simulations using different seeds
results = lapply(1:cur, function(i) sim(seed = seeds[i], n.train = n.train, p = p))

# saveRDS(results, file = "1000NonDir100.RData")

# Load or save results
results = readRDS("1000NonDir100.RData")


# Calculating mean TPR and FDR
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

colMeans(TruePos) * 100
colMeans(TrueNeg) * 100
colMeans(FalDis) * 100
colMeans(PosPred) * 100
colMeans(NegPred) * 100

# Effect Plotting

df.a1 = data.frame(x = results[[1]]$newdata[,1],trueVal = results[[1]]$newdata[,1]**2)
df.a2 = data.frame(x = results[[1]]$newdata[,2],trueVal = 2*tanh(3*results[[1]]$newdata[,2]))
df.a3 = data.frame(x = results[[1]]$newdata[,3],trueVal = cos(results[[1]]$newdata[,3]))

for (i in 1:cur) {
df.a1 = cbind(df.a1,results[[i]]$plt$newX_a1)
df.a2 = cbind(df.a2,results[[i]]$plt$newX_a2)
df.a3 = cbind(df.a3,results[[i]]$plt$newX_a3)
}

names(df.a1)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))
names(df.a2)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))
names(df.a3)[3:(cur + 2)] = paste(c(paste("y",1:cur,"", sep = "")))

df.combined = bind_rows(list(a1 = df.a1, a2 = df.a2, a3 = df.a3),.id = "a")

df.long = df.combined %>% pivot_longer(cols = starts_with("y"),names_to = "var", values_to = "value")


custom_labels = c(a1 = "\u03b11",  
                  a2 = "\u03b12",  
                  a3 = "\u03b13")

CairoPDF("test.pdf", width = 9, height = 5)
ggplot(df.long, aes(x, value)) +
  geom_line(aes(group = var), color = "black", linewidth = 1, alpha = 0.1) +
  geom_line(data = unique(df.long[,1:3]),
            aes(x = x, y = trueVal), color = "red", linewidth = 1.5, linetype = "dotted") +
  facet_grid(~ a, scales = "free_x", labeller = labeller(a = custom_labels)) +
  labs(x = "x", y = "f(x)", title = "") +
  theme_bw() +
  theme(strip.text = element_text(color = "black"))
dev.off()

