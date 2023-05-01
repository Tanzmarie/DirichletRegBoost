#Load dependencies and data

library("DirichletReg")
library("tidyverse")
library("stargazer")

data = read.table("D:/Universität/Master/Masterarbeit/DirichletRegBoost/Application/data_party_replication.tab",header =T, sep = "\t")

#combine election number with country of election
data$idvar <- paste(data$country_iso,data$election_number)


df = data.frame()

for (i in unique(data$idvar)){
  subing = data[data$idvar == i,]
  
  subing2 = subing %>% group_by(parfam) %>% group_by(parfam) %>% summarise(sum(t1_vote), first(wevent), first(wpart_all), first(lpopw_part_all), first(duration), first(imf), first(bailout), first(crisis_election), first(unemployment_change), first(debt_change), first(gdp_change), first(east), first(idvar))
  subing2 = as.data.frame(subing2)
  
  df = rbind(df,subing2)
  
}

datn = reshape(df, timevar = "parfam", idvar = "first(idvar)",  direction = "wide")

#Finding NA's in the votes for nonexistent parties and assigning them to zero
datn[which(is.na(datn[,"sum(t1_vote).2"])), "sum(t1_vote).2"] <- 0
datn[which(is.na(datn[,"sum(t1_vote).3"])), "sum(t1_vote).3"] <- 0
datn[which(is.na(datn[,"sum(t1_vote).4"])), "sum(t1_vote).4"] <- 0
datn[which(is.na(datn[,"sum(t1_vote).5"])), "sum(t1_vote).5"] <- 0
datn[which(is.na(datn[,"sum(t1_vote).6"])), "sum(t1_vote).6"] <- 0
datn[which(is.na(datn[,"sum(t1_vote).7"])), "sum(t1_vote).7"] <- 0

#dividing in order to have the interval [0,1] for proportions instead of percents
datn[,c("sum(t1_vote).1",
        "sum(t1_vote).2",
        "sum(t1_vote).3",
        "sum(t1_vote).4",
        "sum(t1_vote).5",
        "sum(t1_vote).6",
        "sum(t1_vote).7")] <- datn[,c("sum(t1_vote).1",
                                      "sum(t1_vote).2",
                                      "sum(t1_vote).3",
                                      "sum(t1_vote).4",
                                      "sum(t1_vote).5",
                                      "sum(t1_vote).6",
                                      "sum(t1_vote).7")]/100


#creating dirichlet response
Vote = datn[,c("sum(t1_vote).1",
                       "sum(t1_vote).2",
                       "sum(t1_vote).3",
                       "sum(t1_vote).4",
                       "sum(t1_vote).5",
                       "sum(t1_vote).6",
                       "sum(t1_vote).7")]

names(Vote) = c("Social Democratic", "Conservative / Christian democratic", "Liberals", "Greens", "Radical Left", "Radical Rights", "Others")

Vote = DR_data(Vote)


#changing column names
colnames(datn)[1:13] = c("idvar","t_vote","wevent","wpart_all","lpopw_part_all","duration","imf","bailout","crisis","unemployment","debt","gdp","east")
colnames(datn)[3:13] = c("Protests", "Participants","ScaledParticipants", "Duration", "IMF", "Bailout", "Crisis", "Unemployment", "Debt", "GDP", "East")

#omitting nas in wpart_all
datn[which(is.na(datn[,"Participants"])), "Participants"] = 0

#checking for correlation
cor(datn[,3:13])

#scalling participants for numerical stability, prob not necessary since well use scaled events
datn$ScaledParticipants = datn$ScaledParticipants/1000
datn$Participants = datn$Participants/1000

#scaling events by duration

datn$Protests = datn$Protests/datn$Duration

stargazer(datn[,c(3,8:13)])

rm(df,subing,subing2,i)

############ Some plots

#plotting proportions against protests

plot(rep(datn$Protests, 7),
     as.numeric(Vote),
     pch = 21,
     bg = rep(c("red", "black", "yellow", "green", "darkorchid4","blue","darkgrey"), each = 39),
     xlab = "Protests", ylab = "Proportion", ylim = 0:1)

#plotting the proportions in singular plots

par(mfrow = c(2,4))

for (i in 1:7){
plot(datn$Protests,
     as.numeric(Vote[,i]),
     pch = 21,
     main = as.character(i),
     xlab = "Protests",
     ylab = "Proportion",
     ylim = 0:1)
}

###############################################################################################################
# 0:1: Usual Dirichlet Regression via DiriletReg Package
# Variable Selection is performed via the AIC
# We start with the common parametrization and then move on the alternative for better interpretations

# First, try to catch the effects of the Protests on the Proportion only via DirigReg:

margs = DirichReg(Vote ~ Protests, data = datn)

Xnew = data.frame(Protests = seq(min(datn$Protests), max(datn$Protests), length.out = 117))
pred = predict(object = margs, newdata = Xnew)


#Plotting the predicted curves

par(mfrow = c(2,4))

for (i in 1:7){
  plot(datn$Protests,
       as.numeric(Vote[,i]),
       pch = 21,
       main = as.character(i),
       xlab = "Protests",
       ylab = "Proportion",
       ylim = 0:1)
  
  lines(Xnew$Protests, pred[,i], col = "red", lwd = 2)
}

# Only the predictions

for (i in 1:7){
  plot(Xnew$Protests, pred[, i], type = "l",
       main = "Predicted Proportion of Vote 1 vs Protests",
       xlab = "Protests", ylab = "Proportion",
       ylim = c(0, 1))
}

#All in one plot

x11()

plot(rep(datn$Protests, 7),
     as.numeric(Vote),
     pch = 21,
     bg = rep(c("red", "black", "yellow", "green", "darkorchid4","blue","darkgrey"), each = 39),
     xlab = "Protests", ylab = "Proportion", ylim = 0:1)


for (i in 1:7){
  lines(cbind(Xnew$Protests, predict(margs, Xnew)[, i]), col = c("red", "black", "yellow", "green", "darkorchid4","blue","darkgrey")[i], lwd = 2)
}


# Searching for the best model according to AIC
variables <- c("Protests", "Unemployment", "Debt", "GDP", "Crisis", "Bailout", "East")

results = list()


for (i in 1:length(variables)){
  combinations = combn(variables, i, simplify = FALSE)
  
  for (j in 1:length(combinations)){
    
    tr = paste("Vote ~ ", paste(combinations[[j]], collapse = "+"))
    
    results[[tr]] = DirichReg(eval(as.formula(tr)), data = datn)
  }
}


#calculate aic for every model
aic_list = sapply(results, function(x) AIC(x))

#find the minimum
min_index <- which.min(aic_list)

# extract corresponding model formula
best_formula <- names(results)[min_index]

#fit the best model according to aic
mod = DirichReg(eval(as.formula(best_formula)), data = datn)

summary(mod)


################################# Lets repeat for alternative parametriziation

variables <- c("wevents", "unemployment", "debt", "gdp", "crisis", "bailout", "east")

results = list()


for (i in 1:length(variables)){
  combinations = combn(variables, i, simplify = FALSE)
  
  for (j in 1:length(combinations)){
    
    tr = paste("Vote ~ ", paste(combinations[[j]], collapse = "+"), "|", paste(paste(combinations[[j]], collapse = "+")))
    
    results[[tr]] = DirichReg(eval(as.formula(tr)) , data = datn, model = "alternative")
  }
}

#calculate aic for every model
aic_list = sapply(results, function(x) AIC(x))

#find the minimum
min_index <- which.min(aic_list)

# extract corresponding model formula
best_formula <- names(results)[min_index]

#fit the best model according to aic
mod = DirichReg(eval(as.formula(best_formula)), data = datn, model = "alternative")

summary(mod)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#Boosting with the common parametriziation

source("families/septivariateDirichlet.R")

boostmod = glmboostLSS(Vote ~ Protests + Unemployment + Debt + GDP + Crisis + Bailout + East, data = datn,  families = DirichletSV(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL), control = boost_control(trace = TRUE, mstop = 500, nu  = 0.1), method = 'noncyclic')
coef(boostmod, off2int = TRUE)
plot(boostmod)




