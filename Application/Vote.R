#Load dependencies and data
library("DirichletReg")
library("tidyverse")
library("stargazer")
library("gamboostLSS")

data = read.table("D:/Universität/Master/Masterarbeit/DirichletRegBoost/Application/data_party_replication.tab",header =T, sep = "\t")

#combine election number with country of election
data$idvar = paste(data$country_iso,data$election_number)



df = data.frame()

for (i in unique(data$idvar)){
  subing = data[data$idvar == i,]
  
  subing2 = subing %>% group_by(parfam) %>% group_by(parfam) %>% summarise(sum(t1_vote), first(wevent), first(wpart_all), first(lpopw_part_all), first(duration), first(imf), first(bailout), first(crisis_election), first(unemployment_change), first(debt_change), first(gdp_change), first(east), first(idvar))
  subing2 = as.data.frame(subing2)
  
  df = rbind(df,subing2)
  
}

datn = reshape(df, timevar = "parfam", idvar = "first(idvar)",  direction = "wide")

#Finding NA's in the votes for nonexistent parties and assigning them to zero
datn[which(is.na(datn[,"sum(t1_vote).2"])), "sum(t1_vote).2"] = 0
datn[which(is.na(datn[,"sum(t1_vote).3"])), "sum(t1_vote).3"] = 0
datn[which(is.na(datn[,"sum(t1_vote).4"])), "sum(t1_vote).4"] = 0
datn[which(is.na(datn[,"sum(t1_vote).5"])), "sum(t1_vote).5"] = 0
datn[which(is.na(datn[,"sum(t1_vote).6"])), "sum(t1_vote).6"] = 0
datn[which(is.na(datn[,"sum(t1_vote).7"])), "sum(t1_vote).7"] = 0

#dividing in order to have the interval [0,1] for proportions instead of percents
datn[,c("sum(t1_vote).1",
        "sum(t1_vote).2",
        "sum(t1_vote).3",
        "sum(t1_vote).4",
        "sum(t1_vote).5",
        "sum(t1_vote).6",
        "sum(t1_vote).7")] = datn[,c("sum(t1_vote).1",
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

names(Vote) = c("Social democratic", "Conservative / Christian democratic", "Liberals", "Greens", "Radical left", "Radical populist right", "Others")

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

rm(df,subing,subing2,i)

# Boosting Dirichlet regression models with the common parametriziation

source("families/Dirichlet.R")

set.seed(100)

boostmod = glmboostLSS(Vote ~  Protests + Unemployment + Debt + GDP + Crisis + Bailout + East, data = datn,
                       families = Dirichlet(K=7),
                       control = boost_control(trace = TRUE, mstop = 1000, nu  = 0.1), method = 'noncyclic')


cv10 = cv(model.weights(boostmod), type = "subsampling")
cvr = cvrisk(boostmod, folds = cv10, grid = 1:500)
StopIT = mstop(cvr)

boostmod = glmboostLSS(Vote ~  Protests + Unemployment + Debt + GDP + Crisis + Bailout + East, data = datn,
                       families =  Dirichlet(K=7),
                       control = boost_control(trace = TRUE, mstop = StopIT, nu  = 0.1), method = 'noncyclic')
coef(boostmod, off2int = TRUE)

par(mfrow = c(2,4), mar = c(4, 4, 2, 5))

plot(boostmod)

########################################
s = stabsel(boostmod, q = 35, PFER = 6)

x11()
selected(s)
plot(s, main = "Maximum selection frequencies", type = "maxsel", np = 15)
########################################

# Plotting predictions for the Boosting approach

set.seed(100)

newX = data.frame(Protests = seq(min(datn$Protests), max(datn$Protests), length.out = 117),
                  Unemployment = seq(min(datn$Unemployment), max(datn$Unemployment), length.out = 117),
                  Debt = seq(min(datn$Debt), max(datn$Debt), length.out = 117),
                  GDP = seq(min(datn$GDP), max(datn$GDP), length.out = 117),
                  Crisis = sample(c(0,1), size = 117, replace = TRUE),
                  Bailout = sample(c(0,1), size = 117, replace = TRUE),
                  East = sample(c(0,1), size = 117, replace = TRUE))


trueProtest = newX[,1]
trueUnemployment = newX[,2]
trueDebt = newX[,3]
trueGDP = newX[,4]
trueCrisis = newX[,5]
trueBailout = newX[,6]
trueEast = newX[,7]

# Define a function to make predictions
predf = function(parameter, coeff) {
  pred = predict(boostmod, newX, parameter = parameter, which = coeff, type = "response")
  # Check if pred is a single value, repeat it if necessary
  if (length(pred) == 1) {
    pred = rep(pred, times = 117)
  }
  return(pred)
}

parameters = paste0("a", 1:7)

pred_list = lapply(parameters, coeff = 2, predf)
pred.A = do.call(cbind, pred_list)

pred.mu = pred.A / rowSums(pred.A)
pred.mu

colnames(pred.mu) =  colnames(Vote)

# Convert pred.mu to a data frame and set column names based on Vote
pred_mu = as.data.frame(pred.mu)


pred_mu$TrueProtest = trueProtest
# Reshape data to long format, excluding TrueProtest from the reshaping
pred_df = pred_mu %>%
  pivot_longer(cols = -TrueProtest, 
               names_to = "Party", 
               values_to = "Proportion")

# Set the factor levels of Party to match the order of the labels from Vote
pred_df$Party = factor(pred_df$Party, levels = colnames(Vote))

# Create the plot using ggplot2
ggplot(pred_df, aes(x = TrueProtest, y = Proportion)) +
  geom_line() +
  facet_wrap(~ Party, nrow = 2, ncol = 4) +
  ylim(0, 0.3) +
  labs(y = "Expectations of voting proportions", x = "Crisis (no/yes)") +
  #scale_x_continuous(breaks = seq(-50,100, length.out = 7)) +
  theme_bw()

# ggplot(pred_df, aes(x = factor(TrueProtest), y = Proportion)) +
#   geom_boxplot() +
#   facet_wrap(~ Party, nrow = 2, ncol = 4) +
#   ylim(0, 0.3) +
#   labs(y = "Expectations of voting proportions", x = "East (no/yes)") +
#   theme_bw()




############ Plotting Proportions

# Plotting Proportions against year for any country of choice via country_iso == ""
# 
# pres = data %>% filter(country_iso == "GR")
# pres = pres %>%
#   mutate(parfam = recode(parfam,
#                          `1` = "Social democratic",
#                          `2` = "Conservative / Christian democratic",
#                          `3` = "Liberals",
#                          `4` = "Greens",
#                          `5` = "Radical left",
#                          `6` = "Radical populist right",
#                          `7` = "Others"))
# 
# pres = pres %>%
#   mutate(parfam = factor(parfam, levels = c("Social democratic",
#                                             "Conservative / Christian democratic",
#                                             "Liberals",
#                                             "Greens",
#                                             "Radical left",
#                                             "Radical populist right",
#                                             "Others")))
# 
# pres = pres %>% select(c("parfam", "year_t", "wevent", "t_vote"))
# pres = pres %>% 
#   group_by(parfam, year_t) %>%
#   summarise(t_vote = sum(t_vote, na.rm = TRUE) / 100, wevent = first(wevent))
#  
# ggplot(pres, aes(x = year_t, y = t_vote, color = factor(parfam), group = parfam)) +
#   geom_line(size = 1) +     
#   geom_point(size = 2) +     
#   labs(x = "Year", y = "Voting proportions", color = "Party family") + 
#   theme_bw() +          
#   scale_x_continuous(breaks =  seq(min(test$year_t), max(test$year_t), by = 1)) +      
#   theme(legend.position = "bottom")+  
#   annotate("rect", xmin = 2009, xmax = 2012, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
#   annotate("text", x = 2010.5, y = max(test$t_vote), label = "1043 Protests", size = 4, vjust = 0.5, color = "black")
# 
#                       
# 
# par(mfrow = c(2,4))
# 
# for (i in 1:7){
# plot(datn$Protests,
#      as.numeric(Vote[,i]),
#      pch = 21,
#      main = as.character(i),
#      xlab = "Protests",
#      ylab = "Proportion",
#      ylim = 0:1)
# }





###############################################################################################################
# 0:1: Usual Dirichlet Regression via DiriletReg Package
# Variable Selection is performed via the AIC
# We start with the common parametrization and then move on the alternative for better interpretations

# First, try to catch the effects of the Protests on the Proportion only via DirigReg:
# 
# margs = DirichReg(Vote ~  Protests, data = datn)
# 
# Xnew = data.frame(Protests = seq(min(datn$Protests), max(datn$Protests), length.out = 117))
# pred = predict(object = margs, newdata = Xnew)
# 
# 
# #Plotting the predicted curves
# 
# par(mfrow = c(2,4))
# 
# for (i in 1:7){
#   plot(datn$Protests,
#        as.numeric(Vote[,i]),
#        pch = 21,
#        main = as.character(i),
#        xlab = "Protests",
#        ylab = "Proportion",
#        ylim = 0:1)
#   
#   lines(Xnew$Protests, pred[,i], col = "red", lwd = 2)
# }
# 
# 
# #All in one plot
# 
# x11()
# 
# plot(rep(datn$Protests, 7),
#      as.numeric(Vote),
#      pch = 21,
#      bg = rep(c("red", "black", "yellow", "green", "darkorchid4","blue","darkgrey"), each = 39),
#      xlab = "Protests", ylab = "Proportion", ylim = 0:1)
# 
# 
# for (i in 1:7){
#   lines(cbind(Xnew$Protests, predict(margs, Xnew)[, i]), col = c("red", "black", "yellow", "green", "darkorchid4","blue","darkgrey")[i], lwd = 2)
# }
# 
# 
# # Searching for the best model according to AIC
# variables = c("Protests", "Unemployment", "Debt", "GDP", "Crisis", "Bailout", "East")
# 
# results = list()
# 
# 
# for (i in 1:length(variables)){
#   combinations = combn(variables, i, simplify = FALSE)
#   
#   for (j in 1:length(combinations)){
#     
#     tr = paste("Vote ~ ", paste( combinations[[j]], collapse = "+"))
#     
#     results[[tr]] = DirichReg(eval(as.formula(tr)), data = datn)
#   }
# }
# 
# 
# #calculate aic for every model
# aic_list = sapply(results, function(x) AIC(x))
# 
# #find the minimum
# min_index = which.min(aic_list)
# 
# #extract corresponding model formula
# best_formula = names(results)[min_index]
# 
# #fit the best model according to aic
# mod = DirichReg(eval(as.formula(best_formula)), data = datn)
# 
# summary(mod)
# 
# 
# pred.DR = predict.DirichletRegModel(mod, newdata = newX)
# 
# predict(mod)
# 
# nnewX
# 
# par(mfrow = c(2,4))
# for(i in 1:7){
#   plot(x = trueProtest, y = pred.DR[,i], ylim = c(0,0.3), type = "l", ylab = "Proportion", xlab = "Protests", main = i)
# }
# 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


