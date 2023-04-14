library("DirichletReg")
require('gamboostLSS')

source("trivariateDirichlet.R")

head(ArcticLake)

ArcticLake

ArcticLake[, 1:3]

AL <- DR_data(ArcticLake[, 1:3])

#Plots

plot(AL, cex = 0.5, a2d = list(colored = FALSE, c.grid = TRUE))

plot(rep(ArcticLake$depth, 3),
     as.numeric(AL),
     pch = 21,
     bg = rep(c("#E495A5", "#86B875", "#7DB0DD"), each = 39),
     xlab = "Depth (m)", ylab = "Proportion", ylim = 0:1)

# Fitting the data via DirichReg

lake1 <- DirichReg(AL ~ depth, ArcticLake)
lake1

lake2 <- update(lake1, . ~ . + I(depth^2) | . + I(depth^2) | . + I(depth^2))
anova(lake1, lake2)

summary(lake2)

# Boosting approach
ACL = ArcticLake[, 1:3]
y1 = ACL[,1]
y2 = ACL[,2]
y3 = ACL[,3]

model1 = glmboostLSS(formula = cbind(y1,y2,y3) ~ depth + I(depth^2), data = ArcticLake, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 10000, nu = 0.1), method = 'noncyclic')
model2 = glmboostLSS(formula = cbind(y1,y2,y3) ~ depth + I(depth^2), data = ArcticLake, families = DirichletTV(), control = boost_control(trace = TRUE, mstop = 100000, nu = 0.01), method = 'cyclic')

cv10f = cv(model.weights(model1), type = "kfold")
cvr = cvrisk(model1, fold = cv10f)

plot(cvr)
cvr

coef(model1)
coef(model2)

plot(model1, parameter = "alpha3")

# Result plotting

par(mar = c(4, 4, 4, 4) + 0.1)

plot(rep(ArcticLake$depth, 3),
     as.numeric(AL), pch = 21, bg = rep(c("#E495A5", "#86B875", "#7DB0DD"), each = 39), xlab = "Depth (m)", ylab = "Proportion", ylim = 0:1, 
     main = "Sediment Composition in an Arctic Lake")

Xnew <- data.frame(depth = seq(min(ArcticLake$depth), max(ArcticLake$depth), length.out = 100))
for (i in 1:3) lines(cbind(Xnew, predict(lake2, Xnew)[, i]), col = c("#E495A5", "#86B875", 
                                                                     "#7DB0DD")[i], lwd = 2)
legend("topleft", legend = c("Sand", "Silt", "Clay"), lwd = 2, col = c("#E495A5", 
                                                                       "#86B875", "#7DB0DD"), pt.bg = c("#E495A5", "#86B875", "#7DB0DD"), pch = 21, 
       bty = "n")

par(new = TRUE)
plot(cbind(Xnew, predict(lake2, Xnew, F, F, T)), lty = "24", type = "l", ylim = c(0, 
                                                                                  max(predict(lake2, Xnew, F, F, T))), axes = F, ann = F, lwd = 2)
axis(4)
mtext(expression(paste("Precision (", phi, ")", sep = "")), 4, line = 3)
legend("top", legend = c(expression(hat(mu[c] == hat(alpha)[c]/hat(alpha)[0])), expression(hat(phi) == 
                                                                                             hat(alpha)[0])), lty = c(1, 2), lwd = c(3, 2), bty = "n")
