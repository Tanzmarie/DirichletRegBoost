source("families/Dirichlet.R")

set.seed(1234)

Z = sapply(1:4,function(x) runif(1000,-.5,.5))

A = cbind(exp(+Z[,1]+Z[,2]+Z[,3]-Z[,4]),
           exp(+Z[,1]+Z[,2]-Z[,3]-Z[,4]),
           exp(+Z[,1]-Z[,2]-Z[,3]-Z[,4]))


y = rdirichlet(nrow(A),A)

colnames(y) = paste0("y",1:3)

# - Model
y = DR_data(y, trafo = FALSE)

coef(DirichReg(y~Z))

BooMod = glmboostLSS(y ~ Z, families = Dirichlet(K=3),
                      control = boost_control(trace = TRUE, mstop = 500, nu = 0.1), method = 'noncyclic')

custom_labels <- c(a1 = "\u03b11",  
                   a2 = "\u03b12",  
                   a3 = "\u03b13")

CairoPDF("test.pdf", width = 10, height = 4)
# Set up the plotting area
par(mfrow = c(1, 3), mar = c(5, 4, 4, 6))

plot(BooMod, off2int = TRUE, col = c(1, 2, 3, 4, "orange"), parameter = "a1", ylim = c(-1.2, 1.2))
title(main = custom_labels["a1"])
abline(h = c(-1, 1), lty = 4, col = "darkgrey")

plot(BooMod, off2int = TRUE, col = c(1, 2, 3, 4, "orange"), parameter = "a2", ylim = c(-1.2, 1.2))
title(main = custom_labels["a2"])
abline(h = c(-1, 1), lty = 4, col = "darkgrey")

plot(BooMod, off2int = TRUE, col = c(1, 2, 3, 4, "orange"), parameter = "a3", ylim = c(-1.2, 1.2))
title(main = custom_labels["a3"])
abline(h = c(-1, 1), lty = 4, col = "darkgrey")

dev.off()





