set.seed(777)

theta0 <- rnorm(1e2)
beta0 <- rnorm(10)
beta1 <- rnorm(10, 0, 2)
theta1 <- beta0 + beta1


pdf("~/GitHub/Bi-ASH/output/model_enrichment.pdf", width = 5, height = 4)
par(mar = c(5.1, 4.1, 1, 1))

plot(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01)), type = "l", lty = 1, lwd = 2, xlab = "signal", ylab = "density")
lines(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01), 0, sqrt(1 + 2^2)), type = "l", lty = 1, col = "blue", lwd = 2)
points(theta0, rep(0, 1e2), pch = 20)
points(theta1, rep(0.01, 10), pch = 13, col = "blue")
legend("topleft", pch = c(13, 20), col = c("blue", "black"), c("gene set", "other genes"), bty = "n")
text(x = c(0, 0), y = c(0.38, 0.15), c(expression("f"), expression("h")), col = c("black", "blue"))

dev.off()