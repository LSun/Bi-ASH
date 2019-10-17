set.seed(777)

theta0 <- rnorm(1e2)
beta0 <- rnorm(10)
beta1 <- rnorm(10, 0, 2)
theta1 <- beta0 + beta1


pdf("~/GitHub/Bi-ASH/output/enrich.pdf", width = 10, height = 4)
par(mfrow = c(1, 2))

plot(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01)), type = "l", lty = 3, lwd = 2, xlab = "effect", ylab = "density", main = expression("(a)"))
lines(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01), 0, sqrt(1 + 2^2)), type = "l", lty = 3, col = "blue", lwd = 2)
points(theta0, rep(0, 1e2), pch = 20)
points(theta1, rep(0.01, 10), pch = 13, col = "blue")
legend("topleft", pch = c(13, 20), col = c("blue", "black"), c("gene set", "background"), bty = "n")
text(x = c(0, 0), y = c(0.38, 0.15), c(expression("f"), expression("h")), col = c("black", "blue"))

plot(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01)), type = "l", lty = 3, lwd = 2, xlab = "effect", ylab = "density", main = expression("(b)"))
lines(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01), 0, sqrt(5)), type = "l", lty = 3, col = "blue", lwd = 2)
lines(seq(-3.5, 3.5, by = 0.01), dnorm(seq(-3.5, 3.5, by = 0.01), 0, 2), type = "l", lty = 3, col = "red", lwd = 2)
points(theta1, rep(0.02, 10), pch = 13, col = "blue")
points(beta0, rep(0, 10), pch = 1, col = "black")
points(beta1, rep(0.01, 10), pch = 4, col = "red")
legend("topleft", pch = c(13, 1, 4), col = c("blue", "black", "red"), c("gene set", "background part", "enrichment part"), bty = "n")
text(x = c(0, 0, 0), y = c(0.38, 0.15, 0.22), c(expression("f"), expression("h=g*f"), expression("g")), col = c("black", "blue", "red"))

dev.off()