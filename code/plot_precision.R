p1 <- 1e3
p2 <- 1e3

set.seed(777)

theta1 <- rnorm(p1)
alpha2 <- rnorm(p2)

s1 <- 0.5
s2 <- 2

X1 <- rnorm(p1, theta1, s1)
X2 <- rnorm(p2, alpha2, s2)

z1 <- X1 / s1
z2 <- X2 / s2

P1 <- -pnorm(abs(z1)) * 2
P2 <- -pnorm(abs(z2)) * 2


group.ind <- factor(c(rep("Group 1", p1), rep("Group 2", p2)), levels = c("Group 2", "Group 1"))

pdf("~/GitHub/Bi-ASH/output/powerparity_hist_signal.pdf", width = 4, height = 4)
signal.dat <- data.frame(ind = group.ind, dat = c(theta1, alpha2))
lattice::histogram(~ dat | ind, data = signal.dat, layout = c(1, 2), xlab = "True signal", main = expression(paste("(a): Histograms of true signals")))
dev.off()

pdf("~/GitHub/Bi-ASH/output/powerparity_hist_zscore.pdf", width = 4, height = 4)
zscore.dat <- data.frame(ind = group.ind, dat = c(z1, z2))
lattice::histogram(~ dat | ind, data = zscore.dat, layout = c(1, 2), xlab = "z-score", main = expression(paste("(b): Histograms of z-scores")))
dev.off()

pdf("~/GitHub/Bi-ASH/output/powerparity_hist_pvalue.pdf", width = 4, height = 4)
pvalue.dat <- data.frame(ind = group.ind, dat = c(P1, P2))
lattice::histogram(~ dat | ind, data = pvalue.dat, layout = c(1, 2), xlab = "p-value", main = expression(paste("(c): Histograms of p-values")))
dev.off()