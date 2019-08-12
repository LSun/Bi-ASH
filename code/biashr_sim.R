source('~/GitHub/cashr_paper/code/RNAseq_pipeline.R')
source('~/GitHub/Bi-ASH/code/biashr.R')

r <- readRDS("~/GitHub/Bi-ASH/data/liver.rds")
r <- readRDS("~/GitHub/Bi-ASH/data/muscle.rds")

set.seed(777)

num.nctrl <- 5e3
num.ctrl <- 1e3
prop.null <- 0.9
nsamp <- 10
q <- 0.1

Y <- lcpm(r)
subset <- top_genes_index(2e4, Y)
r <- r[sample(subset, num.nctrl + num.ctrl), sample(ncol(r), nsamp)]
r.nctrl <- r[1 : num.nctrl, ]
r.ctrl <- r[-(1 : num.nctrl), ]
pi0hat.biashr <- FDPq.biashr <- c()
pi0hat.ashr <- FDPq.ashr <- c()

for (i in 1 : 1000) {
r.nctrl.signaladded <- seqgendiff::thin_2group(mat = as.matrix(r.nctrl), 
                     prop_null = prop.null,
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 0.5))

theta <- r.nctrl.signaladded$coefmat

summary.stat <- count_to_summary(rbind(r.nctrl.signaladded$mat, as.matrix(r.ctrl)), r.nctrl.signaladded$designmat)

x1 <- summary.stat$X[1 : num.nctrl]
s1 <- summary.stat$s[1 : num.nctrl]
x2 <- summary.stat$X[-(1 : num.nctrl)]
s2 <- summary.stat$s[-(1 : num.nctrl)]

biashr.fit <- biashr(x1, s1, x2, s2)
ashr.fit <- ashr::ash(x1, s1)

pi0hat.biashr[i] <- biashr.fit$g.fitted$pi[1]
pi0hat.ashr[i] <- ashr::get_pi0(ashr.fit)
FDPq.biashr[i] <- mean(theta[biashr.fit$theta.lfdr <= q] == 0)
FDPq.ashr[i] <- mean(theta[ashr::get_lfdr(ashr.fit) <= q] == 0)
}

pdf("~/GitHub/Bi-ASH/output/biashr.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
boxplot(pi0hat.ashr, pi0hat.biashr, names = c("ashr", "biashr"), ylab = expression(hat(pi[0])))
abline(h = prop.null, col = 2, lty = 2, lwd = 2)
boxplot(FDPq.ashr, FDPq.biashr, names = c("ashr", "biashr"), ylab = "FDP")
abline(h = q, col = 2, lty = 2, lwd = 2)
dev.off()