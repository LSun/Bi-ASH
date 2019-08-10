r <- readRDS("~/GitHub/Bi-ASH/data/liver.rds")
source('~/GitHub/cashr_paper/code/RNAseq_pipeline.R')
source('~/GitHub/Bi-ASH/code/biashr.R')

set.seed(777)

nctrl <- 1e3
nnctrl <- 1e3
prop.null <- 0.9
nsamp <- 20
q <- 0.1

Y = lcpm(r)
subset = top_genes_index(nctrl + nnctrl, Y)
r = r[subset, sample(ncol(r), nsamp)]
pi0hat <- FDPq <- c()

for (i in 1 : 600) {
r.signaladded <- seqgendiff::thin_2group(mat = as.matrix(r), 
                     prop_null = (prop.null * nnctrl + nnctrl) / (nctrl + nnctrl), 
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 0.75))

gene.ctrl <- sample(seq(nctrl + nnctrl)[r.signaladded$coefmat == 0], nctrl)
gene.nonctrl <- setdiff(seq(nctrl + nnctrl), gene.ctrl)
theta <- r.signaladded$coefmat[gene.nonctrl]

summary.stat <- count_to_summary(r.signaladded$mat, r.signaladded$designmat)

x1 <- summary.stat$X[gene.nonctrl]
s1 <- summary.stat$s[gene.nonctrl]
x2 <- summary.stat$X[gene.ctrl]
s2 <- summary.stat$s[gene.ctrl]

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