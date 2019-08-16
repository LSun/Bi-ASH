p53 <- readRDS("~/GitHub/Bi-ASH/data/p53.rds")
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z

gene.names <- rownames(p53)
