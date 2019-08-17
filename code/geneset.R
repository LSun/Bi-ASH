source('~/GitHub/Bi-ASH/code/biashr.R')

p53 <- readRDS("~/GitHub/Bi-ASH/data/p53.rds")
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z

p53.dat <- p53[-which(abs(X / s) >= 4), ]

gene.num <- length(X)
gene.names <- rownames(p53.dat)
gene.set.all <- msigdb::read.gmt("~/GitHub/Bi-ASH/data/c2.symbols.gmt")
gene.set.num <- length(gene.set.all$genesets)
gene.set.size.original <- sapply(gene.set.all$genesets, length)
names(gene.set.size.original) <- NULL
gene.set.name <- gene.set.all$geneset.names

gene.set.size <- pi0 <- loglikratio <- converged <- c()

for (i in seq(gene.set.num)) {
  gene.set <- gene.set.all$genesets[[i]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  x1 <- X[gene.set.index]
  s1 <- s[gene.set.index]
  x2 <- X[setdiff(seq(gene.num), gene.set.index)]
  s2 <- s[setdiff(seq(gene.num), gene.set.index)]
  biashr.fit <- biashr(x1, s1, x2, s2)
  gene.set.size[i] <- length(gene.set.index)
  pi0[i] <- biashr.fit$g.fitted$pi[1]
  loglikratio[i] <- biashr.fit$loglikratio
  converged[i] <- biashr.fit$converged
}

gene.set.order <- order(loglikratio, decreasing = TRUE)
res <- cbind.data.frame(
  geneset = gene.set.name[gene.set.order],
  size = gene.set.size.original[gene.set.order],
  size.used = gene.set.size[gene.set.order],
  loglikratio = loglikratio[gene.set.order],
  pi0 = pi0[gene.set.order],
  converged = converged[gene.set.order]
)

saveRDS(res, "~/GitHub/Bi-ASH/output/p53.res.rds")