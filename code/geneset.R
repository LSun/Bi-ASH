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

# p53 <- p53[-which(abs(X / s) >= 4), ]

gene.num <- length(X)
gene.names <- rownames(p53)
gene.set.all <- msigdb::read.gmt("~/GitHub/Bi-ASH/data/c2.symbols.gmt")
gene.set.num <- length(gene.set.all$genesets)
gene.set.size.original <- sapply(gene.set.all$genesets, length)
names(gene.set.size.original) <- NULL
gene.set.name <- gene.set.all$geneset.names

gene.set.size <- pi0 <- loglikratio <- avgloglikratio <- converged <- c()

for (i in seq(gene.set.num)) {
  gene.set <- gene.set.all$genesets[[i]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  x1 <- X[gene.set.index]
  s1 <- s[gene.set.index]
  x2 <- X[setdiff(seq(gene.num), gene.set.index)]
  s2 <- s[setdiff(seq(gene.num), gene.set.index)]
  biashr.fit <- biashr(x1, s1, x2, s2, omega.null.weight = 0)
  gene.set.size[i] <- length(gene.set.index)
  pi0[i] <- biashr.fit$g.fitted$pi[1]
  loglikratio[i] <- biashr.fit$loglikratio
  avgloglikratio[i] <- biashr.fit$loglikratio / length(gene.set.index)
  converged[i] <- biashr.fit$converged
}

gene.set.order <- order(avgloglikratio, decreasing = TRUE)
res <- cbind.data.frame(
  geneset = gene.set.name[gene.set.order],
  size = gene.set.size.original[gene.set.order],
  size.used = gene.set.size[gene.set.order],
  loglikratio = loglikratio[gene.set.order],
  avgloglikratio = avgloglikratio[gene.set.order],
  pi0 = pi0[gene.set.order],
  converged = converged[gene.set.order]
)

saveRDS(res, "~/GitHub/Bi-ASH/output/p53.res.rds")

res <- readRDS("~/GitHub/Bi-ASH/output/p53.res.rds")

pdf("~/GitHub/Bi-ASH/output/genesetsize.pdf", width = 6, height = 4)
hist(gene.set.size, breaks = 50, xlab = "number of genes in a gene set", main = "Histogram of gene set sizes")
dev.off()

topnumber <- 6
res.diff <- c(4, 15, 374, 332)
res.diff.original <- gene.set.order[res.diff]
geneset.col <- scales::hue_pal()(3)
geneset.col <- rep(geneset.col, c(1, 1, 2))
all.geneset.col <- rep(1, gene.set.num)
all.geneset.col[1 : topnumber] <- "red"
all.geneset.col[res.diff] <- geneset.col
all.geneset.pch <- rep(1, gene.set.num)
all.geneset.pch[1 : topnumber] <- 10
all.geneset.pch[res.diff[-1]] <- 19
pdf("~/GitHub/Bi-ASH/output/genesetranking.pdf", width = 5, height = 5)
par(mar=c(5.1, 4.1, 2.1, 2.1))
plot(res$avgloglikratio, ylab = "enrichment score", cex = 1.5,
     pch = all.geneset.pch,
     col = all.geneset.col)
text(setdiff(seq(topnumber), res.diff[1]), res$avgloglikratio[setdiff(seq(topnumber), res.diff[1])], res$geneset[setdiff(seq(topnumber), res.diff[1])], pos = 4, cex = 0.75, col = "red")
text(res.diff[1], res$avgloglikratio[res.diff[1]], res$geneset[res.diff[1]], adj = c(-0.075, 1), cex = 0.75, col = geneset.col[1])
text(res.diff[2], res$avgloglikratio[res.diff[2]], res$geneset[res.diff[2]], pos = 4, cex = 0.75, col = geneset.col[2])
text(res.diff[3], res$avgloglikratio[res.diff[3]], res$geneset[res.diff[3]], adj = c(0.25, -1.2), cex = 0.75, col = geneset.col[3])
text(res.diff[4], res$avgloglikratio[res.diff[4]], res$geneset[res.diff[4]], adj = c(0.85, -1.2), cex = 0.75, col = geneset.col[4])
dev.off()

pdf("~/GitHub/Bi-ASH/output/genesetz.pdf", width = 9, height = 6)
hist(z, prob = TRUE,
     ylim = c(-0.02 * length(res.diff), dnorm(0)),
     xlab = "z-score", main = "Histogram of all z-scores")
legend("topleft", col = "blue", lty = 2, "N(0,1)", bty = "n", cex = 1)
lines(seq(-6, 6, by = 0.1), dnorm(seq(-6, 6, by = 0.1)), lty = 2, col = "blue")
for (i in seq(res.diff.original)) {
  gene.set <- gene.set.all$genesets[[res.diff.original[i]]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  x1 <- X[gene.set.index]
  s1 <- s[gene.set.index]
  z1 <- x1 / s1
  points(z1, rep(-0.02 * i, length(gene.set.index)), col = geneset.col[i], pch = i + 1)
  text(-4.3, -0.02 * i, pos = 4, gene.set.name[res.diff.original[i]], cex = 0.75, col = geneset.col[i])
}
dev.off()
