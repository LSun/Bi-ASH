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

p53 <- p53[-which(abs(X / s) >= 4), ]
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z

gene.num <- length(z)
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

saveRDS(list(res, gene.set.order), "~/GitHub/Bi-ASH/output/p53.res.rds")

res.list <- readRDS("~/GitHub/Bi-ASH/output/p53.res.rds")
res <- res.list[[1]]
gene.set.order <- res.list[[2]]

p53 <- readRDS("~/GitHub/Bi-ASH/data/p53.rds")
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z
which(abs(X / s) >= 4)

p53 <- p53[-which(abs(X / s) >= 4), ]
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z

z.model <- z

p53 <- readRDS("~/GitHub/Bi-ASH/data/p53.rds")
design <- model.matrix(~colnames(p53))
lim <- limma::lmFit(p53, design)
r.ebayes <- limma::eBayes(lim)
p <- r.ebayes$p.value[, 2]
t <- r.ebayes$t[, 2]
z <- -sign(t) * qnorm(p/2)
X <- lim$coefficients[, 2]
s <- X / z

z.plot <- z
z.plot[abs(z) < 4] <- z.model

gene.num <- length(z.plot)
gene.names <- names(z.plot)

pdf("~/GitHub/Bi-ASH/output/genesetintro.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
pcol <- rep(1, length(z.plot))
pcol[abs(z.plot) >= 5] <- "red"
pch <- rep(1, length(z.plot))
pch[abs(z.plot) >= 5] <- 19
qqnorm(z.plot, main = "(a) Normal Q-Q Plot for all z-scores",
       xlim = c(-4, 5.2),
       ylim = c(-4, 5.2),
       pch = pch,
       col = pcol)
abline(0, 1, lty = 2, lwd = 2, col = "blue")
legend("topleft", col = "blue", lty = 2, lwd = 2, "The y = x line")
text(x = -qnorm(mean(z.plot[which(z.plot >= 5)][1] <= z.plot)),
     y = z.plot[which(z.plot >= 5)][1],
     labels = names(z.plot[which(z.plot >= 5)][1]),
     pos = 4,
     col = "red",
     cex = 0.75)
text(x = -qnorm(mean(z.plot[which(z.plot >= 5)][2] <= z.plot)),
     y = z.plot[which(z.plot >= 5)][2],
     labels = names(z.plot[which(z.plot >= 5)][2]),
     pos = 2,
     col = "red",
     cex = 0.75)
# hist(z.plot, prob = TRUE, xlab = "z-score", main = expression("Histogram of z-scores of all genes"))
# lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = "blue", lty= 2, lwd = 2)
# legend("topright", "N(0,1)", lty = 2, col = "blue", lwd = 2)
hist(res$size.used, breaks = 50, xlab = "Number of genes in a gene set",
     # main = expression(paste("Histogram of gene set sizes ", p[1]))
     main = "(b) Histogram of gene set sizes"
     )
dev.off()

topnumber <- 10

pdf("~/GitHub/Bi-ASH/output/genesetranking.pdf", width = 5, height = 5)
par(mar=c(5.1, 4.1, 2.1, 2.1))
plot(res$avgloglikratio, ylab = "enrichment score", cex = 1.5,
     pch = c(rep(10, topnumber), rep(1, 522 - topnumber)),
     col = c(rep("red", topnumber), rep(1, 522 - topnumber)))
dev.off()

pdf("~/GitHub/Bi-ASH/output/genesetz.pdf", width = 10, height = 12)

par(mar=c(5.1, 4.1, 1.1, 2.1))

hist(z.plot, prob = TRUE,
     ylim = c(-0.20 * (topnumber + 2.5), 0.26),
     main = "",
     yaxt = 'n',
     xaxt = "n",
     ylab = "",
     xlab = "z-score",
     xlim = c(-5.5, 5.5))
axis(side = 1, at = seq(-6, 6, by = 1))
lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = "blue", lty= 2, lwd = 2)
legend("topright", "N(0,1)", lty = 2, col = "blue", lwd = 2)

gene.set.chosen.order <- gene.set.order[1 : topnumber]
gene.set.chosen.common <- gene.set.order[c(1, 2, 4, 5, 6, 10)]
gene.set.chosen.have <- gene.set.order[c(3, 7, 8, 9)]
gene.set.chosen.havenot <- gene.set.order[c(326, 369)]
geneset.col <- scales::hue_pal()(3)

for (i in seq(gene.set.chosen.common)) {
  gene.set <- gene.set.all$genesets[[gene.set.chosen.common[i]]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  z1 <- z.plot[gene.set.index]
  points(z1, rep(-0.2 * i, length(gene.set.index)), col = geneset.col[1], pch = i + 1)
  text(-6, -0.2 * i, pos = 4, gene.set.name[gene.set.chosen.common[i]], cex = 1, col = geneset.col[1])
}

for (i in seq(gene.set.chosen.have)) {
  gene.set <- gene.set.all$genesets[[gene.set.chosen.have[i]]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  z1 <- z.plot[gene.set.index]
  points(z1, rep(-0.2 * (i + length(gene.set.chosen.common)), length(gene.set.index)), col = geneset.col[2], pch = i + 1 + length(gene.set.chosen.common))
  text(-6, -0.2 * (i + length(gene.set.chosen.common)), pos = 4, gene.set.name[gene.set.chosen.have[i]], cex = 1, col = geneset.col[2])
}

for (i in seq(gene.set.chosen.havenot)) {
  gene.set <- gene.set.all$genesets[[gene.set.chosen.havenot[i]]]
  gene.set.index <- purrr::discard(match(gene.set, gene.names), is.na)
  z1 <- z.plot[gene.set.index]
  points(z1, rep(-0.2 * (i + length(gene.set.chosen.common) + length(gene.set.chosen.have)), length(gene.set.index)), col = geneset.col[3], pch = i + 1 + length(gene.set.chosen.common + length(gene.set.chosen.have)))
  text(-6, -0.2 * (i + length(gene.set.chosen.common) + length(gene.set.chosen.have)), pos = 4, gene.set.name[gene.set.chosen.havenot[i]], cex = 1, col = geneset.col[3])
}

dev.off()