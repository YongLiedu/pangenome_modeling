#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Script to model pangenome expansion and core genome reduction
#' Authors: Yong LI
#' Date: 2021/06/09
#' Institution: Northeast Agricultural University
#' Usage: Rscript pangenome_modelling.R pav.txt
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library("minpack.lm")
library("parallel")

args = commandArgs(T)
pav <- read.table("melon_MR1_pan_g0_c02.pav", sep = "\t", header = T, row.names = 1)

sample_cnt <- 10000
threads <- 120

cal_gene_cnt <- function(x) {
    pav_samples <- as.matrix(pav[, sort(sample(1:varieties, x))])
    pav_samples_rs <- rowSums(pav_samples)
    core <- length(which(pav_samples_rs == x))
    pan <- length(which(pav_samples_rs > 0))
    return(c(x, core, pan))
}

varieties <- ncol(pav)
gene_cnt <- data.frame(genomes = integer(), core = integer(), pan = integer())
for (i in 1:varieties) {
    cl <- makeCluster(threads, type = "FORK")
    t_gene_cnt <- parLapply(cl, rep(i, sample_cnt), cal_gene_cnt)
    comb_gene_cnt <- do.call("rbind", t_gene_cnt)
    gene_cnt <- rbind(gene_cnt, comb_gene_cnt)
    stopCluster(cl)
}
colnames(gene_cnt) <- c("genomes", "core", "pan")
save.image("pangenome_modelling.RData")  # once an nlsLM error happens below, load this file to avoid re-running above codes.

#' adjust the values of A, B and C to suit your data
#' unsuitable values could cause in nlsLM error
para0p <- c(A = 9000, B = -5, C = 0.1)
para0c <- c(A = 9000, B = -5, C = 0.1)
fitp <- nlsLM(pan ~ A * genomes^B + C, gene_cnt, start = para0p, trace = T)
fitc <- nlsLM(core ~ A * exp(B * genomes) + C, gene_cnt, start = para0c, trace = T)
print(fitp)
print(fitc)
confint(fitp, parm = c("A", "B", "C"), level = 0.95)
confint(fitc, parm = c("A", "B", "C"), level = 0.95)
print(summary(fitp))
print(summary(fitc))
new <- data.frame(xdata = seq(min(gene_cnt$genomes), max(gene_cnt$genomes), len = 2000))

## plot - plot
r <- range(gene_cnt$genomes)
xNew <- seq(r[1], r[2], length.out = 2000)
yNewp <- predict(fitp, list(genomes = xNew))
yNewc <- predict(fitc, list(genomes = xNew))
png("pangenome_modelling.png", width = 2000, height = 1500, res = 300)
ymin <- min(min(gene_cnt$core), min(yNewc))
ymax <- max(max(gene_cnt$pan), max(yNewp))
plot(gene_cnt$genomes, gene_cnt$pan, xlab = "Number of genomes", ylab = "Number of genes", 
    pch = 15, col = "#A6CEE3", ylim = c(ymin, ymax), xaxt = "n", yaxt = "n")
axis(1, at = seq(0, varieties, by = ceiling(varieties/30) * 5))
yticks <- seq(0, ymax, by = 1000)
labels <- format(yticks, big.mark = ",", scientific = FALSE)
axis(2, at = yticks, labels = labels)
points(gene_cnt$genomes, gene_cnt$core, pch = 17, col = "#1F78B4")
lines(xNew, yNewp, lwd = 2)
lines(xNew, yNewc, lwd = 2)
legend("bottomleft", legend = c("pan-genome", "core-genome"), pch = c(15, 17), col = c("#A6CEE3", 
    "#1F78B4"))
dev.off()
save.image("pangenome_modelling.RData")

## plot - ggplot2
library(ggplot2)
library(reshape2)
gene_cnt_gg <- melt(gene_cnt, id.vars=c("genomes"), measure.vars=c("core","pangenome"), variable.name="genome_type", value.name="genes")
p_panmodel <- ggplot(gene_cnt_gg, aes(x=genomes, y=genes, colour=genome_type)) + geom_point() + geom_smooth(method="loess")
ggsave("pangenome_modelling_gg.png", plot=p_panmodel, width = 6.5, height = 5, dpi = 300)
#ggsave("pangenome_modelling_gg.svg", plot=p_panmodel, width = 6.5, height = 5)
save.image("pangenome_modelling.RData")
