args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], header = TRUE, sep = "\t")
x <- d$bin / max(d$bin)
pdf(paste0(args[2], ".coverage.pdf"), width = 6, height = 4)
plot(x, d$plus, type = "l", col = "#377eb8", lwd = 2,
     ylim = range(c(d$plus, d$minus), 0), xlab = "Normalized position",
     ylab = "Mean weighted coverage", main = "Coverage profile")
lines(x, d$minus, col = "#e41a1c", lwd = 2)
legend("topright", legend = c("+", "-"), col = c("#377eb8", "#e41a1c"), lwd = 2, bty = "n")
dev.off()
