args <- commandArgs(trailingOnly = TRUE)
plus_file <- args[1]; minus_file <- args[2]
period <- as.numeric(args[3]); out <- args[4]
panel <- function(file, title, period) {
  if (file.exists(file) && file.info(file)$size > 0) {
    d <- read.table(file, header = TRUE, sep = "\t")
    d <- d[order(d$distance), ]
    plot(d$distance, d$count, type = "l", lwd = 1, col = "black",
         main = title, xlab = "5'-5' distance (nt)", ylab = "Weighted pairs",
         panel.first = grid(lty = "dotted", col = "lightgray"))
    for (m in seq(period, max(d$distance, na.rm = TRUE), by = period))
      abline(v = m, lty = 2, col = "red")
  } else {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(1, 1, "No data")
  }
}
pdf(paste0(out, ".phasing.pdf"), width = 5, height = 6)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 2))
panel(plus_file, "Sense (+)", period)
panel(minus_file, "Antisense (-)", period)
dev.off()
