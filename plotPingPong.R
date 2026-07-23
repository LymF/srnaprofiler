library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], header = TRUE, sep = "\t")
d$fill <- ifelse(d$overlap == 10, "ping-pong (10)", "other")
pdf(paste0(args[2], ".pingpong.pdf"), width = 6, height = 4)
ggplot(d, aes(x = overlap, y = pairs, fill = fill)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ping-pong (10)" = "#e41a1c", "other" = "grey60")) +
  theme_classic(base_size = 16) +
  xlab("5'-5' overlap (nt)") + ylab("Weighted pairs") +
  theme(legend.title = element_blank())
dev.off()
