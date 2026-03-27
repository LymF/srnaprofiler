library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
# args[1] = TSV; args[2] = outprefix; args[3] = y_label

data <- read.csv(args[1], sep="\t", header=TRUE)
data[is.na(data)] <- 0

m <- max(max(rowSums(data)), -min(rowSums(data))) * 1.1
b <- melt(data, id.vars="size")
b1 <- subset(b, value >= 0)
b2 <- subset(b, value < 0)

pdf(paste(args[2], ".pdf", sep=""), width=7, height=4)
ggplot() + theme_classic(base_size=20) +
  geom_bar(data=b1, aes(size, value, fill=variable), stat="identity", width=.8) +
  geom_bar(data=b2, aes(size, value, fill=variable), stat="identity", width=.8) +
  ylim(-m, m) + xlab("small RNA length (nt)") +
  ylab(paste("Number of small RNAs \n(", args[3], ")", sep="")) +
  theme(panel.border=element_blank(), axis.line=element_line(color='black')) +
  scale_fill_manual(values=c("A"="#108041","C"="#3B54A2","G"="#F9A31B","U"="#EA2627"),
                    name="5' base \ndistribution") +
  scale_x_continuous(breaks=seq(15,35,1), labels=c("","16","","18","","20","","22","","24","","26","","28","","30","","32","","34",""))
dev.off()
