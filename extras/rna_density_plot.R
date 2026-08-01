#Single-end
library(ggplot2)

data <- read.csv("pv1negative.tabular",sep="\t",header=F)
d2 <-  data[,c(2,3)]
m<-max(d2$V3)
d2 <-d2[order(d2[,1]),]
ggplot(d2,aes(x=V2,y=V3),size=0.5)+ 
  scale_y_continuous(limits = c(0, m)) +ylab("RNA density") +
  xlab("PV1") + ylim(0,as.integer(m)) +
  theme_classic(base_size=20)+ 
  geom_density(aes(y=V3,col="V3"),adjust=5,fill="#0099CC",colour="#0099CC",size=0.1,stat="identity") + 
  theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) + 
  ggtitle("SRR27896439") + theme(plot.title = element_text(hjust = 0.5)) 

max(d2$V3)


#Paired end
data <- read.csv("pv1both.tabular",sep="\t",header=F)
d2 <-  data[,c(2,3,4)]
m <- max(d2[,3],d2[,2]) + 10
d2 <-d2[order(d2[,1]),]
d2[,3] <- d2[,3] * -1
ggplot(d2,aes(x=V2,y=V3),size=0.5)+ 
  scale_y_continuous(limits = c(-m, m)) +ylab("RNA density") +
  xlab("PV1") + ylim(-as.integer(m),as.integer(m)) +
  theme_classic(base_size=20)+ 
  geom_density(aes(y=V3,col="V3"),adjust=5,fill="#0099CC",colour="#0099CC",size=0.1,stat="identity") + 
  geom_density(aes(y=V4,col="V4"),colour="#996631",fill="#996631",size=0.1,stat="identity")  + 
  theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) + 
  ggtitle("SRR27896439") + theme(plot.title = element_text(hjust = 0.5)) 
