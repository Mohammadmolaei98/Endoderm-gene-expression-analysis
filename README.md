library(pheatmap)
library(ggplot2)
library(reshape)
library(grid)
library(rgl)

theme_complete_bw <- function(base_size = 24, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
      axis.ticks =        element_line(colour = "black"),
      axis.title.x =      element_text(size = base_size, vjust = 0.5),
      axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.margin = unit(0.1, "cm"),
      
      legend.background = element_rect(colour=NA), 
      legend.key =        element_rect(fill =NA, colour = "black", size = 0.25),
      legend.key.size =   unit(1.5, "lines"),
      legend.text =       element_text(size = base_size * 0.7),
      legend.title =      element_text(size = base_size * 0.8),
      legend.position =   "top",
      
      panel.background = element_rect(fill = "white", colour = NA), 
      panel.border =     element_rect(fill = NA, colour = "black", size=2), 
      panel.grid.major = element_line(colour = NA, size = 0.2), #"grey"
      panel.grid.minor = element_line(colour = NA, size = 0.5), #"grey"
      panel.margin =     unit(0.25, "lines"),
      
      strip.background = element_rect(fill = NA, colour = NA), 
      strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
      strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = +90),
      
      plot.background =  element_rect(colour = NA, fill = "white"),
      plot.title =       element_text(size = base_size*.8),
      plot.margin =      unit(c(1, 1, .5, .5), "lines"))
}

setwd("~/Documents/Rproject/pro1/")
x = read.delim("~/Desktop/R Class/Endoderm.txt")
heatmap1 <- pheatmap(x[,-1])
pdf("heatmap1.pdf")
heatmap1
dev.off()


#select rowmanes
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
head(x)
pheatmap2 <- pheatmap(x,fontsize_row=5,border_color=NA)#fontsize
pdf("heatmap2.pdf")
pheatmap2
dev.off()

# normalazie with log2
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)
logpheatmap <- pheatmap(x,fontsize_row=5,border_color=NA)#fontsize
pdf("logpheatmap.pdf")
logpheatmap
dev.off()


#new session
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
#t.test
a <- x[1:3,1]
a
b <- x[4:6,1]
b
t.test(a,b)
t.test(x[1:3,"PAX4"], x[4:6, "PAX4"])
#new session
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
head(x)
x = t(x)# tame data ha bayed ya adad ya harf bashe # dastore t data ro matrix mikonehhh
head(x)
class(x)
x <- data.frame(x)
ggplot1 <- ggplot(x, aes(x=DE.1, y= P5PE.3)) + geom_point()
pdf("ggplot1")
ggplot1
dev.off()
#set label
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
head(x)
x = t(x)

x$Gene = rownames(x)
head(x)
googooli=x
ggplot2 = ggplot(googooli, aes(x=DE.1, y=DE.2, label=Gene)) + geom_point()
ggplot2 = ggplot2 + geom_text()
pdf("ggplot2")
ggplot2
dev.off()
#barplot
#y <- x[,"DE.1"] #tamame esme gene ha hazf shood ________ pass mishe  y <- x[,"DE.1", drop=FALSE]
#ya 

y = x[,c("Gene","DE.1")]
p1 = ggplot(y, aes(x=Gene,y=DE.1,fill=Gene))+geom_bar(stat="identity") #stat="identity >> hamini ke midam ro bekesh 
p2 = p1 + ylab("Expression of genes in Defenitive Endoderm (log2)")

pdf("barplot1")
p1
p2
dev.off()

p3 = ggplot(y, aes(x=Gene,y=DE.1,fill=DE.1))+geom_bar(stat="identity") #stat="identity >> hamini ke midam ro bekesh 
pdf("barplot2")
p3
dev.off()
# session5


x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x.cor <- cor(x)
head(x.cor) #NA data
any(is.na(x)) # T means there are NA data
x <- na.omit(x) # remove NA data
x.cor <- cor(x)
head(x.cor)
pheatmap(x.cor)
pheatmat3 = pheatmap(x.cor)
pdf("pheatmat3.pdf")
pheatmat3
dev.off()

ggplot3 = ggplot(x,aes(NEUROD1, NKX6.1))+geom_point()
ggplot4 = ggplot(x,aes(HLXb9, HHEX))+geom_point()
pdf("ggplot3,4.pdf")
ggplot3
ggplot4
dev.off()


ggplot5 = ggplot(x,aes(NEUROD1, NKX6.1))+geom_point()+geom_smooth(method = "lm")
ggplot6 = ggplot(x,aes(HLXb9, HHEX))+geom_point()+geom_smooth(method = "lm")
pdf("ggplot5,6.pdf")
ggplot5
ggplot6
dev.off()

#sample  correlation
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)

y <- t(x)
y.cor <- cor(y)#error
y.cor <- na.omit(y.cor)
dim(y.cor)
head(y)
y[,2:3]=y[,2:3]+runif(18,min=-0.001,max=+0.001)
dim(y[,2:3])
y.cor <- cor(y)
pheatmap(y.cor)
pheatmap4 = pheatmap(y.cor)
pdf("pheatmap4.pdf")
pheatmap4
dev.off()
####session6####
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
summary(x)
summary(x[,1])
#library(reshape)
x.m <- melt(x)
head(x.m)
ggplot7 <- ggplot(x.m, aes(x= variable,y=value))+geom_boxplot()
ggplot8 <- ggplot(x.m, aes(x= variable,y=value,fill = variable))+geom_boxplot(outlier.size =0)
ggplot9 <- ggplot(x.m, aes(x= variable,y=value,fill = variable))+geom_violin(outlier.size =0)

pdf("ggplot7,8,9.pdf")
ggplot7
ggplot8
ggplot9
dev.off()

####session7####
setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
pheatmap5 <- pheatmap(x)

pheatmap6 <- pheatmap(x, clustering_distance_rows = "correlation"
         , clustering_distance_cols ="correlation") 
#vojod e maghdire yeksan dar SC2 and SC3 ke dar satr e 2 3 ast 
x <- na.omit(x)
head(x)
x <- x[-2:-3 , ]

head(x)

pheatmap6 <- pheatmap(x, clustering_distance_rows = "correlation"
                      , clustering_distance_cols ="correlation") 


pdf("pheatmap5.pdf")
pheatmap5

dev.off()
pdf("pheatmap6.pdf")

pheatmap6
dev.off()

#cor
cor(x[,c("HLXb9", "HHEX")])
cor(x[,c("ISL1", "HHEX")])
ggplot3 = ggplot(x,aes(NEUROD1, NKX6.1))+geom_point()
ggplot4 = ggplot(x,aes(HLXb9, HHEX))+geom_point()
pdf("ggplot3,4.pdf")
ggplot3
ggplot4
dev.off()


ggplot5 = ggplot(x,aes(NEUROD1, NKX6.1))+geom_point()+geom_smooth(method = "lm")
ggplot6 = ggplot(x,aes(HLXb9, HHEX))+geom_point()+geom_smooth(method = "lm")
pdf("ggplot5,6.pdf")
ggplot5
ggplot6
dev.off()



setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)

apply(x,1,min) #min data 1 har satr
#na.rm=T >>>ignore data
#na.omit >>> hazf e kamel

lapply(1:3, function(x) x^2)    #1 ta 3 bordar ya vector ast
sapply(1:3, function(x) x^2) 

t.test(x[4:10,],x[11:32,])

MyTest <- function(i) {
  t.test(x[4:10,i], x[11:32,i])$p.value
}

#MyTest(1) baray e soton e yek

sapply(1:ncol(x),MyTest)



MyTest2 <- function(y) {
  t.test(y[4:10], y[11:32])$p.value
}

apply(x,2,MyTest2)
sapply(1:9, MyTest)
MyTest3 <- function(googool){
  t.test(x[googool,1:5], x[googool, 6:10])$p.value
}


####session8####
setwd("~/Documents/Rproject/pro1/")
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
head(x)

se <- function(x) sd(x)/sqrt(length(x))
#sd(c(3,5,7))
#se(c(3,5,7))
se(x[4:6,1])
sd(x[4:6,1])

####session9####
#principal component analysis
setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)
pc <- prcomp(x) # error not available
pc <- prcomp(x)
head(pc)
plot(pc)

pcx <- data.frame(pc$x) 
head(pcx)

pcx$Sample <- rownames(pcx)
pcx$Sample <- substr(pcx$Sample, 1, nchar(pcx$Sample)-2) #az harf e#
ggplot10 <- ggplot(pcx, aes(PC1,PC2,color=Sample))+geom_point(size=3)+theme_complete_bw()#har nighte sample ast#
ggplot10

pdf("ggplot10.pdf")
ggplot10
dev.off() 

#rotation zarib e ahamiat gene ha dar tafavot e nighat ya variation #
# chand mesal#
barplot(x$PDX1)
barplot(x$NEUROD1)
barplot(x$PAX4)
 
pcr <- data.frame(pc$rotation)
pcr$Gene <- rownames(pcr)
ggplot11 <- ggplot(pcr, aes(PC1,PC2,label=Gene))+geom_text()+theme_complete_bw() #har nighte gene ast#
ggplot11

pdf("ggplot11.pdf")
ggplot11
dev.off() 

pdf("ggplot10,11.pdf")
ggplot10
ggplot11
dev.off() 

colnames(x)
x[,3,drop=F] #data frame bemoneh#
plot3d(x[,1:3])

####session10####
#anova#
setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)
x$Sample <- substr(rownames(x), 1, nchar(rownames(x))-2)
head(x)

y <- x[,c("HLXb9","Sample")]
y


AOV <- function(Gene) {
  y <- x[,c(Gene,"Sample")]
  colnames(y)[1]="Gene"
  anova(aov(Gene ~ Sample, y))[1,5]
}

AOV("HLXb9")
genes <- colnames(x)
head(genes)
genes <- genes[genes!="Sample"]
head(genes)
sapply(genes, AOV)

ggplot12 <- ggplot(x, aes(Sample,HLXb9))+geom_boxplot()
ggplot12
pdf("ggplot12.pdf",width=14,height=14)
ggplot12
dev.off()
####session11####
setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)
#correalation#

head(x.t)
x.t <- t(x)
x.t <- x.t[,-2:-3] #do sate aval NA ast
xc <- cor(x.t)
pheatmap7 <- pheatmap(xc)
pdf("pheatmap7.pdf",width=14,height=14)
pheatmap7
dev.off()

x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)
#correalation#

head(x.t)
x.t <- t(x)
x.t <- x.t[,-2:-3] #do sate aval NA ast
xc <- cor(x.t,method = "spearman")
pheatmap8 <- pheatmap(xc)
pdf("pheatmap8.pdf",width=14,height=14)
pheatmap8
dev.off()



pc <- prcomp(x)
pcx <- data.frame(pc$x)
pcx$Sample <- rownames(pcx)
pcx$Sample <- substr(pcx$Sample, 1, nchar(pcx$Sample)-2)
pdf("ggplot13.pdf",width=10,height=10)
ggplot(pcx, aes(PC1, PC2, color=Sample))+geom_point(size=3)+theme_complete_bw()
dev.off()

#anova
setwd("~/Documents/Rproject/pro1/")
x = read.delim("Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)

x.m <- melt(as.matrix(x))
head(x.m)
colnames(x.m) <- c("Sample","Gene","Exp")
anova(aov(Exp~Gene+Sample,x.m))


