## import CNV varscan data
df<-read.table("varscan.bed",h=T,sep="\t")


## reorganized data (optional)
df$id<-paste0(df$chrom,df$chr_start)
x <- "id"
df<-df[c(x, setdiff(names(df), x))]
df2<-df[,-9]
write.table(df2,file="B9CN.tsv",row.names=F)


## load libraries
library(ggplot2)
library(dplyr)
library(showtext)
library(pacman)
library(extrafont)
library(zoo)


## custom text font
font.families.google()

p_load(showtext)
font_add_google("Gochi Hand", "gochi")
font_add_google("Schoolbell", "bell")
font_add_google("Bad Script")


## chromosome 1

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr1")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 1(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")




## chromosome 2

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr2")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 2(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 3

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr3")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 3(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 4

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr4")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 4(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 5

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr5")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 5(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 6

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr6")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 6(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")

## chromosome 7

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr7")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 7(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 8

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr8")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 8(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 9

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr9")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 9(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 10

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr10")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 10(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 11

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr11")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 11(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")

## chromosome 12

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr12")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 12(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 13

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr13")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 13(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")




## chromosome 14

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr14")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 14(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")

## chromosome 15

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr15")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 15(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 16

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr16")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 16(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 17

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr17")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 17(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")


## chromosome 18

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr18")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 18(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome 19

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chr19")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome 19(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")




## chromosome x

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chrX")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome X(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")



## chromosome M

dev.new(width = 600, height = 300, unit = "px")
df%>%filter(chrom=="chrM")%>%
ggplot(aes(x=chr_start,y=log2_ratio))+
ylim(-5,3)+
labs(x=" MM10 - Chromosome M(bp)",
		y="CNV log2 Ratio")+
	theme(axis.line.x = element_line(size = .5, colour = "black"),
	axis.line.y = element_line(size = .5, colour = "black"),
	axis.text.x = element_text(colour = "black", size = 16),
	axis.text.y = element_text(colour = "black", size = 16),
	legend.key = element_rect(fill = "white", colour = "white"),
	legend.position = "bottom", legend.direction = "horizontal",
	legend.title = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), panel.border = element_blank(),
	panel.background = element_blank(),
	plot.title = element_text(family = "bell",size=16),
	text = element_text(family = "bell",size=16))+
geom_point(colour = "red",size=1,alpha = 0.2)+
geom_line((aes(y=rollmean(log2_ratio, 10, na.pad=TRUE))),colour = "steelblue")+
geom_hline(yintercept=0.5, linetype="dashed", color = "darkgreen")
