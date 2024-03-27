library(circlize)
#chrlen
chrlen<-read.csv("input/chr_length.csv",header=F)

#genedensity
library(dplyr)
library(tidyr)
gene <- read.table("input/evm_gene.gff",sep="\t",header=F)
colnames(gene) <- c("chr","source","type","start","end",
                   "score","strand","phase","attributes")
gene <- gene[!substr(gene$chr,1,3)=='CTG',]
gene <- gene[,c(1,4)]
gene <- arrange(gene,chr,start)
#gene$window <- gene$start%/%90000
#genenub <- gene %>%
#   group_by(chr,window) %>%
#   summarise(nub=length(window))
#genenub$start <- genenub$window*90000
#genenub$end <- genenub$start+90000
#genenub <- genenub[,c(1,4,5,3)]

#value_max <- max(genenub$nub)
#colorsChoice <- colorRampPalette(c('#FA9082', '#FF2933'))
#gene_color_assign <- colorRamp2(breaks = c(0:value_max), col = colorsChoice(value_max + 1))
#"#A05050FF"

#snpdensity 
snpnub<-read.csv("input/snp_chr_pos.csv",header=F)
snpnub<-separate(data = snpnub ,col = V1,into = c("chr","pos"),sep = "[:]+")
snpnub <- snpnub[!substr(snpnub$chr,1,3)=='CTG',]
snpnub <- arrange(snpnub,chr,pos)
#snpnub$window <- as.numeric(snpnub$pos)%/%9000
#snpnub <- snpnub %>%
#   group_by(chr,window) %>%
#   summarise(nub=length(window))
#snpnub$start <- snpnub$window*9000
#snpnub$end <- snpnub$start+9000
#snpnub <- snpnub[,c(1,4,5,3)]
#colnames(snpnub)<-c("chr","start","end","nub")


#value_max <- max(snpnub$nub)
#colorsChoice <- colorRampPalette(c('#FFFFFF', 'blue'))
#snp_color_assign <- colorRamp2(breaks = c(0:value_max), col = colorsChoice(value_max + 1))

#indeldensity 
indelnub<-read.table("input/indel_chr_pos.txt",header=F)
indelnub<-separate(data = indelnub ,col = V1,into = c("chr","pos"),sep = "[:]+")
indelnub <- indelnub[!substr(indelnub$chr,1,3)=='CTG',]
indelnub <- arrange(indelnub,chr,pos)
#indelnub$window <- as.numeric(indelnub$pos)%/%9000
#indelnub <- indelnub %>%
#   group_by(chr,window) %>%
#   summarise(nub=length(window))
#indelnub$start <- indelnub$window*9000
#indelnub$end <- indelnub$start+9000
#indelnub <- indelnub[,c(1,4,5,3)]
#colnames(indelnub)<-c("chr","start","end","nub")


#value_max <- max(indelnub$nub)
#colorsChoice <- colorRampPalette(c('#FFFFFF', 'blue'))
#indel_color_assign <- colorRamp2(breaks = c(0:value_max), col = colorsChoice(value_max + 1))

#pi
pi <- read.table("input//pixy_pi.txt",header = T)
pi <- na.omit(pi)
pi <- pi[,c(2:5)]
pi<-pi[!substr(pi$chromosome,1,3)=='CTG',]



#adaptive SNP

apsnp<-read.table('input//2sd_selected_snp_name.txt')
apsnp<-separate(data = apsnp ,col = V1,into = c("chr","pos"),sep = "[:]+")
apsnp<-apsnp[!substr(apsnp$chr,1,3)=='CTG',]
apsnp$pos <- as.numeric(apsnp$pos)-75000;apsnp$posend <- as.numeric(apsnp$pos)+75000
apsnp$pos[apsnp$pos<0] <- 0

#adaptive SNP
apindel<-read.table('input//2sd_selected_indel_name.txt')
apindel<-separate(data = apindel ,col = V1,into = c("chr","pos"),sep = "[:]+")
apindel<-apindel[!substr(apindel$chr,1,3)=='CTG',]
apindel$pos <- as.numeric(apindel$pos)-75000;apindel$posend <- as.numeric(apindel$pos)+75000
apindel$pos[apindel$pos<0] <- 0

#adaptive gene
adgene <- read.csv("input//maf001_2sd_fdr001_308gene.csv")
adgene<-adgene[!substr(adgene$genename,1,3)=='CTG',]
adgene<-separate(data = adgene ,col = genename,into = c("chr","pos"),sep = "[.]+")
adgene<-adgene[,-2]

#colline 
colline <- read.table("input/GPH_LickedRegion.tab.txt")

icolline <- colline[!substr(colline$V1,1,5)==substr(colline$V4,1,5),] 
icolline <- rbind(icolline,colline[colline$V1==colline$V4,])
icolline$col <- "#FFFFA0"
icolline$h<- 0.2

colline <- colline[substr(colline$V1,6,6)=="A",]
colline <- colline[seq(1,nrow(colline),60),]
colline$col[substr(colline$V4,6,6)=="B"] <- "#E86E55"
colline$col[substr(colline$V4,6,6)=="C"] <- "#1E90BB"
colline$col[substr(colline$V4,6,6)=="D"] <- "#54AC78"
colline$col[is.na(colline$col)] <- "#FFFFA0"

colline$h[substr(colline$V4,6,6)=="B"] <- 0.1
colline$h[substr(colline$V4,6,6)=="C"] <- 0.15
colline$h[substr(colline$V4,6,6)=="D"] <- 0.2
colline$h[is.na(colline$h)] <- 0.2






#plot
pdf("colline.pdf",width = 21*0.8 , height = 29.7*0.8 )
circos.genomicInitialize(chrlen, plotType = c('axis', 'labels'))
#pi
pma <- max(pi$avg_pi)
pmi <- min(pi$avg_pi)
circos.track(track.height = 0.05,  bg.border = NA,ylim=c(pmi,pma),
             factor=pi$chromosome,x = pi$window_pos_1, y = pi$avg_pi,
             panel.fun = function(x, y) {
                circos.lines(x, y, col = "red") #默认折线图
             })
##genedensity
circos.trackHist(track.height = 0.05,bg.border = NA,ylim = c(-20,40),
                 bin.size=90000,area =T,
                 col = "#D2691E",border="#D2691E",
                 factor=gene$chr, x = gene$start
)

#circos.genomicTrackPlotRegion(
#   genenub, track.height = 0.05, stack = TRUE, bg.border = NA,
#   panel.fun = function(region, value, ...) {
#      circos.genomicRect(region, value, col = gene_color_assign(value[[1]]), border = NA, ...)
 #  } )
#snpdensity
circos.trackHist(track.height = 0.05,bg.border = NA,ylim = c(-250,500),
                 bin.size=90000,area =T,
                 factor=snpnub$chr, x = as.numeric(snpnub$pos), 
                 col = "#ff7f00",border="#ff7f00"
)
#circos.genomicTrackPlotRegion(
#   snpnub, track.height = 0.08, stack = TRUE, bg.border = NA,
 #  panel.fun = function(region, value, ...) {
#      circos.genomicRect(region, value, col = snp_color_assign(value[[1]]), border = NA, ...)
#   } )
#indeldensity
circos.trackHist(track.height = 0.05,bg.border = NA,ylim = c(-20,25),
                 bin.size=150000,area =T,
                 factor=indelnub$chr, x = as.numeric(indelnub$pos), 
                 col = "#fdbf6f",border="#fdbf6f"
)
#circos.genomicTrackPlotRegion(
#   indelnub, track.height = 0.08, stack = TRUE, bg.border = NA,
#   panel.fun = function(region, value, ...) {
#      circos.genomicRect(region, value, col = indel_color_assign(value[[1]]), border = NA, ...)
#   } )

#adgene
circos.genomicTrackPlotRegion(
   adgene, track.height = 0.05, stack = TRUE, bg.border = NA,
   panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "red", border = NA, ...)
   } )
#adsnp
circos.genomicTrackPlotRegion(
   apsnp, track.height = 0.05, stack = TRUE, bg.border = NA,
   panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "#FF7F50", border = NA, ...)
   } )
#adindel
circos.genomicTrackPlotRegion(
   apindel, track.height = 0.05, stack = TRUE, bg.border = NA,
   panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "orange", border = NA, ...)
   } )




#colline plot
#circos.track(ylim=c(0,1),track.height = 0.08,bg.border = NA,)
scolline <- colline[colline$col=="#54AC78",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h)


scolline <- colline[colline$col=="#1E90BB",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h)

scolline <- colline[colline$col=="#E86E55",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h)

scolline <- colline[colline$col=="#FFFFA0",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h)

bed1 <- cbind(icolline[,c(1,2,3)])
bed2 <- cbind(icolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=icolline$col,h=0.2)

dev.off()

circos.clear()


######
pdf("colline.pdf",width = 21*0.7 , height = 29.7*0.4 )
#png("colline.png",width = 1200,height = 1000,units = "px")
circos.genomicInitialize(chrlen, plotType = c('axis', 'labels'))


scolline <- colline[colline$col=="#FBD75F50",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h,lwd=0.01,border=NA)


scolline <- colline[colline$col=="#0050B750",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h,lwd=0.01,border=NA)

scolline <- colline[colline$col=="#DB5A6B50",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h,lwd=0.01,border=NA)


scolline <- colline[colline$col=="#808080",]
bed1 <- cbind(scolline[,c(1,2,3)])
bed2 <- cbind(scolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=scolline$col,h=scolline$h,lwd=0.01,border=NA)


bed1 <- cbind(icolline[,c(1,2,3)])
bed2 <- cbind(icolline[,c(4,5,6)])
colnames(bed1) <- colnames(bed2) <- c("chr","start",'end')
circos.genomicLink(bed1,bed2,col=icolline$col,h=icolline$h,lwd=0.5,border=NA)

dev.off()
