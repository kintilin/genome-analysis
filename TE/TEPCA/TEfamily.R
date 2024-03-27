chr <- read.table("chr.list",header=F)
TE <- read.table("total_TEfamily.list",header=F)
gff <- read.table("../genome.fa.mod.EDTA.intact.gff3",header=F)

result <- data.frame(matrix(nrow=nrow(chr),ncol=nrow(TE)))
colnames(result) <- TE$V1
rownames(result) <- chr$V1
#head(result)
name <- lapply(strsplit(gff$V9,split=";"), function(x) substr(x[substr(x,1,4)=="Name"],6,16))
gffchr <- gff$V1
#head(gffchr)
#length(gffchr)
#length(name)
for (i in rownames(result) ){
  for (j in colnames(result)){
    result[i,j] <- sum(gffchr==i & name==j)  
}
}

result <- result[,colSums(result)!=0]
write.csv(result,'TEnub.csv',quote = F)



library(FactoMineR)
library(factoextra)
library(ggplot2)
tenub<-read.csv("TEnub0.csv",header=T,row.names=1)
tenub<-tenub[,colSums(tenub)!=0]
pcaBio<- PCA(tenub,scale.unit = T)
#summary(pcaBio)

kipbio<-pcaBio$eig[,3]
pdata <- data.frame(pcaBio$ind$coord[,1:3])
ggplot(pdata,aes(x=Dim.1,y=Dim.2))+
   geom_text(aes(label=rownames(pdata)))



