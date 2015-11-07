# Figure 4 #
library(pheatmap);
library(doMC);
setwd("/scratchLocal/shl2018/AML/ERRBS/epipoly2/")
sumt=read.fast("~/ntd/errbs/amlc/T1/T1_entropy.txt",header=T,sep="\t")
pheatmap(sumt[which(rowMeans(sumt[,-1:-3])< -30), -1:-3])

files=dir("~/ntd/errbs/amlc/", pattern="_entropy.txt", recursive = T, full.names = T)
tnsum=list()
for(f in files){
  tnsum[[basename(f)]]=read.fast(f,header=T,sep="\t")
}

n=2
pheatmap(tnsum[[n]][which(rowMeans(tnsum[[n]][,-1:-3]) < -60),-1:-3])
unionloci=foreach(n=1:5, .combine=c) %do% {
  loci=paste(tnsum[[n]][,1], tnsum[[n]][,2], tnsum[[n]][,3],sep="_")
  id=which(rowSums(tnsum[[n]][,-1:-3] < -90) >= 11)
  loci[id]
}
uniqueloci=unique(unionloci)
unionm=foreach(n=1:5, .combine=cbind) %do% {
  loci=paste(tnsum[[n]][,1], tnsum[[n]][,2], tnsum[[n]][,3],sep="_")
  tmp=rep(0, length(uniqueloci))
  id=which(rowSums(tnsum[[n]][,-1:-3] < -90) >= 11)
  tmp[uniqueloci %in% loci[id] ]=1
  tmp
}
colnames(unionm)=splitn(names(tnsum),"_",1)
pdf("~/ntd/errbs/amlc/entropy_heatmap.pdf",width=4,height=4)
pheatmap(unionm, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = "binary",
         c("grey","black"))
dev.off()

# use e90_11p_eloci.txt
files=dir("~/ntd/errbs/amlc",pattern="e90_11p_eloci.txt", full.names = T, recursive = T)
elocil=foreach(f = files) %do% {
  read.table(f, header=F,sep="\t",stringsAsFactors = F)$V1
}
unieloci=names(sort(table(unlist(elocil))))
elocim=foreach(i=1:5,.combine=cbind) %do% {
  tmp=rep(0, length(unieloci))
  tmp[unieloci%in%elocil[[i]]]=1
  tmp
}
colnames(elocim)=paste0("T",1:5)
elocim2=elocim[order(elocim[,1], elocim[,2], elocim[,3], elocim[,4], elocim[,5]),]
pdf("~/ntd/errbs/amlc/entropy_heatmap_short.pdf",width=3.5,height=2.5)
pheatmap(elocim2,color=c("cornsilk","deepskyblue"), cluster_rows = F, cluster_cols = F)
dev.off()
