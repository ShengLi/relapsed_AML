# Figure 3 # 
# Eloci shared
print(load("/Volumes/icb2/Awork/AMLpaper/manuscript_materials//epicluster/fig3_14N.shared.aug.rda"))
greens=c("#a1d99b","#41ab5d", "#006d2c")
dat4$variable=gsub("\\.", " ", dat4$variable)
p=ggplot(dat4, aes(ID, epiallele.shift,fill=variable))+ 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values=greens) +
  facet_grid(~cluster,scale="free",labeller=label_both, space="free_x") + gb2.theme("bottom") +
  scale_x_continuous(breaks=c(1, 138))
printtopdf(p, "~/awork/projects/paml/figures/Figure3/fig3.shared_14N_3_epi_green2_oct.pdf", width=8,height=5)

# SNV callings by merging three snv caller results.
library(VariantAnnotation)
getFinalSnvsInfo=function(pid){
  nid=gsub("-[12]$", "-N", pid)
  
  vvcf=readVcf(paste0(pid, "/", pid, "_", nid, ".varscan.snv.fltbyindel.Somatic.hc.vcf"),"hg19")
  mvcf=readVcf(paste0(pid, "/", pid, "_", nid, ".mutect.snv.Somatic.vcf"),"hg19")
  # mutect
  pid2=grep("N", colnames(geno(mvcf)$AD), value=T, invert=T)
  nid2=grep("N", colnames(geno(mvcf)$AD), value=T, invert=F)
  t.alt.count=round(unlist(geno(mvcf)$DP[,pid2])* unlist(geno(mvcf)$FA[,pid2]))
  t.depth=unlist(geno(mvcf)$DP[,pid2])
  t.ref.count=t.depth-t.alt.count
  n.alt.count=round(unlist(geno(mvcf)$DP[,nid2])* unlist(geno(mvcf)$FA[,nid2]))
  n.depth=unlist(geno(mvcf)$DP[,nid2])
  n.ref.count=n.depth-n.alt.count
  freq=unlist(geno(mvcf)$FA[,pid2])
  info=as.data.frame(rowRanges(mvcf))
  msnvs=data.frame(info[,1:2],t.depth, t.ref.count, t.alt.count, n.depth, n.ref.count, n.alt.count, freq, info[,7:8])
  mids=paste(msnvs$seqnames, msnvs$start, msnvs$REF, sep="_")
  # varscan
  t="TUMOR"
  n="NORMAL"
  t.depth=geno(vvcf)$DP[,t]
  t.ref.count=geno(vvcf)$RD[,t]
  t.alt.count=geno(vvcf)$AD[,t]
  n.depth=geno(vvcf)$DP[,n]
  n.alt.count=geno(vvcf)$AD[,n]
  n.ref.count=geno(vvcf)$RD[,n]
  freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ)[,t])/100
  info=as.data.frame(rowRanges(vvcf))
  vsnvs=data.frame(info[,1:2],t.depth, t.ref.count, t.alt.count, n.depth, n.ref.count, n.alt.count, freq, info[,7:8])
  vids=paste(vsnvs$seqnames, vsnvs$start, vsnvs$REF, sep="_")
  
  # merge varscan + mutect
  merged.vcf=rbind(msnvs, vsnvs[which(! vids %in% mids),])
  mergedids=paste(merged.vcf$seqnames, merged.vcf$start, merged.vcf$REF, sep="_")
  
  # use of somatic loci shared by at least two callers.
  snvloci=read.table(paste0(pid, "/", pid, ".somatic.snv.bed"), header=F, sep="\t", stringsAsFactors = F)
  somaticids=paste(snvloci$V1, snvloci$V2, snvloci$V4, sep="_")
  
  finalsnv=merged.vcf[mergedids %in% somaticids,]
  
  write.table(finalsnv, file = paste0(pid, "/", pid, ".somatic.snv.info.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  save(finalsnv,file=paste0(pid, "/", pid, ".somatic.snv.info.rda"))
  finalsnv
}

library(doMC)
registerDoMC(cores=5);

snvs=foreach(id = dir(pattern="-[12]$")) %dopar% {
  snv=getFinalSnvsInfo(id);
}
names(snvs)=dir(pattern="-[12]$")
save(snvs, file="merged.snvs.rda")

# high confidence indels filterings
vindelvcf=readVcf(paste0(pid, "/", pid, "_", nid, ".varscan.indel.Somatic.hc.vcf"),"hg19")
vindelvcf1=vindelvcf[unique(findOverlaps(rowRanges(vindelvcf), ref.exon.gr)@queryHits)]
# varscan
t="TUMOR"
n="NORMAL"
t.depth=geno(vindelvcf1)$DP[,t]
t.ref.count=geno(vindelvcf1)$RD[,t]
t.alt.count=geno(vindelvcf1)$AD[,t]
n.depth=geno(vindelvcf1)$DP[,n]
n.alt.count=geno(vindelvcf1)$AD[,n]
n.ref.count=geno(vindelvcf1)$RD[,n]
freq=as.numeric(gsub("%", "", geno(vindelvcf1)$FREQ)[,t])/100
info=as.data.frame(rowRanges(vindelvcf1))
vindels=data.frame(info[,1:2],t.depth, t.ref.count, t.alt.count, n.depth, n.ref.count, n.alt.count, freq, info[,7:8])
vids=paste(vindels$seqnames, vindels$start, vindels$REF, sep="_")

# counts snvs:
setwd("/Volumes/hippo/home/aml/Exome/Sep2015/pairedAML")
pids=unique(gsub("-[12N]$", "",dir(pattern="-")))
pids1=pids[!pids %in% c("SGUR4", "SGUP4", paste0("PAML-", c(77,80, 89, 96, 105)))]
library(reshape)
library(doMC)
df1=foreach(i = pids1, .combine=rbind)%do%{
  loci=foreach(stage=1:2) %do% {
    snvloci=read.table(paste0(i, "-", stage, "/", i, "-", stage, ".somatic.snv.bed"), stringsAsFactors=F);
    somaticids=paste(snvloci$V1, snvloci$V2, snvloci$V4, sep="_")
  };
  d=length(setdiff(loci[[1]], loci[[2]]));
  r=length(setdiff(loci[[2]], loci[[1]]));
  dr=length(intersect(loci[[1]], loci[[2]]));
  data.frame(ID=i, type='snv', D=d, DR=dr, R=r)
}

indeldf1=foreach(i = pids1, .combine=rbind)%do%{
  loci=foreach(stage=1:2) %do% {
    snvloci=read.table(paste0(i, "-", stage, "/", i, "-", stage, "_", i, "-N", ".varscan.indel.Somatic.hc.exon.vcf"), stringsAsFactors=F);
    somaticids=paste(snvloci$V1, snvloci$V2, snvloci$V4, sep="_")
  };
  d=length(setdiff(loci[[1]], loci[[2]]));
  r=length(setdiff(loci[[2]], loci[[1]]));
  dr=length(intersect(loci[[1]], loci[[2]]));
  data.frame(ID=i, type='indel',D=d, DR=dr, R=r)
}

mutationdf=ddply(rbind(df1, indeldf1), .(ID),function(X){colSums(X[,3:5])})
mutationdf1=melt(mutationdf)
ids=gsub("SG", "",gsub("PAML","HU", as.character(mutationdf$ID)))
i=grep("HU",invert=T, ids)
ids2=c(ids[-i], paste(substr(ids[i], 1,2), substr(ids[i], 3, 4), sep="-"))
names(ids2)=mutationdf$ID
mutationdf$cluster=cl[ids2]
mutationdf1$cluster=cl[ids2[as.character(mutationdf1$ID)]]
mutationdf2=cast(mutationdf1)

# Order as epigenetic
ids=unique(dat4$sample)
ids2=gsub("HU-","PAML-", ids)
ids2[grep("PAML", invert=T, ids2)]=gsub("-","", paste0("SG",ids2[grep("PAML", invert=T, ids2)]))
ordered.id=ids2[ids2 %in% unique(as.character(mutationdf1$ID))]

mutationdf1$ID=factor(as.character(mutationdf1$ID), levels = ordered.id)
mutationdf1$types="Diagnosis-specific"
mutationdf1$types[mutationdf1$variable=="DR"]="Shared"
mutationdf1$types[mutationdf1$variable=="R"]="Relapse-specific"
mutationdf1$types=factor(mutationdf1$types, levels=c("Diagnosis-specific", "Shared", "Relapse-specific"))


blues=c("#74a9cf","#0570b0", "#034e7b")
p=ggplot(mutationdf1,aes(x=as.integer(ID), y=value, fill=types)) + geom_bar(stat='identity', position="fill") +
  facet_grid(~cluster,scale="free_x",labeller=label_both, space="free_x")  +
  scale_fill_manual(values=blues) + gb2.theme("bottom") + xlab("ID")+
  ylab("proportion of somatic mutations") +scale_x_continuous(breaks=c(1, 48))
printtopdf(p,"~/awork/projects/paml/figures/Figure3/mutations.burden.snv.indel2.Oct.48.pdf", width=8,height=5)
save(mutationdf1, df1, indeldf1, file="mutationdf1.48.rda")
save(cl, file="~/awork/projects/paml/figures/Figure3/cl.rda")
