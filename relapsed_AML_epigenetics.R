# Figure 1 #

# survival analysis
# EPM
load("dn_epm_RTIME.rda")
library(survival)
survdiff(Surv(RTIME) ~ (EPM > median(EPM)), dn.epm.rtime) 

par(mar=c(4,4,4,1))
plot(survfit(Surv(RTIME/365) ~ (EPM > median(EPM)), dn.epm.rtime), col=c("black","red"),
     xlab="year",ylab="relapse-free probability", xlim=range(dn.epm.rtime$RTIME/365))
legend("topright", col=c("black","red"), c("Low EPM","High EPM"), lty=1, box.lwd=0, bty="n")

# Somatic Mutations
par(mar=c(4,4,4,1))
plot(survfit(Surv(rtime/365) ~ (D > median(D)), b138.num.octrtime), col=c("green","blue"),
     xlab="year",ylab="relapse-free probability", main="All AML (n=48)\np-value = 0.688", xlim=range(b138.num.octrtime$rtime/365))
legend("topright", col=c("green","blue"), c("Low MUT","High MUT"), lty=1, box.lwd=0, bty="n")

# Figure 2 # 

setwd("/Volumes/icb2/Awork/AMLpaper/epiloci/")
load("~/Google Drive/MasonLab/AML/manuscript_materials/b138.aug_rtimeupdate_137patient.rda")
AML147DR=read.table("survival/DR_eloci_147_Aug18.txt",header=F,sep="\t",stringsAsFactors=F)
AML147DN=read.table("survival/AML147_DN.txt",header=F,sep="\t",stringsAsFactors=F)
AML147RN=read.table("survival/AML147_RN.txt",header=F,sep="\t",stringsAsFactors=F)

AML144DN9=read.table("/Volumes/aml2/home/ntd/errbs/epipoly/AML144_DN_9N_Sep6.txt",header=F,sep="\t",stringsAsFactors=F)
AML144RN9=read.table("/Volumes/aml2/home/ntd/errbs/epipoly/AML144_RN_9N_Sep6.txt",header=F,sep="\t",stringsAsFactors=F)

print(load("../manuscript_materials/epicluster/epi_cl_aug.rda"))

header=c("ID", paste0("e", seq(90,50,length=5)),"loci")
dat=rbind(AML147DN, AML147RN, AML144DN9, AML144RN9)
colnames(dat)=header
dat$sample=paste(splitn(dat[,1], "-", 1), splitn(dat[,1],"-",2), sep="-" )
dat$stage=splitn(dat[,1], "-|_", 3)
dat$N=splitn(dat[,1], "_", 2)
dat$epm=dat$e90/dat$loc * 10^6
#dat1=ddply(dat, .(sample, stage), function(x)data.frame(epm90=mean(x$epm)))
dat1=ddply(dat, .(sample, stage), function(x)data.frame(epm90=mean(x$epm), epm80=mean(x$e90/x$loc*10^6),  epm70=mean(x$e70/x$loc*10^6)))
dat2=dat1[dat1$sample %in% c(b138.aug$sample, "UR-2"),]
dat2$cl=paste("Cluster", cl[dat2$sample])
save(dat2, file="~/awork/projects/paml/figures/Figure2/epmByCl.rda")
p=ggplot(dat2, aes(x=ifelse(stage==1, "D vs. N", "R vs. N"), y=log10(epm90), fill=stage)) + geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+ facet_grid(.~cl)+
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[2:3]) + xlab(NULL) + ggtitle("All AML (n=138)") 
printtopdf(p, "../manuscript_materials/aug2015/Figure2_log10epm90_epicluster.pdf")

p=ggplot(dat2, aes(x=ifelse(stage==1, "D vs. N", "R vs. N"), y=epm90, fill=stage)) + geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+ facet_grid(.~cl)+
  gb2.theme("none") + ylab("EPM") + scale_fill_manual(values=cbPalette[2:3]) + xlab(NULL) + ggtitle("All AML (n=138)") 
printtopdf(p, "../manuscript_materials/aug2015/Figure2_epm90_epicluster.pdf")
pdf("../manuscript_materials/aug2015/Figure2_epm90_epicluster_pval.pdf")
pval.by.stage=ddply(dat2, .(stage), function(X){data.frame(not.cl=1:3, pval=sapply(1:3, function(i) wilcox.test(epm90~cl, X[gsub("Cluster ","",X$cl)!=i,])$p.value))})
pval.by.cl=ddply(dat2, .(cl), function(X){data.frame(pval=wilcox.test(epm90~stage, X)$p.value)})
textplot(pval.by.stage)
textplot(pval.by.cl)
dev.off()
DN.epm90.aug=dat2$epm90[dat2$stage==1]
names(DN.epm90.aug)=dat2$sample[dat2$stage==1]
b138.aug$DN.epm90.aug=DN.epm90.aug[b138.aug$sample]

b138.aug$cl.aug=cl[b138.aug$sample]
summary(coxph(Surv(relapsetime.aug)~ cl.aug, b138.aug)) # pvalue=0.365
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=1),]) #p= 0.589 
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=2),]) #p= 0.317 
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=3),]) #p= 0.704 

# cluster and age/gender
anova(lm(dage~as.factor(cl.aug), b138.aug)) # pvalue=0.4413

p=ggplot(b138.aug, aes(x=factor(cl.aug), y=dage, fill=factor(cl.aug))) + geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+
  gb2.theme("none") + ylab("Age") + scale_fill_manual(values=cbPalette[5:7]) + xlab("Cluster") + ggtitle("All AML (n=138)\np-value=0.4413") 
printtopdf(p, pdf="aug2015/epm_vs_age.pdf",width=4,height=4)

anova(lm(sex~as.factor(cl.aug), b138.aug)) # pvalue=0.4413
chisq.test(table(b138.aug$sex, b138.aug$cl.aug))

dat.aug=ddply(b138.aug, .(cl.aug), function(X)data.frame(table(X$sex)))
p=ggplot(dat.aug, aes(x=Var1, y=Freq,fill=factor(cl.aug))) + 
  geom_bar(stat="identity", position="fill") + 
  gb2.theme("none") + ylab("Proportion of patients\nfrom three clusters") + scale_fill_manual(values=cbPalette[5:7]) + xlab("Gender") + ggtitle("All AML (n=138)\np-value=0.01048") 
printtopdf(p, pdf="aug2015/epm_vs_sex.pdf",width=4,height=4)


# add dr
colnames(AML147DR)=c("ID", paste0("e", seq(90,50,length=5)),"loci")
dat3=data.frame(sample=AML147DR$ID, stage=3, epm90=1e6*AML147DR$e90/AML147DR$loci, cl=paste("Cluster", cl[AML147DR$ID]))
dat4=dat3[dat3$sample %in% unique(dat2$sample),]
drn.dat=rbind(dat2, dat3[dat3$sample %in% unique(dat2$sample),])

# get the order by D vs NBM
dn.epm90=dat2$epm90[which(dat2$stage==1)]
names(dn.epm90)=dat2$sample[which(dat2$stage==1)]
od.sample=names(dn.epm90)[order(dn.epm90)]
dat2$ID=factor(dat2$sample, levels=od.sample)
dat4$ID=factor(dat4$sample, levels=od.sample)
dat4$epm90[is.na(dat4$epm90)]=0
p=ggplot( dat4, aes(x=as.numeric(ID), y=log10(epm90+1), fill=as.factor(stage))) + 
  geom_bar(stat="identity") + xlab("sample index")  +scale_x_continuous(breaks=c(1, 50,100, 138))+
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[4]) + ylim(c(0,5.1))
printtopdf(p,"~/awork/projects/paml/figures/Figure2/figure2d.pdf", width=5,height=2.5)

p=ggplot( dat2[dat2$stage==2,], aes(x=as.numeric(ID), y=log10(epm90+1), fill=as.factor(stage))) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[3]) + ylim(c(0,5.1))
printtopdf(p,"~/awork/projects/paml/figures/Figure2/figure2b.pdf", width=5,height=2.5)

p=ggplot( dat2[dat2$stage==1,], aes(x=as.numeric(ID), y=log10(epm90+1), fill=as.factor(stage))) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[2]) + ylim(c(0,5.1))
printtopdf(p,"~/awork/projects/paml/figures/Figure2/figure2a.pdf", width=5,height=2.5)

wilcox.test(dat2[dat2$stage==1,"epm90"],dat2[dat2$stage==2,"epm90"]) #p-value = 0.00463

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

# Supp Figure 6 #

# Clonal analysis using sciClone

library(sciClone)
pids=unique(gsub("-[12N]$", "",dir(pattern="-")))
pids1=pids[!pids %in% paste0("PAML-", c(80, 89, 96, 105))]
library(sciClone)
registerDoMC(cores=5)
foreach(p = pids1) %dopar%{
  Tsnvs=foreach(i = 1:2) %do% {
    sample.id=paste(p,i,sep="-");
    snvf=paste0(sample.id,"/", sample.id, ".somatic.snv.info.rda");
    load(snvf);
    finalsnv[,c(1,2,4,5,9)]
  }
  names(Tsnvs)=paste(p, 1:2, sep="-")
  acnvs=foreach(i = 1:2) %do% {
    sample.id=paste(p,i,sep="-");
    cnvf=paste0(sample.id,"/", sample.id, "_mergedcnv.rda");
    load(cnvf)
    cnv
  }
  names(acnvs)=paste(p, 1:2, sep="-")
  xids=paste(p,1:2,sep="-");
  print(p)
  #Tsnvs=gm.asnvs.exons[xids]
  # filling the SNVs gaps by making psudo zero vaf postion
  Tsnvsr=do.call(rbind, Tsnvs)
  Tsnvsu=unique(Tsnvsr[,1:2])
  Tsnvsu$refn=50
  Tsnvsu$altn=0
  Tsnvsu$freq=0
  rownames(Tsnvsu)=paste(Tsnvsu$seqnames, Tsnvsu$start,sep=".")
  Tsnvs2=lapply(Tsnvs, function(X){
    x=setdiff(rownames(Tsnvsu), paste(X$seqnames, X$start, sep="."))
    colnames(Tsnvsu)=colnames(X)
    rbind(X, Tsnvsu[x,])
  })
  for(i in 8){
    tryCatch({
      sc2=sciClone(Tsnvs2, copyNumberCalls = acnvs[xids],sampleNames =xids,cnCallsAreLog2 = F,minimumDepth = 50, maximumClusters = i)
      sc.plot2d(sc2,outputFile = paste0(xids,"/",xids,"_sc2clo.exons.max", i, ".pdf"))
      save(sc2, file=paste0(xids[1],"/",xids[1],"_sc2clo.exons.max", i, ".rda"))
      writeClusterTable(sc2, paste0(xids[1],"/",xids[1],"_sc2clo.exons.max", i, ".txt"))
    }, error=function(e){})
  }
}

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

# Figure 5 # 
getgr=function(x, n=8){
  x$start=as.numeric(as.character(x$start))
  x$end=as.numeric(as.character(x$end))
  gr=GRanges(Rle(x$chr),IRanges(x$start, x$end), value=x[,n])
  gr
}

elociexpr2=function(id, prpkmdiff2, en=-80, nen=-20, p1=refseq.promoters, p2=refseq.promoters2){
  #/Volumes/hippo/zenodotus/dat01/melnick_lab_scratch_paired_aml/shl2018
  e=read.table(paste0("~/.aml/ERRBS/epipoly/Jan302014/", id,"/",id,".entropy.gz"),header=T,sep="\t",stringsAsFactors=F)
  egr=getgr(e[e$entropy < en,], 8)
  o=findOverlaps(egr,p1)
  egs=(values(p1)[o@subjectHits,"gene"])
  useg=intersect(rownames(prpkmdiff2), unique(egs))
  ent=data.frame(gene=egs, ent=values(egr)[o@queryHits,1])
  entt=ddply(ent[ent$gene %in% useg,], .(gene), summarize, value=mean(ent))
  expr=prpkmdiff2[as.character(entt$gene),id]
  
  egr2=getgr(e, 8)
  o2=findOverlaps(egr2, p2)
  locigenes=data.frame(gene=values(p2)$gene[o2@subjectHits], entropy=e$entropy[o2@queryHits])
  locigenes2=ddply(locigenes, .(gene), summarize, value=min(entropy))
  useg2=intersect(rownames(prpkmdiff2), locigenes2$gene[locigenes2$value > nen])
  locigenes3=locigenes2[locigenes2$gene%in%useg2,]
  expr2=prpkmdiff2[as.character(locigenes3$gene),id]
  
  dat=rbind(data.frame(id=id, type="eloci", expr=expr, entt ),
            data.frame(id=id, type="noneloci", expr=expr2,locigenes3))
  return(dat)
}

me=foreach(id = colnames(prpkmdiff2)[-30], .combine=rbind) %dopar% {
  df=elociexpr2(id, prpkmdiff2, -70, -10, refseq.promoters, refseq.promoters); df
}

gmme=foreach(id = colnames(pgmrpkmdiff2)[c(-2,-17)], .combine=rbind) %dopar% {
  elociexpr2(id, pgmrpkmdiff2, -70, -10, refseq.promoters, refseq.promoters)
}
ref=c("UP-12", "UP-13", "UP-16")
meee=rbind(me[!as.character(me$id) %in% ref,],gmme)
meee$ID=newID2[as.character(meee$id)]

meee2=meee[meee$value < -80 | meee$value > -2, ]
meee3=meee2[meee2$ID%in% names(which(table(meee2$ID[meee2$type=="eloci"]) > 30) ), ]
p=ggplot(meee3, aes(x=ID, y=abs(expr), colour=type)) + geom_boxplot(outlier.size = 0) + ylim(c(0, 4)) +
  gb3.theme()  + ylab("absolution gene expression log fold change")
printtopdf(p, "eloci_expr_75_5.pdf",width=8,height=6)

p=ggplot(meee3[meee3$ID %in% c("AML_08", "AML_31", "PAML_02", "PAML_64"),], aes(x=expr, colour=type))+
  geom_rect(aes(xmin = -1, xmax = 1, ymin = 0, ymax = 1), fill="grey88", colour="grey88", alpha = 0.05)+
  geom_density() + gb2.theme() + xlab("gene expression log fold change") + facet_wrap(~ID)
printtopdf(p,"AML_4p.density.pdf",width=6,height=6)

p=ggplot(ddply(meee3,.(ID, type), summarize, value=var((expr))), aes(x=type,y=value, colour=type))+
  geom_boxplot(outlier.size = 0)+ ylab("variance of gene expression\nlog fold change") +
  gb2.theme("none") + ggtitle(expression("p-value = 5.96 x 10"^-8))
printtopdf(p, "eloci30_expr_75_5_var_box.pdf",width=4,height=4)
