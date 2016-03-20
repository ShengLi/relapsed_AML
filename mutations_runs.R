# somatic mutations and indels analysis
# calculate the number of somatic mutations based on updated information on exome-capture inform.
source("mutations_functions.R")

## see methods section for calling somatic mutations by varscan/mutect/somaticSniper
# sample name example: AML-100-1 is patient ID AML-100 at diagnosis; 
# sample name example: AML-100-2 is patient ID AML-100 at relapse;
# sample name example: AML-100-N is patient ID AML-100 matched normal.
# varscan somatic snvs files extensions are ".varscan.snv.fltbyindel.Somatic.hc.eff.vcf.gz"
# example: AML-100-1/AML-100-1_AML-100-N.varscan.snv.fltbyindel.Somatic.hc.eff.vcf.gz
# mutect somatic snvs files extensions are ".mutect.snv.Somatic.eff.vcf.gz"
# example: AML-100-1/AML-100-1_AML-100-N.mutect.snv.Somatic.eff.vcf.gz
# somaticsniper somatic snvs files extensions are ".somaticsniper.snv.Somatic.eff.vcf.gz"
# example: AML-100-1/AML-100-1_AML-100-N.somaticsniper.snv.Somatic.eff.vcf.gz
merge.snv.full=function(pid, exome.gr){
  nid=gsub("-[12]$", "-N", pid)
  varscan.vcffile=paste0(pid, "/", pid, "_", nid, ".varscan.snv.fltbyindel.Somatic.hc.eff.vcf.gz")
  mutect.vcffile=paste0(pid, "/", pid, "_", nid, ".mutect.snv.Somatic.eff.vcf.gz")
  ss.vcffile=paste0(pid, "/", pid, "_", nid, ".somaticsniper.snv.Somatic.eff.vcf.gz")
  
  # get snvs and filter by n.freq < 0.05, n.depth & t.depth >= 10
  varscan.snvs=read.varscan.vcf(varscan.vcffile,exome.gr)
  varscan.snvs.flt=varscan.snvs[varscan.snvs$n.freq < 0.05 & 
                                  varscan.snvs$t.depth >=10 & 
                                  varscan.snvs$n.depth >=10,]
  mutect.snvs=read.mutect.vcf(mutect.vcffile,exome.gr)
  mutect.snvs.flt=mutect.snvs[mutect.snvs$n.freq < 0.05 & 
                                mutect.snvs$t.depth >=10 & 
                                mutect.snvs$n.depth >=10,]
  ss.snvs=read.ss.vcf(ss.vcffile,exome.gr)
  ss.snvs.flt=ss.snvs[ss.snvs$n.freq < 0.05 
                      & ss.snvs$t.depth >=10 
                      & ss.snvs$n.depth >=10,]
  
  v.id=get.snvs.id(varscan.snvs.flt)
  m.id=get.snvs.id(mutect.snvs.flt)
  s.id=get.snvs.id(ss.snvs.flt)
  id.freq=table(c(s.id,m.id,v.id))
  somaticids=names(which(id.freq>1))
  
  # merge varscan + mutect
  merged.vcf=rbind(mutect.snvs.flt, 
                   varscan.snvs.flt[which(v.id%in% setdiff( v.id, m.id)),],
                   ss.snvs.flt[which(s.id%in% setdiff( s.id, unique(c(m.id, v.id)))),])
  mergedids=paste(merged.vcf$seqnames, merged.vcf$start, merged.vcf$REF, merged.vcf$ALT, sep="_")
  merged.vcf$freq= id.freq[mergedids]
  merged.vcf$varscan = 0
  merged.vcf$varscan[which(mergedids %in% v.id)] = 1
  merged.vcf$mutect = 0
  merged.vcf$mutect[which(mergedids %in% m.id)] = 1
  merged.vcf$somaticsniper = 0
  merged.vcf$somaticsniper[which(mergedids %in% s.id)] = 1
  
  merged.vcf
  
}

# call snvs 
library(doMC)
registerDoMC(cores=25); 
p.snvs.full=foreach(id = dir(pattern="-[12]$")) %dopar% {
  snv=merge.snv.full(id, exome.gr);
}
names(p.snvs.full)=dir(pattern="-[12]$")

# convert snvs list into one data.frame
p.snvs.full.df=foreach(id=names(p.snvs.full), .combine=rbind) %do% {
  i=p.snvs.full[[id]]; 
  i$id=id; 
  i
}
p.snvs.full.df$pid=gsub("-[12]$","",p.snvs.full.df$id)
p.snvs.full.df$sample=p.snvs.full.df$pid
p.snvs.full.df$pub.id=p.snvs.full.df$pid
p.snvs.full.df$stage=substr(p.snvs.full.df$id, nchar(p.snvs.full.df$id), nchar(p.snvs.full.df$id))

# indels
p.indels.full=foreach(id = dir(pattern="-[12]$"), .combine=rbind) %dopar% {
  nid=gsub("-[12]$", "-N", id)
  vindels.vcffile=paste0(id, "/", id, "_", nid,".varscan.indel.Somatic.hc.eff.vcf")
  indels=read.vindels.vcf(vindels.vcffile, exome.gr);
  data.frame(id=id, indels)
}
save(p.indels.full, file="p.indels.full.rda")

p.indels.full$pid=gsub("-[12]$","",p.indels.full$id)
p.indels.full$sample=p.indels.full$pid
p.indels.full$pub.id=p.indels.full$pid
p.indels.full$stage=substr(p.indels.full$id, nchar(p.indels.full$id), nchar(p.indels.full$id))

# combine snvs and indels locations
cov.loc=unique(rbind(p.snvs.full.df[,c(1,2,2)], p.indels.full[,c(2,3,3)]))
write.table(cov.loc, quote=F, sep="\t", row.names=F, col.names=F, 
            file="p.snvs.indels.loc.txt")

# use "bedtools intersect" to calculate the coverage of each locus in p.snvs.indels.loc.txt file for each exome-seq bam file.
# the extension of the output of bedtools file is ".muts.bed.gz"
# read the coverage based bedtools output for each exome-seq bam file.
cov.files=dir("./",pattern="muts.bed.gz", full.names = T, recursive = T)
cov.df=foreach(f = cov.files, .combine=cbind) %do% {
  read.table(f, header=F, sep="\t", stringsAsFactors = F)$V5
}
colnames(cov.df)=basename(dirname(cov.files))
cov.f1=read.table(cov.files[1], header=F, sep="\t", stringsAsFactors = F)
rownames(cov.df)=paste(cov.f1$V1,cov.f1$V2,sep="_")
# filter loci based on the number of exome-seq bam file with the coverage equal or greater than 10 reads.
cov.f1.2=cov.f1[rowSums(cov.df >=10 ) >= 108,1:2]
colnames(cov.f1.2)=c("seqnames","start")
# filter snvs and indels based on coverage.
p.indels.full.bycov=merge(p.indels.full, cov.f1.2)
p.indels.full.df.bycov=ddply(p.indels.full.bycov, .(sample, pub.id, stage), summarize, indels=length(sample))
p.snvs.full.df.bycov=merge(p.snvs.full.df[p.snvs.full.df$freq>1 & p.snvs.full.df$t.freq>=0.05,], cov.f1.2)
p.snvs.full.df1.bycov=ddply(p.snvs.full.df.bycov, .(sample, pub.id, stage), summarize, snvs=length(sample))
p.s.i.full.bycov=merge(p.snvs.full.df1.bycov, p.indels.full.df.bycov, all=T)
# get epigenetic cluster information based on the Kmeans clustering
# input of Kmeans clustering is proportion of eloci that are 
# 1. specific to diagnosis stage
# 2. shared by both stages
# 3. specific to relapse stage
# save the cluster information in object "cl"
# snvs and indels numbers after filter by coverage 
# stage change on VAF
col.id=c("sample", "pub.id","seqnames","start","stage", "t.freq")
p.mut.full.bycov=rbind(data.frame(type="indels", p.indels.full.bycov[,col.id]),
                       data.frame(type="snvs", p.snvs.full.df.bycov[,col.id]))
p.mut.full.bycov$cluster=cl[p.mut.full.bycov$sample]

# get the bed file 'overlapping_third.bed' with the regions overlap with the third exome-capture targeted file.
e3=read.table('overlapping_third.bed', header=F, sep="\t", stringsAsFactors = F)
e3$V1=paste("chr", e3$V1, sep="")
e3.gr=getgr(e3)
p.mut.full.bycov1=rbind(data.frame(type="indels", p.indels.full.bycov[,c(col.id,"genes")]),
                        data.frame(type="snvs", p.snvs.full.df.bycov[,c(col.id,"genes")]))
p.mut.full.bycov1$cluster=cl[p.mut.full.bycov1$sample]
p.mut.full.bycov.gr=getgr(p.mut.full.bycov1[, c(4,5,5)])
e3.o=findOverlaps(p.mut.full.bycov.gr, e3.gr)
p.mut.bycov.e3=p.mut.full.bycov1[unique(e3.o@queryHits),]
# uncovered by e3 but with at least 140 samples covered were also kept.
p.mut.bycov.e3.uncov=p.mut.full.bycov1[-unique(e3.o@queryHits),]
uncov.ids=paste(p.mut.bycov.e3.uncov$seqnames,p.mut.bycov.e3.uncov$start, sep="_")
p.mut.bycov.e3.uncov.140.ids = which(rowSums(cov.df[uncov.ids,] >= 10) >= 140)
p.mut.bycov.e3.uncov.140=p.mut.bycov.e3.uncov[p.mut.bycov.e3.uncov.140.ids,]
# calculate the num of mutations by: 
# fall into the all three files 
# or fall into the first two files from before, and 
# also cover by at least 140 patients with at least 10 reads in third exome bed files
p.mut.e3=rbind(p.mut.bycov.e3,p.mut.bycov.e3.uncov.140)

p.mut.num.e3=ddply(p.mut.e3, .(pub.id, stage, cluster), summarize, value=length(t.freq))
# supp table 6
ddply(p.mut.num.e3, .(stage, cluster), summarize, medianmut=median(value))

# figure 3c
#p.mut.full.bycov.num=ddply(p.mut.full.bycov, .(pub.id, cluster, stage), summarize, value=length(t.freq))
p=ggplot(subset(p.mut.num.e3, stage ==1), aes(x=paste("cluster",cluster), y=log10(value+1), fill=as.factor(cluster))) +
  geom_violin() + stat_summary(fun.y=median.quartile,geom='point') + xlab(NULL) +  gb2.theme("none")+ 
  scale_fill_manual(values=cbPalette[5:7]) +ylab("number of somatic mutations (log10)")
printtopdf(p, pdf="~/awork/projects/paml/data4reviewerNov30/figures/figure3c.pdf",width=5,height=5)
wilcox.test(value~ cluster, subset(p.mut.num.e3, cluster != 2 &stage ==1), exact=F)
# figure 3d 
p=ggplot(subset(p.mut.num.e3, stage ==2), aes(x=paste("cluster",cluster), y=log10(value+1), fill=as.factor(cluster))) + 
  geom_violin() + stat_summary(fun.y=median.quartile,geom='point') + xlab(NULL) + gb2.theme("none")+ 
  scale_fill_manual(values=cbPalette[5:7]) +ylab("number of somatic mutations (log10)")
printtopdf(p, pdf="~/awork/projects/paml/data4reviewerNov30/figures/figure3d.pdf",width=5,height=5)
wilcox.test(value~ cluster, subset(p.mut.num.e3, cluster != 2 &stage ==2), exact=F)

p.mut.bycov.cast=cast(p.mut.e3, sample+pub.id+ seqnames+start + type + cluster  ~ stage, value = "t.freq")
# figure 3b
nd=7; nr=nd+1
colnames(p.mut.bycov.cast)[nd:nr]=c("Diagnosis","Relapse")
p.mut.bycov.cast.dsr=ddply(p.mut.bycov.cast, .(pub.id, sample, cluster), 
                           function(x){
                             S=sum(!is.na(x[,nd]) & !is.na(x[,nr]))
                             D=sum(!is.na(x[,nd]) & is.na(x[,nr]))
                             R=sum(is.na(x[,nd]) & !is.na(x[,nr]))
                             data.frame(D=D, S=S, R=R)})
dsr.melt=melt(p.mut.bycov.cast.dsr, id.vars=c("pub.id","sample","cluster"))
# significance for figure 3b
with(subset(p.mut.bycov.cast.dsr, cluster==1),wilcox.test(D/(D+S+R), R/(D+S+R), paired = T)) 
with(subset(p.mut.bycov.cast.dsr, cluster==2),wilcox.test(D/(D+S+R), R/(D+S+R), paired = T)) 
with(subset(p.mut.bycov.cast.dsr, cluster==3),wilcox.test(D/(D+S+R), R/(D+S+R), paired = T)) 

# figure 3e
p.mut.bycov.cast.1=p.mut.bycov.cast
p.mut.bycov.cast.1[is.na(p.mut.bycov.cast.1[,"Diagnosis"]),"Diagnosis"]=0
p.mut.bycov.cast.1[is.na(p.mut.bycov.cast.1[,"Relapse"]),"Relapse"]=0
# Order as epigenetically determined clusters
print(load("fig3_14N.shared.aug.rda"))
ids=unique(dat4[,c(1,8)])
ids2=ids[,2]
names(ids2)=ids[,1]
dsr.melt$ID=ids2[dsr.melt$sample]
dsr.melt$ID=factor(as.character(dsr.melt$ID), levels = sort(unique(dsr.melt$ID)))
dsr.melt$types="Diagnosis-specific"
dsr.melt$types[as.character(dsr.melt$variable)=="S"]="Shared"
dsr.melt$types[as.character(dsr.melt$variable)=="R"]="Relapse-specific"
dsr.melt$types=factor(dsr.melt$types, levels=c("Diagnosis-specific", "Shared", "Relapse-specific"))

blues=c("#74a9cf","#0570b0", "#034e7b")
p=ggplot(dsr.melt,aes(x=as.integer(ID), y=value, fill=types)) + geom_bar(stat='identity', position="fill") +
  facet_grid(~cluster,scale="free_x",labeller=label_both, space="free_x")  + 
  scale_fill_manual(values=blues) + gb2.theme("bottom") + xlab("ID")+ 
  ylab("proportion of somatic mutations") +scale_x_continuous(breaks=c(1, 48))
printtopdf(p,"figure3b.pdf", width=8,height=5)

# supp fig3
# dncli2.Feb2016 is a data.frame contain the 
# 1. PubName (sample name), 2. epm90 (EPM with ∆S < -90), and 3. rtime (days from diagnosis to relapse)
p.mut.full.bycov.num.epm=merge(p.mut.num.e3, dncli2.Feb2016[,c("PubName","epm90", "rtime")], 
                               by.x="pub.id", by.y="PubName")

# median epm90 of 48 AML patients
ggplot(p.mut.full.bycov.num.epm, 
         aes(ifelse(epm90 > median(epm90), "high EPM", "low EPM"), 
             y=log10(value+1), fill=ifelse(epm90 > median(epm90), "high.EPM", "log.EPM"))) +
  geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+ ylim(c(0.9,4.15))+
  gb2.theme("none") + ylab("somatic mutations (log10)") + scale_fill_manual(values=c("#f1a340", "#998ec3")) + xlab(NULL)
wilcox.test(value ~ epm90 > median(epm90), p.mut.full.bycov.num.epm, exact=F)

# figure 1c
mut.rtime=subset(p.mut.full.bycov.num.epm, stage==1)
library(survival)
plot(survfit(Surv(rtime/365) ~ (value > median(value)), mut.rtime), col=c("green","blue"),
     xlab="year",ylab="relapse-free probability", main="All AML (n=48)\np-value = 0.272", 
     xlim=range(mut.rtime$rtime/365))
legend("topright", col=c("green","blue"), c("Low MUT","High MUT"), lty=1, box.lwd=0, bty="n")
survdiff(Surv(rtime/365) ~ (value > median(value)), mut.rtime) # p= 0.272 
# figure3e
# increase/descreased proprotion of snvs
# coverage of all the loci with somatic mutations:
# vcf files generated by GATK 

cov.df1=cov.df
rownames(cov.df1)=gsub("_",":", rownames(cov.df))
readRD1=cov.df1[paste(misscov$seqnames,misscov$start, sep=":"),colnames(readRD)]

# coverage for reads count support REF and ALT alleles from gatk UnifiedGenotyper, 
# bedtools coverage results to count for the uncovered 191 snvs by gatk.
gatkcov=readVcf("paml.raw.variants", "hg19")
refRead=foreach(i=1:ncol(gatkcov), .combine=cbind) %do% sapply(geno(gatkcov)$AD[,i], function(x)x[1])
colnames(refRead)=gsub("B", "", colnames(gatkcov))
rownames(refRead)=paste(seqnames(rowRanges(gatkcov)), start(rowRanges(gatkcov)), sep=":")
refRead1=cov.df1[setdiff(paste(misscov$seqnames,misscov$start, sep=":"),rownames(refRead)),]
refRead2=rbind(refRead,refRead1)
altRead=foreach(i=1:ncol(gatkcov), .combine=cbind) %do% sapply(geno(gatkcov)$AD[,i], function(x)x[2])
colnames(altRead)=colnames(refRead)
rownames(altRead)=rownames(refRead)
altRead1=cov.df1[rownames(refRead1),]
altRead1[altRead1!=0]=0
altRead2=rbind(altRead,altRead1)
rFREQ2=altRead2/(altRead2+refRead2)

# run sciClone sciClone_1.0.7

snvs.for.sciclone=merge(p.snvs.full.df.bycov, unique(subset(p.mut.e3, type=="snvs")[,4:5]))
pids=unique(snvs.for.sciclone$pub.id)
library(sciClone)
registerDoMC(cores=20)
a=foreach(p1 = pids) %dopar%{
  Tsnvs=foreach(i = 1:2) %do% {
    # use all somatic snvs
    subset(snvs.for.sciclone, pub.id == p1 & stage == i)[, c(1,2,4,5,6)]
  }
  names(Tsnvs)=paste(p1, 1:2, sep="-")
  p2=subset(snvs.for.sciclone, pub.id==p1)[1,"pid"]
  acnvs=foreach(i = 1:2) %do% {
    sample.id=paste(p2,i,sep="-");
    cnvf=paste0(sample.id,"/", sample.id, "_mergedcnv.rda"); 
    load(cnvf)
    cnv
  }
  names(acnvs)=paste(p1, 1:2, sep="-")
  xids=paste(p1,1:2,sep="-");
  print(p1)
  #Tsnvs=gm.asnvs.exons[xids]
  # filling the SNVs gaps by making psudo zero vaf postion
  Tsnvsr=do.call(rbind, Tsnvs)
  Tsnvsu=unique(Tsnvsr[,1:2])
  rownames(Tsnvsu)=paste(Tsnvsu$seqnames, Tsnvsu$start,sep=":")
  Tsnvs2=lapply(1:2, function(i){
    X=Tsnvs[[i]]
    X$t.freq=100*X$t.freq
    x=setdiff(rownames(Tsnvsu), paste(X$seqnames, X$start, sep=":"))
    X2=Tsnvsu[x,]
    X2$t.ref.count=refRead2[x,paste(p2, i, sep="-")]
    X2$t.alt.count=altRead2[x,paste(p2, i, sep="-")]
    X2$t.freq=100* rFREQ2[x,paste(p2, i, sep="-")]
    colnames(X2)=colnames(X)
    rbind(X, X2)
  })
  for(i in c(6)){
    for(cluMethods in c("bmm")){ # 
      tryCatch({
        yids=paste(p2, 1, sep="-")
        system(paste("mkdir -p ", yids,"/sciclone/", sep=""))
        pfix=paste0(yids,"/sciclone/",yids,"_", p1,"_", cluMethods, "_sc2clo_Mar92016_max_allsomaticsnvs", i)
        if(!file.exists(paste0(pfix,".rda"))){
          sc2=sciClone(Tsnvs2, copyNumberCalls = acnvs[xids],sampleNames =xids,cnCallsAreLog2 = F,minimumDepth = 50, maximumClusters = i)
          sc.plot2d(sc2,outputFile = paste0(pfix, ".pdf"))
          save(sc2, file=paste0(pfix, ".rda"))
          writeClusterTable(sc2, paste0(pfix, ".txt"))
        }
      }, error=function(e){})
    }
  }
}

# summary sciClone
sc2=foreach(f = files1, .combine=rbind)%do%{
  x=read.table(f,header=T,sep="\t",stringsAsFactors = F)
  x1=x[which(!is.na(x$cluster) & x$cluster!=0),c(1:2,5,11,16)]
  x1$ID=gsub(".[12].vaf","", colnames(x1)[3])
  colnames(x1)[3:5]=c("Diagnosis","Relapse", "sciCluster")
  x1  
}
sc2.cl=merge(sc2, p.mut.bycov.cast.dsr[, c(1,3)], by.x="ID", by.y="pub.id")
library(extrafont)
font_import(pattern="[A/a]rial")
loadfonts(device="pdf")
for(ei in 1:3){
  p=ggplot(sc2.cl[which(sc2.cl$cluster==ei  & sc2.cl$ID %in% names(which(table(sc2.cl$ID) >=6))),], 
           aes(x=Diagnosis, y=Relapse, colour=as.factor(sciCluster), shape=as.factor(sciCluster)))+ 
    geom_point(size=0.75) + #xlim(c(0,100)) + ylim(c(0,100))+
    scale_shape(solid = FALSE) + gb2.theme() + scale_color_manual(values=cbbPalette)+facet_wrap(~ID, ncol=5)+
    xlab("Variant allele frequency at diagnosis") + ylab("Variant allele frequency at relapse")  +
    scale_x_continuous(breaks=c(0, 75)) + scale_y_continuous(breaks=c(0, 75)) 
  h=ifelse(ei %in% 1:2, 5,6)
  printtopdf(p, paste0("Suppfig7_", ei, ".pdf"),width=7,height=h)
}
# Figure 3e 3f
# increase/descreased proprotion of snvs
mut.merged.freq=foreach(p1 = pids, .combine=rbind) %do%{
  Tsnvs=foreach(i = 1:2) %do% {
    subset(snvs.for.sciclone, pub.id == p1 & stage == i)[, c(1,2,4,5,6)]
  }
  names(Tsnvs)=paste(p1, 1:2, sep="-")
  # filling the SNVs gaps by making psudo zero vaf postion
  Tsnvsr=do.call(rbind, Tsnvs)
  Tsnvsu=unique(Tsnvsr[,1:2])
  rownames(Tsnvsu)=paste(Tsnvsu$seqnames, Tsnvsu$start,sep=":")
  Tsnvs2=lapply(1:2, function(i){
    X=Tsnvs[[i]]
    X$t.freq=100*X$t.freq
    x=setdiff(rownames(Tsnvsu), paste(X$seqnames, X$start, sep=":"))
    X2=Tsnvsu[x,]
    X2$t.ref.count=refRead2[x,paste(p2, i, sep="-")]
    X2$t.alt.count=altRead2[x,paste(p2, i, sep="-")]
    X2$t.freq=rFREQ2[x,paste(p2, i, sep="-")]
    colnames(X2)=colnames(X)
    rbind(X, X2)
  })
  ecl=subset(p.mut.e3, pub.id==p1)[1,"cluster"]
  mmsnvs=data.frame(pub.id=p1, cluster=ecl, merge(Tsnvs2[[1]], Tsnvs2[[2]], by=c("seqnames","start")))
}
loci=rownames(cov.df1)
p.mut.bycov.cast.2=foreach(p1=unique(p.snvs.full.df.bycov$pub.id), .combine=rbind) %do% {
  p2=subset(p.snvs.full.df.bycov, pub.id==p1)[1,"pid"]
  keep.loci.pid=paste(p1, loci[rowSums(cov.df1[,paste(p2, c(1,2,"N"), sep="-")]>=10)==3], sep=":")
  cast.df=subset(p.mut.bycov.cast.1, paste(pub.id, seqnames, start, sep=":") %in% keep.loci.pid)
}

sc3=ddply(p.mut.bycov.cast.2, .(cluster, pub.id), summarize, inc=sum(Diagnosis < Relapse - 0.1)/length(Diagnosis), dec=sum(Relapse < Diagnosis - 0.1)/length(Diagnosis))
p=ggplot(sc3, aes(x=paste("cluster",cluster), y=inc, fill=as.factor(cluster))) + 
  geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+ ylim(c(0,0.75))+
  gb2.theme("none") + ylab("somatic mutations' proportion") + scale_fill_manual(values=cbPalette[5:7]) + xlab(NULL) 
printtopdf(p, "Figure3e_inc_48.pdf", width=4,heigh=4)
p=ggplot(sc3, aes(x=paste("cluster",cluster), y=dec, fill=as.factor(cluster))) + 
  geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+ ylim(c(0,0.75))+
  gb2.theme("none") + ylab("somatic mutations' proportion") + scale_fill_manual(values=cbPalette[5:7]) + xlab(NULL) 
printtopdf(p, "Figure3f_dec_48.pdf", width=4,heigh=4)
# supp table 5
ddply(sc3, .(cluster), summarize, med.inc=median(inc), med.dec=median(dec))
wilcox.test(dec~cluster,sc3[sc3$cluster!=2,]) 
wilcox.test(inc~cluster,sc3[sc3$cluster!=2,]) 

# Supp fig 6e
mut.rtime$mutcl=ifelse(mut.rtime$value>median(mut.rtime$value),"high MUT", "low MUT")
num.dat.d=ddply(mut.rtime, .(cluster, mutcl), summarize, value=length(value))
p=ggplot(num.dat.d, aes(x=mutcl, y=value, fill=as.factor(cluster)))+
  geom_bar(stat="identity", position="fill") +
  gb2.theme("right") + ylab("eloci cluster proportion") + scale_fill_manual(values=cbPalette[5:7]) + xlab(NULL)
printtopdf(p, "~/awork/projects/paml/data4reviewerNov30/figures/suppfig_6e.pdf", width=4,height=4)
chisq.test(matrix(num.dat.d$value, byrow=F, ncol=3)) 
# relapse
mut.rtimer=subset(p.mut.full.bycov.num.epm, stage==2)
mut.rtimer$mutcl=ifelse(mut.rtimer$value>median(mut.rtimer$value),"high MUT", "low MUT")
num.dat.r=ddply(mut.rtimer, .(cluster, mutcl), summarize, value=length(value))
chisq.test(matrix(num.dat.r$value, byrow=F, ncol=3)) # p-value = 0.05378
chisq.test(with(mut.rtimer, table(cluster, mutcl)))
p=ggplot(num.dat.r, aes(x=mutcl, y=value, fill=as.factor(cluster)))+
  geom_bar(stat="identity", position="fill") +
  gb2.theme("none") + ylab("eloci cluster proportion") + scale_fill_manual(values=cbPalette[5:7]) + xlab(NULL)
printtopdf(p, "suppfig_6e_forrelapse.pdf", width=4,height=4)

# Supp Figure 8
sc30=ddply(sc2.cl, .(ID, cluster, sciCluster), summarize, meanD=mean(Diagnosis), meanR=mean(Relapse))
sc31=ddply(sc30, .(cluster, ID), function(x){data.frame(D=sum(x$meanD>=5), R=sum(x$meanR>= 5))})
sc32=ddply(sc31, .(cluster, ID), function(x){
  if(x$D>x$R){type="decreasing"} else if (x$D< x$R){type="increasing"}else type="stable"})
sc33=ddply(sc32, .(cluster), function(x)data.frame(table(x$V1)))
p=ggplot(sc33, aes(x=cluster, y=Freq, fill=Var1)) + geom_bar(stat='identity', position="fill") +
  gb2.theme("right") + ylab("Proportion of patients")  + xlab("Cluster")
printtopdf(p, "suppfig_8_cloevol_inc_dec_48.pdf",width=5,height=4)
# percentage of patients with increased/decreased and stable clonality for manuscript.
sc34=ddply(sc33, .(Var1), summarize, num=sum(Freq))
100*sc34$num/sum(sc34$num)
chisq.test(with(sc32, table(cluster, V1)), simulate.p.value = T, B=1e6) 
