## function section
# get eff from vcf file | general function
eff.from.vcf=function(vcf){
  eff=info(vcf)$EFF
  idx=which(sapply(eff, function(x)any(grepl("HIGH|MODERATE",x))) )
  eff2=foreach(i=eff, .combine=rbind) %do% {
    i1=grep("HIGH|MODERATE", i, value=T)
    biotype=paste(unique(splitn(i1,"\\(",1)), collapse=",")
    genes=paste(unique(splitn(i1,"\\|",6)), collapse=",")
    aa=paste(unique(splitn(i1,"\\|",4)), collapse=",")
    
    data.frame(biotype = biotype, genes=genes, aa=aa)
  }
}
# read varscan vcf file
## filter by: somatic pvalue < 0.05 & somatic status  
read.varscan.vcf=function(varscan.vcffile, exome.gr) {
  vvcf=readVcf(varscan.vcffile,genome = "hg19")
  o=findOverlaps(rowRanges(vvcf), exome.gr)
  idx=unique(o@queryHits)
  vvcf = vvcf[idx]
  idx=which(info(vvcf)$SPV < 0.05 & info(vvcf)$SS == 2)
  vvcf = vvcf[idx]
  
  t="TUMOR"
  n="NORMAL"
  t.depth=geno(vvcf)$DP[,t]
  t.ref.count=geno(vvcf)$RD[,t]
  t.alt.count=geno(vvcf)$AD[,t]
  t.freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ)[,t])/100
  
  n.depth=geno(vvcf)$DP[,n]
  n.alt.count=geno(vvcf)$AD[,n]
  n.ref.count=geno(vvcf)$RD[,n]
  n.freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ)[,n])/100
  v.info=as.data.frame(rowRanges(vvcf))
  ALT = unlist(sapply(v.info[,8], function(x)paste(as.character(x), collapse=",")))
  vsnvs=data.frame(v.info[,1:2],t.depth, t.ref.count, t.alt.count,t.freq, n.depth, 
                   n.ref.count, n.alt.count, n.freq, REF=v.info[,7], ALT)
  vgenes=data.frame(vsnvs, eff.from.vcf(vvcf))
}
# read varscan indel vcf file
## filter by: somatic pvalue < 0.05 & somatic status  
read.vindels.vcf=function(vindels.vcffile, exome.gr) {
  indel=readVcf(vindels.vcffile, genome="hg19")
  o=findOverlaps(rowRanges(indel), exome.gr)
  idx=unique(o@queryHits)
  vvcf = indel[idx]
  idx=which(info(vvcf)$SPV < 0.05 & info(vvcf)$SS == 2)
  vvcf = vvcf[idx]
  
  t="TUMOR"
  n="NORMAL"
  t.depth=geno(vvcf)$DP[,t]
  t.ref.count=geno(vvcf)$RD[,t]
  t.alt.count=geno(vvcf)$AD[,t]
  t.freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ)[,t])/100
  
  n.depth=geno(vvcf)$DP[,n]
  n.alt.count=geno(vvcf)$AD[,n]
  n.ref.count=geno(vvcf)$RD[,n]
  n.freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ)[,n])/100
  v.info=as.data.frame(rowRanges(vvcf))
  ALT = unlist(sapply(v.info[,8], function(x)paste(as.character(x), collapse=",")))
  vsnvs=data.frame(v.info[,1:2],t.depth, t.ref.count, t.alt.count,t.freq, n.depth, 
                   n.ref.count, n.alt.count, n.freq, REF=v.info[,7], ALT)
  vgenes=data.frame(vsnvs, eff.from.vcf(vvcf))
}
# read mutect vcf file
read.mutect.vcf=function(mutect.vcffile, exome.gr) {
  mvcf=readVcf(mutect.vcffile,"hg19")
  o=findOverlaps(rowRanges(mvcf), exome.gr)
  idx=unique(o@queryHits)
  mvcf = mvcf[idx]
  
  # mutect
  pid2=grep("N", colnames(geno(mvcf)$AD), value=T, invert=T)
  nid2=grep("N", colnames(geno(mvcf)$AD), value=T, invert=F)
  t.alt.count=round(unlist(geno(mvcf)$DP[,pid2])* unlist(geno(mvcf)$FA[,pid2]))
  t.depth=unlist(geno(mvcf)$DP[,pid2])
  t.ref.count=t.depth-t.alt.count
  t.freq=unlist(geno(mvcf)$FA[,pid2]) 
  n.alt.count=round(unlist(geno(mvcf)$DP[,nid2])* unlist(geno(mvcf)$FA[,nid2]))
  n.depth=unlist(geno(mvcf)$DP[,nid2])
  n.ref.count=n.depth-n.alt.count
  n.freq=unlist(geno(mvcf)$FA[,nid2]) 
  m.info=as.data.frame(rowRanges(mvcf))
  ALT = unlist(sapply(m.info[,8], function(x)paste(as.character(x), collapse=",")))
  msnvs=data.frame(m.info[,1:2],t.depth, t.ref.count, t.alt.count, t.freq, 
                   n.depth, n.ref.count, n.alt.count, n.freq, REF=m.info[,7], ALT)
  mgenes=data.frame(msnvs, eff.from.vcf(mvcf))
}
# read somaticsniper vcf file
## filter by: average mapping quality > 40 & somatic score > 40 
read.ss.vcf=function(ss.vcffile,exome.gr){
  svcf=readVcf(ss.vcffile,"hg19")
  # filter by exome
  o=findOverlaps(rowRanges(svcf), exome.gr)
  idx=unique(o@queryHits)
  svcf = svcf[idx]
  # filter by average mapping quality > 40 & somatic score > 40
  idx=which(geno(svcf)$TMQ[,"TUMOR"] > 40 & 
              geno(svcf)$SSC[,"TUMOR"] > 40)
  svcf = svcf[idx]
  if (length(svcf) > 0) {
    # somatic sniper
    pid2=grep("N", colnames(geno(svcf)$AD), value=T, invert=T)
    nid2=grep("N", colnames(geno(svcf)$AD), value=T, invert=F)
    t.alt.count=rowSums(geno(svcf)$DP4[,"TUMOR",3:4])
    t.depth=rowSums(geno(svcf)$DP4[,"TUMOR",])
    t.ref.count=rowSums(geno(svcf)$DP4[,"TUMOR",1:2])
    t.freq=t.alt.count/t.depth
    
    n.alt.count=rowSums(geno(svcf)$DP4[,"NORMAL",3:4])
    n.depth=rowSums(geno(svcf)$DP4[,"NORMAL",])
    n.ref.count=rowSums(geno(svcf)$DP4[,"NORMAL",1:2])
    n.freq=n.alt.count/n.depth
    s.info=as.data.frame(rowRanges(svcf))
    ALT = unlist(sapply(s.info[,8], function(x)paste(as.character(x), collapse=",")))
    ssnvs=data.frame(s.info[,1:2],t.depth, t.ref.count, t.alt.count,t.freq, n.depth, 
                     n.ref.count, n.alt.count, n.freq, REF=s.info[,7],ALT)
    mgenes=data.frame(ssnvs, eff.from.vcf(svcf))
  } else {
    mgenes = NULL
  }
}
# get snvs id using chr and start
get.snvs.id=function(snvs){
  paste(snvs$seqnames, snvs$start, as.character(snvs$REF), as.character(snvs$ALT), sep="_")
}
# filter snvs to only have snvs with genes 
snvs.w.genes=function(snvs){
  snvs[snvs$biotype!="",]
}
# merge three tools of snvs 
## filter: 1. allele frequency of normal < 0.05; 2. read coverage (normal and tumor) >=10
merge.snv=function(pid, exome.gr){
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
  
  v.id=get.snvs.id(snvs.w.genes(varscan.snvs.flt))
  m.id=get.snvs.id(snvs.w.genes(mutect.snvs.flt))
  s.id=get.snvs.id(snvs.w.genes(ss.snvs.flt))
  id.freq=table(c(s.id,m.id,v.id))
  somaticids=names(which(id.freq>1))

  # merge varscan + mutect
  merged.vcf=rbind(snvs.w.genes(mutect.snvs.flt), 
                   snvs.w.genes(varscan.snvs.flt)[which(v.id%in% setdiff( v.id, m.id)),],
                   snvs.w.genes(ss.snvs.flt)[which(s.id%in% setdiff( s.id, unique(c(m.id, v.id)))),])
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
# get the genomic ranges from data.frame
getgr=function(x){
  colnames(x)=c("chr","start","end")
  x$start=as.numeric(as.character(x$start))
  x$end=as.numeric(as.character(x$end))
  gr=GRanges(Rle(x$chr),IRanges(x$start, x$end))
  gr
}
