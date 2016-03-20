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

load("b138.aug_rtimeupdate_137patient.rda")
AML147DR=read.table("DR_eloci_147_Aug18.txt",header=F,sep="\t",stringsAsFactors=F)
AML147DN=read.table("AML147_DN.txt",header=F,sep="\t",stringsAsFactors=F)
AML147RN=read.table("AML147_RN.txt",header=F,sep="\t",stringsAsFactors=F)

AML144DN9=read.table("AML144_DN_9N_Sep6.txt",header=F,sep="\t",stringsAsFactors=F)
AML144RN9=read.table("AML144_RN_9N_Sep6.txt",header=F,sep="\t",stringsAsFactors=F)

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
printtopdf(p, "Figure2_epm90_epicluster.pdf")
pdf("Figure2_epm90_epicluster_pval.pdf")
pval.by.stage=ddply(dat2, .(stage), function(X){data.frame(not.cl=1:3, pval=sapply(1:3, function(i) wilcox.test(epm90~cl, X[gsub("Cluster ","",X$cl)!=i,])$p.value))})
pval.by.cl=ddply(dat2, .(cl), function(X){data.frame(pval=wilcox.test(epm90~stage, X)$p.value)})
textplot(pval.by.stage)
textplot(pval.by.cl)
dev.off()
DN.epm90.aug=dat2$epm90[dat2$stage==1]
names(DN.epm90.aug)=dat2$sample[dat2$stage==1]
b138.aug$DN.epm90.aug=DN.epm90.aug[b138.aug$sample]

b138.aug$cl.aug=cl[b138.aug$sample]
summary(coxph(Surv(relapsetime.aug)~ cl.aug, b138.aug)) 
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=1),])  
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=2),]) 
survdiff(Surv(relapsetime.aug)~ cl.aug, b138.aug[which(b138.aug$cl.aug !=3),])  

# cluster and age/gender
anova(lm(dage~as.factor(cl.aug), b138.aug)) 

p=ggplot(b138.aug, aes(x=factor(cl.aug), y=dage, fill=factor(cl.aug))) + geom_violin() + stat_summary(fun.y=median.quartile,geom='point')+
  gb2.theme("none") + ylab("Age") + scale_fill_manual(values=cbPalette[5:7]) + xlab("Cluster") + ggtitle("All AML (n=138)\np-value=0.4413") 
printtopdf(p, pdf="aug2015/epm_vs_age.pdf",width=4,height=4)

anova(lm(sex~as.factor(cl.aug), b138.aug)) 
chisq.test(table(b138.aug$sex, b138.aug$cl.aug))

dat.aug=ddply(b138.aug, .(cl.aug), function(X)data.frame(table(X$sex)))
p=ggplot(dat.aug, aes(x=Var1, y=Freq,fill=factor(cl.aug))) + 
  geom_bar(stat="identity", position="fill") + 
  gb2.theme("none") + ylab("Proportion of patients\nfrom three clusters") + scale_fill_manual(values=cbPalette[5:7]) + xlab("Gender") + ggtitle("All AML (n=138)\np-value=0.01048") 
printtopdf(p, pdf="epm_vs_sex.pdf",width=4,height=4)


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
printtopdf(p,"figure2d.pdf", width=5,height=2.5)

p=ggplot( dat2[dat2$stage==2,], aes(x=as.numeric(ID), y=log10(epm90+1), fill=as.factor(stage))) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[3]) + ylim(c(0,5.1))
printtopdf(p,"figure2b.pdf", width=5,height=2.5)

p=ggplot( dat2[dat2$stage==1,], aes(x=as.numeric(ID), y=log10(epm90+1), fill=as.factor(stage))) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[2]) + ylim(c(0,5.1))
printtopdf(p,"~/awork/projects/paml/figures/Figure2/figure2a.pdf", width=5,height=2.5)

wilcox.test(dat2[dat2$stage==1,"epm90"],dat2[dat2$stage==2,"epm90"]) #p-value = 0.00463

# Figure 3 # 
# Eloci shared
print(load("fig3_14N.shared.aug.rda"))
greens=c("#a1d99b","#41ab5d", "#006d2c")
dat4$variable=gsub("\\.", " ", dat4$variable)
p=ggplot(dat4, aes(ID, epiallele.shift,fill=variable))+ 
  geom_bar(stat="identity", position="fill") + scale_fill_manual(values=greens) +
  facet_grid(~cluster,scale="free",labeller=label_both, space="free_x") + gb2.theme("bottom") +
  scale_x_continuous(breaks=c(1, 138))
printtopdf(p, "fig3.shared_14N_3_epi_green2_oct.pdf", width=8,height=5)

# Figure 4 #
library(pheatmap);
library(doMC);

sumt=read.fast("T1_entropy.txt",header=T,sep="\t")
pheatmap(sumt[which(rowMeans(sumt[,-1:-3])< -30), -1:-3])

files=dir(".", pattern="_entropy.txt", recursive = T, full.names = T)
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
pdf("entropy_heatmap.pdf",width=4,height=4)
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
