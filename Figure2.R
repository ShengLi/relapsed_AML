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
