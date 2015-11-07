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
