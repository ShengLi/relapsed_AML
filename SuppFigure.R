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
