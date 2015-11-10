# Figure 2 # 

brary(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# barplot 1
# define a data.frame "dn.epm"
# Column 1
# Sample
# Column 2
# EPM: Calculate the EPM value from methclone (number of eloci per million loci sequenced) between diagnosis vs. normal bone mrrow
# Average of all 14 NBM

# get the order by Diagnosis vs NBM
od.sample=dn.epm$Sample[order(dn.epm$EPM)]
dn.epm$ID=factor(dn.epm$Sample, levels=od.sample)

# plot for D vs N
ggplot( dn.epm, aes(x=as.numeric(ID), y=log10(EPM+1), fill="D")) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[2]) + ylim(c(0,5.1))


# barplot 2
# define a data.frame "rn.epm"
# Column 1
# Sample 
# Column 2
# EPM: Calculate the EPM value from methclone (number of eloci per million loci sequenced) between relapse vs. normal bone mrrow
# Average of all 14 NBM

# get the order by Diagnosis vs NBM
od.sample=dn.epm$Sample[order(dn.epm$EPM)]
rn.epm$ID=factor(rn.epm$Sample, levels=od.sample)

# plot for R vs N
ggplot( rn.epm, aes(x=as.numeric(ID), y=log10(EPM+1), fill="R")) +
  geom_bar(stat="identity") + xlab("sample index") +scale_x_continuous(breaks=c(1, 50,100, 138)) +
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[3]) + ylim(c(0,5.1))

# barplot 3 
# define a data.frame "rd.epm"
# Column 1
# Sample
# Column 2
# EPM: Calculate the EPM value from methclone (number of eloci per million loci sequenced) between relapse vs. diagnosis

# get the order by Diagnosis vs NBM
rd.epm$ID=factor(rd.epm$Sample, levels=od.sample)

# plot for R vs D
ggplot( rd.epm, aes(x=as.numeric(ID), y=log10(EPM+1), fill="RD")) + 
  geom_bar(stat="identity") + xlab("sample index")  +scale_x_continuous(breaks=c(1, 50,100, 138))+
  gb2.theme("none") + ylab("EPM (log10)") + scale_fill_manual(values=cbPalette[4]) + ylim(c(0,5.1))


