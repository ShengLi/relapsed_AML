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

