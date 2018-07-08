pmu<-read.table("test/pmu.dat")
par(mar=c(1,1,1,1));np=5;layout(matrix(1:(4*np),nrow=4));for(i in 1:4) for(j in 1:np) hist(pmu[(j+np*(500:(nrow(pmu)/np-1))),i],breaks=seq(0,1,0.01),main="",xlab="",ylab="")


s<-read.table("s_matrix.dat")
b=s!=0
tab<-tabulate(apply(b,1,function(x) sum(c(1,2,4,8)[x])))

pdf("results.pdf",12,5)
par(mar=c(5,5,3,3))
conf<-read.table("test/configurations.dat")
ds<-t(apply(conf,1,function(x) tabulate(x)))
boxplot(ds[500:1000,],outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()




