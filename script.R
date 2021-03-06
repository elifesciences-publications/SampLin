pmu<-read.table("test3/pmuk.dat")
par(mar=c(1,1,1,1));
np=10;
layout(matrix(1:(4*np),nrow=4));
for(i in 1:4){
	for(j in 1:np){
		hist(pmu[(j+np*(500:(nrow(pmu)/np-1))),i],breaks=seq(0,1,0.01),main="",xlab="",ylab="")
	}
}


s<-read.table("s_matrix.dat")
b=s!=0
tab<-tabulate(apply(b,1,function(x) sum(c(1,2,4,8)[x])))

pdf("results.pdf",12,5)
par(mar=c(5,5,3,3))
conf<-read.table("test3/configurations.dat")
ds<-t(apply(conf,1,function(x) tabulate(x)))
boxplot(ds[500:1000,],outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

pdf("results_comb.pdf",5,5)
tab_comb=c(tab[1]+tab[2]+tab[3],tab[4]+tab[8]+tab[12],tab[5]+tab[6]+tab[7]+tab[9]+tab[10]+tab[11]+tab[13]+tab[14]+tab[15])
ds_comb<-t(apply(ds,1,function(x) c(x[1]+x[2]+x[3],x[4]+x[8]+x[12],x[5]+x[6]+x[7]+x[9]+x[10]+x[11]+x[13]+x[14]+x[15])))
boxplot(ds_comb[500:1000,],outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab_comb,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

posiz=expand.grid(1:4,1:3)
pm<-as.numeric(apply(as.matrix(t)[P$V1==3,],2,mean)[1:12])
pdf("map+values.pdf",5,4)
image(x=1:4,y=1:3,z=matrix(pm,4,3),col=grey.colors(100),yaxt='n',xlab="",ylab="")
text(x=posiz$Var1,y=posiz$Var2,labels=pm)
dev.off()

