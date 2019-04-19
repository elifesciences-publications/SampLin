# ./bin/gibbs_data 200000 0 100 1 1 s_matrix.dat test3


pmu<-read.table("test3/pmuk.dat")
P<-read.table("test3/P.dat")
par(mar=c(1,1,1,1));
np=10;

pdf("prob_P.pdf",5,5)
barplot(table(P$V1)/2000)
dev.off()

pm<-as.numeric(apply(as.matrix(pmu)[P$V1==3,],2,mean)[1:12])
pdf("pmu.pdf",8,6)
image(x=1:4,y=1:3,z=matrix(pm,4,3)[,c(2,1,3)],col=grey.colors(100),axes=F,xlab="",ylab="")
dev.off()

s<-read.table("s_matrix.dat")
b=s!=0
tab<-tabulate(apply(b,1,function(x) sum(c(1,2,4,8)[x])))

pdf("results.pdf",12,5)
par(mar=c(5,5,3,3))
conf<-read.table("test3/configurations.dat")
ds<-t(apply(conf,1,function(x) tabulate(x)))
bp<-boxplot(ds[500:2000,],outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

pdf("results_comb.pdf",5,5)
tab_comb=c(tab[1]+tab[2]+tab[3],tab[4]+tab[8]+tab[12],tab[5]+tab[6]+tab[7]+tab[9]+tab[10]+tab[11]+tab[13]+tab[14]+tab[15])
ds_comb<-t(apply(ds,1,function(x) c(x[1]+x[2]+x[3],x[4]+x[8]+x[12],x[5]+x[6]+x[7]+x[9]+x[10]+x[11]+x[13]+x[14]+x[15])))
bp_comb<-boxplot(ds_comb[500:1000,],outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab_comb,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

posiz=expand.grid(1:4,1:3)
pdf("map+values.pdf",10,8)
image(x=1:4,y=1:3,z=matrix(pm,4,3),col=grey.colors(100),yaxt='n',xlab="",ylab="")
text(x=posiz$Var1,y=posiz$Var2,labels=pm)
dev.off()

