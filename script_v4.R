# Using s_matrix_all. This is the command used:
# ./bin/gibbsmod_data 100000 1000 10 10 2 s_matrix_all.dat results_040519

folder="test_2205"

pmu<-read.table(paste(folder,"/qmuk.dat",sep=""))
P<-read.table(paste(folder,"/P.dat",sep=""))
occup<-as.matrix(read.table(paste(folder,"/occupancy.dat",sep="")))
mem_traj<-read.table(paste(folder,"/membership_traj.dat",sep=""))
PMAX=4
pmu2<-pmu[P$V1==PMAX,]
pdf(paste(folder,"/prob_P.pdf",sep=""),5,5)
barplot(table(P$V1)/length(P$V1))
dev.off()

pm<-as.numeric(apply(as.matrix(pmu)[P$V1==PMAX,],2,mean,na.rm=T)[1:(4*PMAX)])
pdf(paste(folder,"/pmu.pdf",sep=""),8,6)
image(x=1:4,y=1:PMAX,z=matrix(pm,4,PMAX),col=grey.colors(100),axes=F,xlab="",ylab="")
dev.off()

pdf(paste(folder,"/patprob.pdf",sep=""))
layout(matrix(1:PMAX,PMAX))
par(mar=c(2,4,2,2))
for(k in 1:PMAX){
	boxplot(pmu2[,1:4+(k-1)*4],
		ylim=c(0,0.25),
		xlab="layers", ylab="Occupancy probability",
		names=c("II/III","IV","V","VI"), outline=F)
}
dev.off();

s<-read.table("s_matrix_all.dat")
b=s!=0
tab<-tabulate(apply(b,1,function(x) sum(c(1,2,4,8)[x])))

pdf(paste(folder,"/results.pdf",sep=""),12,5)
par(mar=c(5,5,3,3))
conf<-read.table(paste(folder,"/configurations.dat",sep=""))
ds<-t(apply(conf,1,function(x) tabulate(x,15)))
bp<-boxplot(ds,outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

pdf(paste(folder,"/results_comb.pdf",sep=""),5,5)
tab_comb=c(tab[1]+tab[2]+tab[3],
		   tab[4]+tab[8]+tab[12],
		   tab[5]+tab[6]+tab[7]+tab[9]+tab[10]+tab[11]+tab[13]+tab[14]+tab[15])
ds_comb<-t(apply(ds,1,function(x) c(x[1]+x[2]+x[3],x[4]+x[8]+x[12],x[5]+x[6]+x[7]+x[9]+x[10]+x[11]+x[13]+x[14]+x[15])))
bp_comb<-boxplot(ds_comb,outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab_comb,col="red",pch=7)
legend("topleft",legend=c("Data","Model"),pch=c(7,NA),lty=c(NA,1),col=c("red","black"),lwd=c(0,2))
dev.off()

posiz=expand.grid(1:4,1:PMAX)
pdf(paste(folder,"/map+values.pdf",sep=""),10,8)
image(x=1:4,y=1:PMAX,z=matrix(pm,4,PMAX),col=grey.colors(100),yaxt='n',xlab="",ylab="")
text(x=posiz$Var1,y=posiz$Var2,labels=pm)
dev.off()

plot_clonalsize <- function(burn=100){
	cop = apply(occup[burn:nrow(occup),],1, function(x) tabulate(colSums(matrix(as.numeric(x),4)),20))
    boxplot(t(cop))	
    points(1:20,tabulate(rowSums(s),nbins=20),col="red")

}

plot_layerocc <- function(){
	layers=c("II/III", "IV", "V", "VI")
	layout(matrix(1:4,2,2))
	par(mar=c(3,3,1,1))
	for(l in 1:4){
		cop = apply(occup,1, function(x) tabulate(matrix(as.numeric(x),4)[l,],20))
    	boxplot(t(cop),main=layers[l])	
    	points(1:20,tabulate(s[,l],nbins=20),col="red",pch=7)
	}
}

pdf(paste(folder,"/layer_occupancy.pdf",sep=""),10,8)
plot_layerocc()
dev.off()
dist_test<-function(n=100,write=F){
	Gen=matrix(NA,n,103)
	for(i in 1:n){
	   	for(k in 1:103){
	   		Gen[i,k] = rbinom(n=1,size=20,prob=0.10432247) + rbinom(n=1,size=20,prob=0.10781719) + rbinom(n=1,size=20,prob=0.05195175) + rbinom(n=1,size=20,prob=0.08644125)
		}
	}

	if(write) pdf(paste(folder,"/distribution.pdf",sep=""))   
	boxplot(t(apply(Gen,1,tabulate,nbins=20)))
    points(1:20,tabulate(rowSums(s),nbins=20),col="red")
	if(write) dev.off()
}

plot_s_ord <- function(ind){
	par(mar=c(4,4,2,2))
	order_membership=order(as.numeric(mem_traj[ind,]))
	image(x=1:4,y=1:nrow(s),z=t(as.matrix(s)[order_membership,]),col=grey.colors(100),axes=FALSE)
	axis(2,at=by(1:nrow(s),sort(as.numeric(mem_traj[ind,])),mean),labels=0:max(mem_traj[ind,]))
}

