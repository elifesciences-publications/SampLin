# Using s_matrix_all. This is the command used:
# ./bin/gibbsmod_data 100000 1000 10 10 2 s_matrix_all.dat results_040519

folder="results_150619"

pmu<-read.table(paste(folder,"/pmuk.dat",sep=""))
P<-read.table(paste(folder,"/P.dat",sep=""))
occup<-as.matrix(read.table(paste(folder,"/occupancy.dat",sep="")))
mem_traj<-read.table(paste(folder,"/membership_traj.dat",sep=""))
PMAX=1
pmu2<-pmu[P$V1==PMAX,]

mem_traj_tmp<-mem_traj[P$V1==PMAX,]
mem <- rep(NA,ncol(mem_traj))
mem[1] = which.max(tabulate(mem_traj_tmp[,1]+1,10))-1

for(i in 2:ncol(mem_traj)){
	mem_traj_tmp <- mem_traj_tmp[mem_traj_tmp[,i-1]==mem[i-1],]
	mem[i] = which.max(tabulate(mem_traj_tmp[,i]+1,10))-1
}

pdf(paste(folder,"/prob_P.pdf",sep=""),5,5)
barplot(table(P$V1)/length(P$V1))
dev.off()

if(nrow(mem_traj_tmp)>1){
	pm<-as.numeric(apply(as.matrix(pmu)[as.numeric(rownames(mem_traj_tmp)),],2,mean,na.rm=T)[1:(4*PMAX)])
} else {
	pm<-as.numeric(as.matrix(pmu)[as.numeric(rownames(mem_traj_tmp)),][1:(4*PMAX)])
}

pdf(paste(folder,"/pmu.pdf",sep=""),8,6)
image(x=1:4,y=1:PMAX,z=matrix(pm,4,PMAX),col=grey.colors(100),axes=F,xlab="",ylab="")
dev.off()

pdf(paste(folder,"/patprob.pdf",sep=""))
layout(matrix(1:PMAX,PMAX))
par(mar=c(2,4,2,2))
for(k in 1:PMAX){
	boxplot(pmu[as.numeric(rownames(mem_traj_tmp)),1:4+(k-1)*4],
		ylim=c(0,.6),
		xlab="layers", ylab="Occupancy probability",
		names=c("II/III","IV","V","VI"), outline=F)
}
dev.off();

s<-read.table("s_matrix.dat")
stop=rep(0,nrow(s))
for(i in 1:nrow(s)){
   k=0; while(s[i,k+1]==0) k=k+1;
   stop[i]=k;
}


b=s!=0
tab<-tabulate(apply(b,1,function(x) sum(c(1,2,4,8)[x])))

pdf(paste(folder,"/results.pdf",sep=""),12,5)
par(mar=c(5,5,3,3))
clon_sz=apply(occup,1,function(x) colSums(matrix(x,4)))
conf<-read.table(paste(folder,"/configurations.dat",sep=""))
conf_cs<-cbind(as.matrix(conf),t(clon_sz))
ds<-t(apply(conf_cs,1,function(x) tabulate(x[1:103][x[104:206]<=12],15)))
bp<-boxplot(ds,outline=FALSE,cex.axis=2,cex.lab=2,xlab="Configurations",ylab="Frequencies")
points(tab,col="red",pch=7)
ds_mean=apply(ds,2,mean)
ds_sd=apply(ds,2,sd)
z_conf=mean(apply(ds,1,function(x) sum((x-tab)^2/ds_sd^2)/15))
p_conf=pchisq(z_conf,1,lower.tail=F)

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
ds_comb_mean=apply(ds_comb,2,mean)
ds_comb_sd=apply(ds_comb,2,sd)
z_comb=mean(apply(ds_comb,1,function(x) sum((x-tab_comb)^2/ds_comb_sd^2)/3))
p_comb=pchisq(z_comb,1,lower.tail=F)
dev.off()

# calculate the z-score of the translaminar ones
trans_laminar=c(5,6,7,9,10,11,13,14,15)
tab_tl=tab[trans_laminar]
ds_tl=apply(ds[,trans_laminar],2,mean)
ds_tl_sd=apply(ds[,trans_laminar],2,sd)
z_tl=mean(apply(ds[,trans_laminar],1,function(x) sum((x-tab_tl)^2/ds_tl_sd^2)/length(trans_laminar)))
p_tl=pchisq(z_tl,1,lower.tail=F)


posiz=expand.grid(1:4,1:PMAX)
pdf(paste(folder,"/map+values.pdf",sep=""),10,8)
image(x=1:4,y=1:PMAX,z=matrix(pm,4,PMAX),col=grey.colors(100),yaxt='n',xlab="",ylab="")
text(x=posiz$Var1,y=posiz$Var2,labels=pm)
dev.off()

plot_clonalsize <- function(burn=100){
	cop = apply(occup[burn:nrow(occup),],1, function(x) tabulate(colSums(matrix(as.numeric(x),4)),20)[3:12])
    cop_data=tabulate(rowSums(s),nbins=20)[3:12]
	cop_var=apply(cop,1,var)
    val=list(means=apply(cop,1,mean),stddevs=apply(cop,1,sd))
	b=barplot(t(cbind(cop_data,val$means)),beside=T,ylim=c(0,18),names.arg=3:12)
	arrows(x0=b[2,],x1=b[2,],y0=val$means-val$stddevs,y1=val$means+val$stddevs,angle=90,code=3,length=0)
	z1=apply(cop,2,function(x) sum((x-cop_data)^2/cop_var/9))
	return(list(z=mean(z1),pval=pchisq(mean(z1),1,lower.tail=F)))

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

	clonal_size=matrix(NA,n,103)
	Gen2=array(NA,dim=c(n,103,4))
	pm2=matrix(pm,4,PMAX)
	for(i in 1:n){
	   	for(k in 1:103){
			for(l in 1:4) Gen2[i,k,l] = rbinom(n=1,size=20,prob=pm2[l,mem[k]+1])
			clonal_size[i,k] = sum(Gen2[i,k,(1:4) > stop[k]])
		}
	}

	if(write) pdf(paste(folder,"/distribution.pdf",sep=""))   
	boxplot(t(apply(clonal_size,1,tabulate,nbins=20)))
    points(1:20,tabulate(rowSums(s),nbins=20),col="red")
	if(write) dev.off()
}

compare_correlation <- function(){
	cc=rep(NA,nrow(occup))
	for(i in 1:nrow(occup)){
		tmp=apply(matrix(occup[i,],4),2,function(x) c(x[1]+x[2],x[3]+x[4]))
		tmp=tmp[,colSums(tmp)<=12 & colSums(tmp)>=3]
		cc[i]=cor(t(tmp),method="spearman")[2,1]
	}

	uu=cor(x=s[,1]+s[,2],y=s[,3]+s[,4],method="spearman")
	plot(x=cc,y=1:length(cc),pch=19,cex=0.2)
	abline(v=uu,lty=2,col="green",lwd=2)
	z=mean((cc-uu)^2/var(cc))
	return(list(zscore=mean((cc-uu)^2/var(cc)),pval=pchisq(z,1,lower.tail=F)))
}


plot_s_ord <- function(ind=NA){

	par(mar=c(4,4,2,2))
	if(is.na(ind)) {
		order_membership=order(mem)
		mem2=mem
	} else {
		order_membership=order(as.numeric(mem_traj[ind,]))
		mem2=as.numeric(mem_traj[ind,])
	}

	image(x=1:4,y=1:nrow(s),z=t(as.matrix(s)[order_membership,]),col=grey.colors(100),axes=FALSE)
	axis(2,at=by(1:nrow(s),sort(as.numeric(mem2)),mean),labels=0:max(mem2))
}

cortical_layer_occupancy<- function(){
	# occupancy from the data
	lo_data = colSums(s)/sum(s)
	lo=apply(occup,1,function(x) rowSums(matrix(x,4))/sum(x))
	lo_var=apply(lo,1,var)
    val=list(means=apply(lo,1,mean),stddevs=apply(lo,1,sd))
	z1=apply(lo,2,function(x) sum((x-lo_data)^2/lo_var/4))
	b=barplot(t(cbind(val$means,lo_data)),beside=T,ylim=c(0,max(val$means+val$stddevs)))
	arrows(x0=b[2,],x1=b[2,],y0=val$means-val$stddevs,y1=val$means+val$stddevs,angle=90,code=3,length=0)
	return(list(z=mean(z1),pval=pchisq(mean(z1),1,lower.tail=F)))
}


