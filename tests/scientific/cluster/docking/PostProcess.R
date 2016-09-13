scorelist = read.table("pdblist.txt", header = T)
ntarget = length(scorelist$filename)

mean_n5 = array(0,c(ntarget))
sd_n5 = array(0,c(ntarget))
p_value = array(0,c(ntarget))
mean_low_rmsd=array(100,c(ntarget))

for(a in 1:ntarget){
	print(scorelist$filename[a])
	d1 = read.table(paste(scorelist$filename[a]), header = T, fill=T,skip=1)

	n_cycles=5000
	n_decoy=1000
	nstruct=length(d1$rms)
	
	N4=array(0,c(n_cycles))
	low_rmsd = array(100.0,c(n_cycles))
	
	for (i in 1:n_cycles){
		rand_decoy=sample(1:nstruct,n_decoy,replace=T)
		Irmsd=d1$Irms[rand_decoy]
		Isc=d1$I_sc[rand_decoy]
		
		sorted = order(Isc)
		low_rmsd[i] = min(Irmsd[sorted[1:5]]) 
		for(k in 1:5){
			if (Irmsd[sorted[k]] <= 4.0){
				N4[i]=N4[i]+1
				}
			}
		}
	
	bins=c(0,1,2,3,4,5)
	h1=hist(N4, breaks=bins)
	mean_n5[a] = mean(N4)
	sd_n5[a] = sd(N4)
	mean_low_rmsd[a] = mean(low_rmsd)
	p_value[a] = (h1$counts[1]+h1$counts[2]+h1$counts[3])/sum(h1$counts)
	}	

output_table = data.frame(scorelist$filename,mean_n5, sd_n5, mean_low_rmsd, p_value)
write.table(output_table, file= "output/results.txt", quote = FALSE, sep = ";", row.names = FALSE)
	
postscript(file="funnels.ps",paper="letter", horizontal=T)
#layout(rbind(c(1, 2, 3), c(0, 3, 0))) 
par(mfrow=c(3,4))
par(cex.axis=1.5)
par(cex.lab=1.5)
par(cex.main=1.5)

myplot = function(x,y,...){

x.all = x
y.all = y
y.cutoff=y.all[order(y.all)[round(length(y.all)*0.80)]]
subset = y.all<=y.cutoff
x.all = x.all[subset]
y.all = y.all[subset]

plot( x.all,y.all,...)

}

dataset=read.table("pdblist.txt",header=T)

for(i in 1:dim(dataset)[1]){ 
	data=read.table(paste(dataset$filename[i]),header=T,skip=1)
	myplot( data$rms, data$I_sc,xlim=c(0, 30), main=paste(dataset$filename[i]), xlab="RMSD", ylab="Interface Score", cex.main=1.5, cex.axis=1.5, cex.lab=1.5, pch=20, col="black")
	myplot( data$Irms, data$I_sc,xlim=c(0, 15), main=paste(dataset$filename[i]), xlab="I_RMSD", ylab="Interface Score", cex.main=1.5, cex.axis=1.5, cex.lab=1.5, pch=20, col="black")
}