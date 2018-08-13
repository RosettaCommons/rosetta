#file_list = read.table("filelist", header = T)
file_list = read.table("filelist")
#file_number = length(file_list$filename)
file_number = length(file_list[,1])
#filename = file_list$filename 
filename = file_list[ ,1] 
print(filename)
#SCIENTIFIC BENCHMARK METRICS 

RMS_cutoff = 5.0

Nrms5 = array(0, c(file_number)) #used internally
Nrms10 = array(0, c(file_number)) #used internally

AvgErms5 = array(0, c(file_number))  #average energy of structures within 5 Lrmsd
AvgErms10 = array(0, c(file_number)) #average energy of structures within 10 Lrmsd
rmsdist = array(0, c(file_number, 12)) #histogram of rms distribution for the decoy set
top_rms5 = array(0, c(file_number))  #number among top 5 scorers that meets RMS_cutoff
top_rms10 = array(0, c(file_number)) #number among top 10 scorers that meets RMS_cutoff
ndecoys = array(0, c(file_number))

for (j in 1:file_number){
	filepath = paste("./",filename[j],"/XT",filename[j],".fasc", sep="")
	print(filepath)
	data1 = read.table(filepath, header=T, skip=1)
	ndecoys[j] = length(data1$rms)
	sorted = order(data1$total_score)
	for (i in 1:ndecoys[j]){
		if(data1$rms[i] <= 5.0){
			Nrms5[j] = Nrms5[j]+1
			AvgErms5[j] = AvgErms5[j]+data1$total_score[i]
			}
		if(data1$rms[i] <= 10.0){
			Nrms10[j] = Nrms10[j]+1
			}
		if(data1$rms[i] > 10.0){
			AvgErms10[j] = AvgErms10[j]+data1$total_score[i]
			}
	#rms distribution calculation
		bin_size = 5
		bin_num = ((data1$rms[i]-data1$rms[i]%%bin_size)/bin_size)+1
		if (bin_num > 12){
			bin_num = 12
			}
		rmsdist[j,bin_num] = rmsdist[j,bin_num]+1

		}
	#top scorer results
	for (i in 1:5){
		if(data1$rms[sorted[i]] <= 5.0){ 
			top_rms5[j] = top_rms5[j]+1
			}
		}
	for (i in 1:10){
		if(data1$rms[sorted[i]] <= 5.0){ 
			top_rms10[j] = top_rms10[j]+1
			}
		}

# convert to percentages
	top_rms5[j] = top_rms5[j]/5
	top_rms10[j] = top_rms10[j]/10
	AvgErms5[j] = AvgErms5[j]/Nrms5[j]
	AvgErms10[j] = AvgErms10[j]/(ndecoys[j]-Nrms10[j])
	}
BE = AvgErms5 - AvgErms10
output_table = data.frame(filename, Nrms5, Nrms10, BE, top_rms5, top_rms10, ndecoys, AvgErms5 )
write.table(output_table, file = "benchmark_output.txt", quote = FALSE, sep = ";", row.names=FALSE)
