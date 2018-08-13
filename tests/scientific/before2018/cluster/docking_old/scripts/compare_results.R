data1 = read.table("benchmark_output.txt", header = T, sep =";")
data2 = read.table("benchmark_output_old.txt", header = T, sep = ";")

ntargets = length(data1$filename)
for (i in 1:ntargets){
	if (data1$Nrms5[i] == data2$Nrms5[i]){ 			#if values equal, then passed				
		print(paste(data1$filename[i]," PASSED"))
		} else if (data1$Nrms5[i]*data2$Nrms5[i] == 0){	#if values are unequal, and one of them is 0, then failed	
			print(paste(data1$filename[i]," FAILED"))
			} else {
				Nrms5_ratio = data1$Nrms5[i]/data2$Nrms5[i]	#otherwise pass/fail depends on ratio
				if (Nrms5_ratio > 0.5 | Nrms5_ratio < 2.0){
					print(paste(data1$filename[i]," PASSED"))
					} else {
						print(paste(data1$filename[i]," FAILED"))
						}
				}
			
		
		} 
