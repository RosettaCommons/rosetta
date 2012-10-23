import os
import sys
import re

qsubtemp = "path"
tempscripts = qsubtemp+'/tempscripts'
qsub_exec_path=[]
queue="nmf"
jobs="10"
stru="3200"
jobname="my_job"
offset 
jd2 = True


config_file_path = sys.argv[1]

jobs=int(jobs)
stru=int(stru)
if not os.path.exists(qsubtemp):
	os.mkdir(qsubtemp)
if not os.path.exists(tempscripts):
	os.mkdir(tempscripts)
if not os.path.exists('/common/madsci/Modeling/qsubtemp/tempscripts/'+jobname):
	os.mkdir('/common/madsci/Modeling/qsubtemp/tempscripts/'+jobname)
if not os.path.exists('/common/madsci/Modeling/qsubtemp/'+jobname):
	os.mkdir('/common/madsci/Modeling/qsubtemp/'+jobname)

if len(sys.argv)<2:
	print "No Config file specified. "
	os.abort()

protocol = open(config_file_path, 'r')
config = []
for line in protocol:
	config.append(line)
protocol.close()
outpath = None
for option in config:
	if re.search('-out:path', option):
		ind = config.index(option)
		outpath = config.split()[0]
		config.pop(ind)
rootjob = jobname
for i in range(1, jobs+1):
	script = '/common/madsci/Modeling/qsubtemp/tempscripts/'+rootjob+'/'+jobname+'_'+repr(i)+'.in'
	JOBSCRIPT = open(script, 'w')
	jran = str(1000000+offset)
	suff = str(i)
	config.append('-out:nstruct '+repr(stru)+"\n")
	config.append('-constant_seed -jran '+jran+"\n")
	config.append('-out:prefix '+jobname+"\n")
	config.append('-out:suffix '+suff+"\n")
	JOBSCRIPT.write('#PBS -l nodes=1:ppn=1\n')
	JOBSCRIPT.write('mkdir /scratch/'+jobname+'_'+repr(i)+'\n')
	JOBSCRIPT.write('cd /scratch/'+jobname+'_'+repr(i)+'\n')
	JOBSCRIPT.write(prot+'\n')
	JOBSCRIPT.write('mkdir '+outpath+'\n')
	JOBSCRIPT.write('mkdir '+outpath+'/raw'+'\n')
	JOBSCRIPT.write('cp -r /scratch/'+jobname+'_'+repr(i)+' '+outpath+'/raw'+'\n')
	JOBSCRIPT.write('rm -r /scratch/'+jobname+'_'+repr(i)+'\n')
	JOBSCRIPT.close()
	print 'Kicking Job number '+repr(i)
	offset = offset+10
	run = '/usr/local/bin/qsub -d /common/madsci/Modeling/qsubtemp/'+rootjob+' -N '+jobname+' -V -q '+queue+' '+script
	os.system(run)
