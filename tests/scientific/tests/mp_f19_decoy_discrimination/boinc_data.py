from collections import namedtuple
import os

def read_boinc_dat( ):
	BoincJob = namedtuple('BoincJob','job_id, job_number') 
	boinc_dat = { } # keys = type of run, values are named tuples of job_id and job_number
	boinc_dat_file = open('boinc.dat')
	for line in boinc_dat_file:
		data = line.strip().split()
		if data:
			boinc_dat[ data[2] ] = BoincJob( data[0], data[1] )
	return boinc_dat
		
def write_boinc_dat( arg_boinc ):
	if not os.path.exists('boinc.dat'):
		out = open( 'boinc.dat', 'w' )
		out.write("%s\t%s\t%s" % args.boinc )
	else:
		boinc_data = read_boinc_dat()
		if not arg_boinc[2] in boinc_data:
			out = open( 'boinc.dat', 'a' )
			out.write("\n%s\t%s\t%s" % args.boinc )
