import os
import glob
import argparse

ROSETTA_SCAL_FACT = 0.38
NN_SCAL_FACT = 0.28

def main( args ):

	# First find all the wt scores
	# Find score in 0 dir (wt)
	start_dir = os.getcwd()
	run_dir = os.path.abspath( args.run_dir )
	os.chdir( run_dir + '/0' )
	complex_score, wt_protein_score, rna_score = get_scores( 0.0 )
	dG_wt_0 = ROSETTA_SCAL_FACT*complex_score - ROSETTA_SCAL_FACT*wt_protein_score - NN_SCAL_FACT*rna_score
	ddG = 0.0
	with open('ddG_score.txt', 'w') as df:
	# write: ddG dG complex_score protein_score rna_score
		df.write('%f %f %f %f %f\n' %(ddG, dG_wt_0, complex_score, wt_protein_score, rna_score))

	print "####################################"
	print "RESULTS:\n"
	print "WT ddG: %0.2f kcal/mol" %(ddG)

	# Loop through all wt_* dirs
	for d in glob.glob( run_dir + '/wt_*' ):
		#if 'wt_23' in d or 'wt_24' in d or 'wt_25' in d or 'wt_26' in d: continue
		if not os.path.isdir( d ):
			continue
		os.chdir( d )
		complex_score, protein_score, rna_score = get_scores( wt_protein_score )
		dG = ROSETTA_SCAL_FACT*complex_score - ROSETTA_SCAL_FACT*protein_score - NN_SCAL_FACT*rna_score
		with open('dG_score.txt', 'w') as df:
		# write dG complex_score protein_score rna_score
			df.write( '%f %f %f %f\n' %(dG, complex_score, protein_score, rna_score))
		print "%s ddG: %0.2f kcal/mol" %(d, ddG)

	# parse the sequences
	seqs = {}
	index = 0
	for line in open(start_dir + '/' + args.seq_file):
		index += 1
		seqs[ str(index) ] = line.replace('\n','')

	# Next find all the other mutant scores
	# Loop through remaining directories
	for d in os.listdir( run_dir ):
		#if d == '36': continue
		#if d == 'test_protein_repacking': continue
		if "wt_" in d: continue
		if d == '0': continue
		if not os.path.isdir( run_dir + '/' + d ): #only look in directories
			continue
		#go into the directory
		os.chdir( run_dir + '/' + d )
		complex_score, protein_score, rna_score = get_scores( wt_protein_score )

		# Calculate the dG value
		dG = ROSETTA_SCAL_FACT*complex_score - ROSETTA_SCAL_FACT*protein_score - NN_SCAL_FACT*rna_score

		# Then check if there is a wt_comparison file (should exist for high res)
		if os.path.exists('wt_comparison.txt'):
			with open('wt_comparison.txt','r') as wfil:
				wt_dir = wfil.readline().replace('\n','')
			# Go to the wt_dir and read the dG
			with open('../%s/dG_score.txt' %(wt_dir), 'r') as gwf:
				dG_wt = float( gwf.readline().split()[0].replace('\n',''))
				
		else:
			dG_wt = dG_wt_0

		ddG = dG - dG_wt
		# Write file
		with open('ddG_score.txt', 'w') as gf:
		# write: ddG dG complex_score protein_score rna_score
			gf.write( '%.5f %.5f %.5f %.5f %.5f\n' %(ddG, dG, complex_score, protein_score, rna_score))

		if d not in seqs.keys():
			print "Warning: problem with your sequence file..."
			print "continuing."
			continue
		print "%s ddG: %0.2f kcal/mol" %(seqs[d], ddG)

	print "\n####################################"
	print "Results are also written to the ddG_score.txt files in "
	print "each mutant directory in %s" %(args.run_dir)
	print "The format is:"
	print "ddG dG complex_score protein_score rna_score"
	print "####################################"


def get_scores( wt_protein_score ): # works in the current directory
	# Get the complex score
	os.system('grep "Average of top" get_complex_scores.log > complex_score.txt')
	with open('complex_score.txt', 'r') as cfil:
		line = cfil.readline()
		if len(line) < 2: 
			print "NO COMPLEX score in %s" %(os.getcwd())
			return 0.0, 0.0, 0.0
		complex_score = float(line.split()[-1].replace('\n',''))
	# Get the protein score
	os.system('grep "Average unbound protein score" get_complex_scores.log > protein_scores.txt')
	protein_scores = []
	with open('protein_scores.txt','r') as pfil:
		for line in pfil:
			pscore = line.split()[-1].replace('\n','')
			protein_scores.append( float(pscore) )
	if len( protein_scores ) == 0: protein_scores.append( wt_protein_score )
	# get the average of average protein scores
	total_protein_score = 0.0
	n_protein_scores = len( protein_scores )
	for pscore in protein_scores:
		total_protein_score += pscore
	avg_protein_score = total_protein_score / n_protein_scores
	#avg_protein_score = np.mean( np.array(protein_scores))

	# get the RNA score (this is the free energy)
	#rna_score = np.genfromtxt('rna_nn_E.txt', skip_header=1, usecols=0)
	for line in open('rna_nn_E.txt'):
		if "FE" in line: continue # skip the header
		rna_score = float( line.split()[0] )
		break

	return complex_score, avg_protein_score, rna_score



if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Get the ddG scores")
	parser.add_argument('-r', '--run_dir', type=str, default="", help='Directory that contains all mutant directories')
	parser.add_argument('--seq_file', type=str, default="", help='Text file containing sequences to be modeled (the file provided to general_RNP_setup_script)')
	args = parser.parse_args()
	main( args )
