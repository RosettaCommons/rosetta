import os
import glob
import argparse
import numpy as np

def main( args ):

	# First find all the wt scores
	# Find score in 0 dir (wt)
	run_dir = os.path.abspath( args.run_dir )
	os.chdir( run_dir + '/0' )
	complex_score, wt_protein_score, rna_score = get_scores( 0.0 )
	dG_wt_0 = 0.5*complex_score - 0.5*wt_protein_score - rna_score
	ddG = 0.0
	with open('ddG_score.txt', 'w') as df:
	# write: ddG dG complex_score protein_score rna_score
		df.write('%f %f %f %f %f\n' %(ddG, dG_wt_0, complex_score, wt_protein_score, rna_score))

	# Loop through all wt_* dirs
	for d in glob.glob( run_dir + '/wt_*' ):
		if not os.path.isdir( d ):
			continue
		os.chdir( d )
		complex_score, protein_score, rna_score = get_scores( wt_protein_score )
		dG = 0.5*complex_score - 0.5*protein_score - rna_score
		with open('dG_score.txt', 'w') as df:
		# write dG complex_score protein_score rna_score
			df.write( '%f %f %f %f\n' %(dG, complex_score, protein_score, rna_score))

	# Next find all the other mutant scores
	# Loop through remaining directories
	for d in os.listdir( run_dir ):
		if "wt_" in d: continue
		if d == '0': continue
		if not os.path.isdir( run_dir + '/' + d ): #only look in directories
			continue
		#go into the directory
		os.chdir( run_dir + '/' + d )
		complex_score, protein_score, rna_score = get_scores( wt_protein_score )

		# Calculate the dG value
		dG = 0.5*complex_score - 0.5*protein_score - rna_score

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



def get_scores( wt_protein_score ): # works in the current directory
	# Get the complex score
	os.system('grep "Average of top" get_complex_scores.log > complex_score.txt')
	with open('complex_score.txt', 'r') as cfil:
		complex_score = float(cfil.readline().split()[-1].replace('\n',''))
	# Get the protein score
	os.system('grep "Average unbound protein score" get_complex_scores.log > protein_scores.txt')
	protein_scores = []
	with open('protein_scores.txt','r') as pfil:
		for line in pfil:
			pscore = line.split()[-1].replace('\n','')
			protein_scores.append( float(pscore) )
	if len( protein_scores ) == 0: protein_scores.append( wt_protein_score )
	# get the average of average protein scores
	avg_protein_score = np.mean( np.array(protein_scores))

	# get the RNA score (this is the free energy)
	rna_score = np.genfromtxt('rna_nn_E.txt', skip_header=1, usecols=0)

	return complex_score, avg_protein_score, rna_score



if __name__ == '__main__':
	print "\n#################################################"
	print "WARNING: Please use updated versions of these scripts in apps/public/rnp_ddg/"
	print "#################################################\n"
	parser=argparse.ArgumentParser(description="Get the ddG scores")
	parser.add_argument('-r', '--run_dir', type=str, default="", help='Directory that contains all mutant directories')
	args = parser.parse_args()
	main( args )
