import argparse
import os
import glob

def main( args ):
	if not os.path.exists( args.relax_dir ):
		print "ERROR: %s does not exist!" %(args.relax_dir)
		return
	scores_and_dirs = []

	# find the scores for each relaxed structure
	dirs_to_check = os.listdir( args.relax_dir )
	for rdir in dirs_to_check:
		full_d = args.relax_dir + '/' + rdir
		if not os.path.isdir( full_d ): continue
		if rdir == 'input_lists': continue
		if len( glob.glob( '%s/min_again*pdb' %(full_d) ) ) < 1: continue
		pdb = glob.glob( '%s/min_again*pdb' %(full_d) )[0]
		if not os.path.exists( '%s/relax_2.log' %(full_d)): continue
		score = 0
		for line in open( '%s/relax_2.log' %(full_d)):
			if "Score of the complex" in line: 
				score = float( line.split()[-1].replace('\n',''))
				break
		scores_and_dirs.append( [ score, pdb ] )

	nmodels_to_print = min( len( scores_and_dirs ), args.nmodels )
	lowest_scores_and_dirs = sorted( scores_and_dirs, key=lambda(x): x[0] )[0:nmodels_to_print]
	print "The lowest scoring %d models:" %(nmodels_to_print)
	for i,m in enumerate( lowest_scores_and_dirs ):
		print "Rank %d, Score %0.3f: %s" %(i, m[0], m[1] )

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Get the ddG scores")
	parser.add_argument('-r', '--relax_dir', type=str, default="", help='Directory that contains the results of the relax run')
	parser.add_argument('-n', '--nmodels', type=int, default=20, help='Number of low scoring models to print')
	args = parser.parse_args()
	main( args )
