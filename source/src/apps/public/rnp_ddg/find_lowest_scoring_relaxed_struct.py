import os
import argparse

def main( args ):
	os.system('grep "Score of the complex" %s/*/relax_2.log > %s/relaxed_complex_scores.txt' %(args.relax_dir, args.relax_dir))
	with open('%s/relaxed_complex_scores.txt' %(args.relax_dir), 'r') as fil:
		min_score = 100000.
		min_dir = ''
		for line in fil.readlines():
			score = float(line.split()[5])
			if score < min_score:
				min_score = score
				min_dir = line.split(':')[0].replace('relax_2.log','')
	print "Min struct is in %s" %(min_dir)
	## Copy the min struct
	os.system('cp %s/min_again_*_bound.pdb %s/min_relaxed_%s.pdb' %(min_dir,args.relax_dir,args.start_struct_pref))

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Get the lowest scoring relaxed structure")
	parser.add_argument('-r', '--relax_dir', type=str, default="", help='Relax directory')
	parser.add_argument('-s', '--start_struct_pref', type=str, default="", help='Starting starting structure file prefix')
	args = parser.parse_args()
	main( args )
