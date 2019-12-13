import argparse
import DRRAFTER_util

# fix the numbering of the final models 
def finalize( args ):

	resnum_tag = ''
	for line in open( args.fasta ):
		if line.startswith( '>' ):
			resnum_tag_list = line.replace('> ', '' ).replace('>', '' ).split(' ')[1:]
			for x in resnum_tag_list:
				resnum_tag += x + ' '
	
	resnums = []
	chains = []
	segids = []
	for x in resnum_tag.split(' '):
		DRRAFTER_util.get_resnum_chain( x, resnums, chains, segids )
		
	for i in range( 1, 11 ):
		pdb_name = args.out_pref + '_all_models_all_fits_' + args.final_round + '.out.' + str( i ) + '.pdb'
		DRRAFTER_util.renumber_pdb( [pdb_name], resnums, chains, segids )
	print( "Done finalizing models" )
		
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="")
	parser.add_argument( '-fasta', type=str, required=True )
	parser.add_argument( '-out_pref', type=str, required=True )
	parser.add_argument( '-final_round', type=str, required=True )
	args = parser.parse_args()
	finalize( args )
