from util import read_pdb
import os
import argparse

def combine_pdbs(full_pdb, subset_pdbs, output_pdb):

	fout = open( output_pdb, 'w')
	(coords, pdb_lines, sequence, chains, residues) = read_pdb( full_pdb )
	# concatenate all the subset PDBs together
	string_of_subset_pdbs = ''
	for s in subset_pdbs:
		string_of_subset_pdbs+=s
		string_of_subset_pdbs+=' '
	os.system('cat %s > tmp_subsets.pdb' %(string_of_subset_pdbs))
	
	(s_coords, s_pdb_lines, s_sequence, s_chains, s_residues) = read_pdb( 'tmp_subsets.pdb' )

	#loop through the full pdb
	full_chains = pdb_lines.keys()
	full_chains.sort()
	for c in full_chains:
		resis = sorted(pdb_lines[c].keys())
		for r in resis:
			write_subset_resi = False
			atoms_to_compare = []
			if c in s_pdb_lines.keys():
			 if r in s_pdb_lines[c].keys():
			  if sequence[c][r]!=s_sequence[c][r]: write_subset_resi = True
			  else:
			   for s_atom in s_pdb_lines[c][r].keys():
				  atoms_to_compare.append(s_atom)
			if write_subset_resi:
				for atom in s_pdb_lines[c][r].keys():
					fout.write( s_pdb_lines[c][r][atom]+'\n' )
			else:
				for atom in pdb_lines[c][r].keys():
					# check if this exists in the subset pdb
					if atom in atoms_to_compare:
						fout.write( s_pdb_lines[c][r][atom]+'\n' )
					else:
						fout.write( pdb_lines[c][r][atom]+'\n' )

	os.remove('tmp_subsets.pdb')
	fout.close()
	
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Combine a full pdb with rebuilt subsets")
	parser.add_argument('-f', '--full_pdb', type=str, help='The full/original pdb file')
	parser.add_argument('-s', '--subset_pdbs', type=str, nargs='+', help='The subset files to overwrite atoms in the full pdb')
	parser.add_argument('-o', '--output_pdb', type=str, help='Name for the output pdb file')
	args = parser.parse_args()
	combine_pdbs(args.full_pdb, args.subset_pdbs, args.output_pdb)
