Procedure
1. Renumber the pdb from 1 by executing the script $ROSETTA_SCRIPTS/pdb_renumber.py -n 1 --keep-table --norestart input.pdb output.pdb
	a) The option --preserve can be specified to keep insertion codes and heteroflags, as these options are ignored by rosetta.

Timing

step 1 takes approximately 0.5 seconds per pdb file. 

Troubleshooting

If pdb_renumber.py terminates with an error that either the input file could not be opened, ensure that you have properly typed the file path, and that you have permission to read the file.  If an error that the output file could not be opened occurs, ensure that you have permission to write to the desired directory. 

Anticipated Results

The script produces a PDB file starting from 1 with only the ATOM and TER lines retained.