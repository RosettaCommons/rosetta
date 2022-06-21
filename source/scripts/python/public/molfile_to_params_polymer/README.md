##################
#
# Refactored scripts related to creating params files for polymer
# residues (e.g., non-canonical amino acids)
#
# Original materials written by Ian Davis
#
# Modified by:
# Oanh Vu
# Benjamin P. Brown
# 
#################

NOTES
- The original "molfile_to_params_polymer.py" script was divided into smaller, more compartmentalized scripts:
atom_functions.py
atom_functions.pyc
bond_functions.py
bond_functions.pyc
fragment_functions.py
fragment_functions.pyc
IO_functions.py
IO_functions.pyc
polymer_functions.py
polymer_functions.pyc

- The run script was renamed and modified to use the above scripts:
Rosetta/main/source/scripts/python/public/molfile_to_params_polymer.py

- New functionality was added to:
	- Take SDF conformational ensemble as input and generate PDB rotamers
	- Check that partial charge sums are equivalent to formal charges after charge correction.
	- Added option to use pre-calculated partial charges
	- Minor bug changes throughout
	- Write explicit aromatic bond types or single/double bond orders
	- Assigned mm atom types to halogens based on fa_standard CHARMM params
