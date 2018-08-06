curl -o ./components.cif ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
python split_components.py
rm components.cif
