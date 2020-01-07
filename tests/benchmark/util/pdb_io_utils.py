

import numpy as np

THREE2ONE = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



def ca_pdb_reader(pdb_lines):
    def read_line(line):
        atom_type = line[11:16]
        OLC = THREE2ONE[line[17:20]]
        xyz = np.array(
            (
                float(line[26:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        )
        residue_number = int(line[22:26])
        return residue_number, atom_type, xyz, OLC
    CAs = {}
    CA_olc = {}
    for line in pdb_lines:
        if not line or not line.startswith("ATOM "):
            continue
        residue_number, atom_type, xyz, OLC = read_line(line)
        if ' CA ' not in atom_type:
            continue
        CAs[residue_number] = xyz
        CA_olc[residue_number] = OLC
    return CAs, CA_olc
