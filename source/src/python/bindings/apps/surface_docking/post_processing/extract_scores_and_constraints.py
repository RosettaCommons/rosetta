#!/usr/bin/python2.6
# Script to extract scores from score table and constraint scores appended to end of PDB
# Emily Koo

import os
from sys import argv

def parse(pdbs):
    savefile = pdbs.name + "_constraints"
    save = open(savefile, 'a')
    for pdb in pdbs:
        decoy = pdb.strip()
        file = open(decoy, 'r')
        print decoy
        for line in file:
            if line.startswith('Constraint Energies'):
                newline = file.next()
                # Atom pair 4
                energy = float(newline.split()[3])
                # Dihedral 3
                newline = file.next()
                energy += float(newline.split()[2])
                save.write(decoy + "\t" + str(energy) + "\n")
                break
            
    save.close()

state = argv[1]

ads = 'ls AdsState*.pdb > ads_pdb'
sol = 'ls SolState*.pdb > sol_pdb'
get = "grep 'Total weighted' SolState*.pdb | awk '{print $1, $5}' | sed s/://g | sort > SolState.sorted"
get2 = "grep 'Total weighted' AdsState*.pdb | awk '{print $1, $5}' | sed s/://g | sort > AdsState.sorted"
join = 'join sol_pdb_constraints SolState.sorted > sol_energy_cmp'
join2 = 'join ads_pdb_constraints AdsState.sorted > ads_energy_cmp'

if state == "ads":
    os.system(ads)
    ads_pdb = open('ads_pdb', 'r')
    parse(ads_pdb)
    ads_pdb.close()
    os.system(get2)
    os.system(join2)
elif state == "sol":
    os.system(sol)
    sol_pdb = open('sol_pdb', 'r')
    parse(sol_pdb)
    sol_pdb.close()
    os.system(get)
    os.system(join)
else:
    os.system(sol)
    os.system(ads)
    sol_pdb = open('sol_pdb', 'r')
    ads_pdb = open('ads_pdb', 'r')
    parse(sol_pdb)
    parse(ads_pdb)
    sol_pdb.close()
    ads_pdb.close()

    os.system(get)
    os.system(get2)
    os.system(join)
    os.system(join2)
