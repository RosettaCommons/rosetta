#!/usr/bin/env python 

#author: Jordan Willis
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
import sys
import glob

def usage():
    print "Pulls fasta file of a specific chain up to a specified output file"
    print "python get_fasta_from_pdb_by_chain.py <pdb> <chain> <output>"
    
    
#import warnings

#def fxn():
#    warnings.warn("deprecated", PDBConstructionWarning)

#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#    fxn()



if (len(sys.argv) < 3):
   usage()
   exit()

parser = PDBParser(PERMISSIVE=1)

#if sys.argv[3]:
#    stop = int(sys.argv[3])
#else:
stop = ""

open(sys.argv[3], 'w').write("")

def get_three_letter(pdb, chain_id):
    fasta_three = []
    for model in pdb.get_list():
       for chain in model.get_list():
            if chain.get_id() == chain_id:
               for resi in chain.get_list():
                   fasta_three.append(resi.get_resname())
    return fasta_three
def get_one_letter(list_of_three):
    fasta_one=[]
    for x in list_of_three:
       x=three_to_one(x)
       fasta_one.append(x)
    return fasta_one


def output_to_file(one_letter, file):
    file_handle = open(sys.argv[3], 'a')
    file_handle.write(">" + file.split('.')[0] + "\n")
    file_handle.write("".join(one_letter) + "\n")
    file_handle.close()


#for file in open(sys.argv[1]).readlines():
    
handle = open(sys.argv[1])
struct = parser.get_structure("random",handle)
three_letter = get_three_letter(struct, sys.argv[2])
one_letter = get_one_letter(three_letter)
output_to_file(one_letter,sys.argv[1])

