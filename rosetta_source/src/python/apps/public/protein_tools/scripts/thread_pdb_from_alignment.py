#!/usr/bin/env python2.5
from Bio import AlignIO
import Bio.PDB
from optparse import OptionParser
import sys
import array
from rosettautil.protein import util
from rosettautil.protein import alignment
from rosettautil.protein import pdbStat
from rosettautil.util import fileutil


def list_to_generator(list):
    for x in list:
        yield x


usage = "%prog [options] alignment_file.aln template.pdb output.pdb"
parser=OptionParser(usage)
parser.add_option("--template",dest="template",help="name of the template sequence", default="template")
parser.add_option("--target", dest ="target",help="name of the target sequence",default="target")
parser.add_option("--chain",dest="chain",help="chain to thread pdb around",default="A")
parser.add_option("--align_format",dest="align_format",help="alignment file format, choose from clustal, emboss, fasta, fasta-m10,ig,nexus,phylip,stockholm.  See http://biopython.org/wiki/AlignIO for details",default="clustal")
(options,args)= parser.parse_args()
if len(args) != 3:
    parser.error("you must specify an alignment file, template pdb, and output pdb")

#read in our input files
alignment_file = fileutil.universal_open(args[0],'rU')
alignment_data = AlignIO.read(alignment_file,options.align_format)
alignment_file.close()
template_struct = util.load_pdb(args[1])

#if len(alignment_data)  != 2:
#    sys.exit("alignment file must have exactly 2 sequences!") 

#find all the gaps, get numeric IDs from the string tags in the alignment file
try:
    template_gaps = alignment.find_gaps(alignment_data,options.template)
except LookupError:
    sys.exit("could not find "+options.template+" in alignment file")
try:    
    target_gaps = alignment.find_gaps(alignment_data,options.target)
except LookupError:
    sys.exit("could not find "+options.target+" in alignment file")
try:
    template_id = alignment.get_id_from_tag(alignment_data,options.template)
except LookupError:
    sys.exit("could not find "+options.template+" in alignment file")
try:
    target_id = alignment.get_id_from_tag(alignment_data,options.target)
except LookupError:
    sys.exit("could not find "+options.target+" in alignment file")

#there might be missing density in our pdb file.  Align the template sequence to the pdb file
#and return a gapped sequence.  We will use this in conjunction with the template and alignment sequence
template_sequence = alignment_data[alignment.get_id_from_tag(alignment_data,options.template)]
gapped_template = pdbStat.find_gaps(template_struct,template_sequence,options.chain)

#if you have an alignment gap thats larger than 1 aa, this script won't work
for gap in target_gaps:
    if abs(gap[1]-gap[0]) > 1:
		print "WARNING: gap of size "+ str(gap[0]-gap[1])+" in target sequence. " 
		print "We cannot completely thread this protein in an automatic way, "
		print "manual inspection and adjustment of loop files will be required."

#we need to make a new structure, then a new model, then a new chain, then we fill the chain with residues, and atoms
output_structure_builder = Bio.PDB.StructureBuilder.StructureBuilder()
output_structure_builder.init_structure(args[2]) 
output_structure_builder.init_model(1) #there is only one model
output_structure_builder.init_chain(options.chain) #there is only one chain, same ID as the template
output_structure_builder.init_seg("")



#thats it for the initialization stuff, now we go through add residues, and atoms.
template_residues = None
chain_found = False
for chain in template_struct.get_chains():
    if chain.get_id() == options.chain:
        chain_found = True
        template_residues= list_to_generator(chain.get_list())
        break

if not chain_found:
	sys.exit("ERROR: You specified chain "+options.chain+" but this chain does not exist in the pdb file "+ args[1])

#template_residues = template_struct.get_chains()
sequence_num = 1 #the pdb sequence number
atom_num = 1 #the atom id
for align_resn, temp_resn,gap_temp_resn in zip(alignment_data[target_id],alignment_data[template_id],gapped_template):
    #print align_resn, temp_resn
    if align_resn == '-' and temp_resn == '-':  #this shouldn't happen, but it is safe to ignore
        continue
    elif align_resn != '-' and (temp_resn == '-' or gap_temp_resn == '-'): #gap in the template (or pdb), not in the alignment, build a loop
        align_name3 = Bio.PDB.Polypeptide.one_to_three(align_resn)
        output_structure_builder.init_residue(align_name3," ",sequence_num," ")
        zero_triplet = array.array('f',[0.0,0.0,0.0])
        output_structure_builder.init_atom("N",zero_triplet,0.0,-1.0," "," N  ", atom_num, "N")
        atom_num += 1
        output_structure_builder.init_atom("CA",zero_triplet,0.0,-1.0," "," CA ",atom_num,"C")
        atom_num += 1
        output_structure_builder.init_atom("C",zero_triplet,0.0,-1.0," "," C  ",atom_num,"C")
        atom_num += 1
        output_structure_builder.init_atom("O", zero_triplet,0.0,-1.0," "," O  ",atom_num,"O")
        atom_num += 1
        sequence_num += 1
    elif align_resn == '-' and temp_resn != '-': #gap in the alignment, not in the template, skip the residue
        template_residues.next() #pull a residue out of the pdb and throw it away
        continue
    elif align_resn != '-' and temp_resn != '-': #we're aligned, copy backbone from old pdb to new, if the sidechain is identical, copy that too
        align_name3 = Bio.PDB.Polypeptide.one_to_three(align_resn)
        temp_name3 = Bio.PDB.Polypeptide.one_to_three(temp_resn)
        current_res = template_residues.next() #pull the next residue out of the pdb file 
        if(current_res.get_resname() != temp_name3):  #if the current residue from the pdb isnt the same type as the current from the template, something's broken
            print current_res.get_resname(),temp_name3
            sys.exit("Residue mismatch between alignment and PDB, check that PDB sequence and alignment sequence are identical")
        output_structure_builder.init_residue(align_name3," ",sequence_num," ")
        if(align_name3 == temp_name3): #if we have an exact alignment, copy all the atoms over, including the sidechain
            for atom in current_res:
                coords = atom.get_coord()
                name = atom.get_name()
                element=""
                if name[0].isdigit(): #if the first char of the atomname is a digit, its an H, otherwise its whatever the first char is
                    element = "H"
                else:
                    element = name[0]
                fullname = atom.get_fullname()
                output_structure_builder.init_atom(name,coords,0.0,1.0," ",fullname,atom_num,element)
                atom_num += 1
        else: #we don't have an exact alignment, just copy the backbone coordinates
            #by definition, the first 4 atoms that come out of a residue are N, CA, C, O.  How convenient...
            for atom_name in ["N","CA","C","O"]:
                atom = current_res[atom_name]
                coords = atom.get_coord()
                name = atom.get_name()
                element=""
                if name[0].isdigit(): #if the first char of the atomname is a digit, its an H, otherwise its whatever the first char is
                    element = "H"
                else:
                    element = name[0]
                fullname = atom.get_fullname()
                output_structure_builder.init_atom(name,coords,0.0,1.0," ",fullname,atom_num,element)
                atom_num += 1
        sequence_num += 1

#now, output the structure
output_struct = output_structure_builder.get_structure()
pdb_io = Bio.PDB.PDBIO()
pdb_io.set_structure(output_struct)
outfile = fileutil.universal_open(args[2],'w')
pdb_io.save(outfile)
outfile.close()
