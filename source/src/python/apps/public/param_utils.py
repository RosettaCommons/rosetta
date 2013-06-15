#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
'''
Various functions for manipulating params files.  

Author: Sam DeLuca
'''
import os
def get_name_from_params(path,database):
    '''Given the path to a params file, return the IO_STRING value.  if it doesn't exist, return None'''
    path=path.strip()
    #if path is commented, ignore it
    if len(path) == 0:
        return 0
    if path[0] == "#":
        return 0
    fullpath = database+path
    #print fullpath
    path = path.split("/")
    #is path a params file?
    type = path[-1].split(".").pop()
    if type == "params":
        if os.path.exists(fullpath):
            paramsfile = open(fullpath,'r')
            for line in paramsfile:
                line = line.split()
                if len(line) >0:
                    if line[0] == "IO_STRING":
                        return line[1]
        else:
            return None
    else:
        return None
        
def get_disallowed_ligands(database):
    '''Return a set of 3 letter names which are already assigned to ligands in the database'''
    residue_type_set_path = database+"chemical/residue_type_sets/"
    residue_set_list = os.listdir(residue_type_set_path)
    disallowed_ligands = set()
    for residue_set in residue_set_list:
        if residue_set[0] != ".":
            current_residue_set_path = residue_type_set_path+residue_set+"/"
            residue_types = open(current_residue_set_path+"residue_types.txt",'r')
            for line in residue_types:
                name = get_name_from_params(line, current_residue_set_path)
                #print name
                if name != None:
                    disallowed_ligands.add(name)
    return disallowed_ligands
    
def rename_param_file(param_path,new_name,new_conformer_path):
    '''Rename a param file residue and update the conformer file path'''
    param_file = open(param_path,'r')
    param_lines = [x.rstrip() for x in param_file]
    param_file.close()
    
    for index,line in enumerate(param_lines):
        fields = line.split()
        if fields[0] == "NAME":
            param_lines[index] = "NAME "+new_name
            continue
        if fields[0] == "IO_STRING":
            param_lines[index] = "IO_STRING "+new_name+" "+fields[2]
            continue
        if fields[0] == "PDB_ROTAMERS":
            param_lines[index] = "PDB_ROTAMERS "+new_conformer_path
    
    param_file = open(param_path,'w')
    for line in param_lines:
        param_file.write(line+"\n")
        
        
def rename_pdb_file(pdb_path,new_name):
    '''Renames all the HETATM resnames in the specified pdb to new_name'''
    pdb_file = open(pdb_path,'r')
    pdb_lines = [x.rstrip() for x in pdb_file]
    pdb_file.close()
    assert(len(new_name) == 3)
    for index,line in enumerate(pdb_lines):
        if len(line) < 6:
            continue
        if line[0:6] == "HETATM":
            pdb_lines[index] = line[:17]+new_name+line[20:]
    
    pdb_file = open(pdb_path,'w')
    for line in pdb_lines:
        pdb_file.write(line+"\n")
