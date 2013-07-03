#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/pdbs.py
## @brief  general PDB functions for the tookit.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports

#Python Imports
import os
import re
import copy

#Tkinter Imports
import tkFileDialog

#Toolkit Imports
from window_main import global_variables


class PythonPDB:
    def __init__(self, pdb_file_path):
        """
        Lightweight PDB class specifically for use in CleanPDB window, but functions in quickly manipulating a PDB.
        For more PDB tools, see PDB module.
        """

        self.elements = ("id", "atom_number", "atom_name", "alternate_location", \
                         "three_letter_code", "chain", "residue_number", "i_code", "x", "y", "z", \
                         "occupancy", "b_factor", "element", "charge")
        self.pdb_file_path = pdb_file_path
        self.pdb_map = dict(); #[int line]:[string element]:[string value]
        self.read_pdb_into_map()
        
    def set_pdb_map(self, pdb_map):
        self.pdb_map = pdb_map
        
    def get_pdb_map(self):
        return self.pdb_map

    def read_pdb_into_map(self):
        """
        Reads PDB file path into a basic PDB map.  All data is held as strings.
        """
        
        FILE = open(self.pdb_file_path, 'r')
        i = 1
        for line in FILE:
            line = line.strip('\n')
            if not (re.search("ATOM", line[0:6]) or re.search("HETATM", line[0:6])):
                continue
            
            if not self.pdb_map.has_key(i):
                self.pdb_map[i]=dict()
                
            self.pdb_map[i]["id"]=line[0:6].strip()
            self.pdb_map[i]["atom_number"]=line[6:11].strip();     self.pdb_map[i]["atom_name"] = line[12:16]
            self.pdb_map[i]["alternate_location"]=line[16];        self.pdb_map[i]["three_letter_code"] = line[17:21].strip()
            self.pdb_map[i]["chain"] = line[21].strip();           self.pdb_map[i]["residue_number"]= line[22:26].strip()
            self.pdb_map[i]["i_code"] = line[26];                  self.pdb_map[i]["x"] = line[27:38].strip()
            self.pdb_map[i]["y"]= line[38:46].strip();             self.pdb_map[i]["z"]= line[46:54].strip()
            self.pdb_map[i]["occupancy"] = line[54:60].strip();    self.pdb_map[i]["b_factor"]=line[60:66].strip()
            self.pdb_map[i]["element"]=line[66:78].strip();        self.pdb_map[i]["charge"]=line[78:79].strip()

            i +=1
        FILE.close()
    
    def remove_hetatm_atoms(self):
        temp_pdb_map = self.pdb_map.copy()
        for num in temp_pdb_map:
            if self.pdb_map[num]["id"]=="HETATM":
                del self.pdb_map[num]
                
                
    def remove_element_column(self):
        """
        Removes the extra stuff in the element column, but not the element itself.
        """
        for i in range(1, len(self.pdb_map)+1):
            ele = self.pdb_map[i]["element"]
            e = ele[11]
            self.pdb_map[i]["element"]="           "+e
        print "Extra stuff in Element Columns Removed"
        return self.pdb_map
    
    def remove_waters(self):
        """
        Removes waters from pdb_map
        """
        #codes = ["HOH","TP3","TP5","TIP3","TIP5"]
        temp_pdb_map = self.pdb_map.copy(); #This is to pop elements
        for num in temp_pdb_map:
            if self.pdb_map[num]["three_letter_code"] in ["HOH","TP3","TP5","TIP3","TIP5"]:
                #self.pdb_map.pop(num)
                del self.pdb_map[num]
                
    def remove_alternate_residues(self):
        """
        Removes any alternate residue codes and renumbers by renumbering from 1 and integrating any inserts. 
        """
        
        def get_residue_num(num): return int(self.pdb_map_copy[num]["residue_number"])
        def set_residue_num(num, resnum): self.pdb_map[num]["residue_number"]=str(resnum)
        def get_chain(num):return self.pdb_map_copy[num]["chain"]
        def get_i_code(num):return self.pdb_map_copy[num]["i_code"]
        
        def check_id(num):
            if self.pdb_map_copy[num]['id']=="ATOM":
                return True
            else:
                return False
            
        def check_new_residue(old_num, num, insert_residue=False, pdb_map = False):
            if insert_residue:
                if get_i_code(old_num)==get_i_code(num):
                    return False
                else:
                    return True
            else:
                if get_residue_num(old_num)==get_residue_num(num):
                    return False
                else:
                    return True
        
        def check_new_chain(old_num, num):
            if get_chain(old_num)==get_chain(num):
                return False
            else:
                return True
        
        def check_insertion(num):
            if not get_i_code(num)==" ":
                return True
            else:
                return False
        
        def renumber_from_one(chain_only, start_num):
            resnum = 1
            for num in sorted(chain_only):
                
                insert = check_insertion(num)
                
                #print repr(get_residue_num(num))+":"+repr(insert)
                
                #This is so we don't check if it's a new residue with num-1 - Which won't actually be part of the chain!
                if num==start_num:
                    set_residue_num(num, resnum)

                    
                
                #Iterate resnum if new residue
                elif check_new_residue(num-1, num, insert):
                    resnum+=1
                    set_residue_num(num, resnum)

                else:
                    set_residue_num(num, resnum)
            
            #Set i code at the end, so we can tell if we have new residues or not.
            for num in sorted(chain_only):
                self.pdb_map[num]["i_code"]=" "
                
        def renumber_from_insert(chain_only, start_num):
            pass
        
        self.pdb_map_copy = copy.deepcopy(self.pdb_map)
        
        #Get chains with insertion codes - Now renumbers all chains. Will be an option later.
        chains_with_inserts = dict(); 
        for num in sorted(self.pdb_map):
            #if get_i_code(num)==" ":
            chains_with_inserts[get_chain(num)]=True

        
        #Iterate through all lines/atoms
        #Initialize for scope
        start_residue=0;
        new_start=False
        for chain in chains_with_inserts:
            print "Renumbering chain "+chain
            chain_only=dict()
            for num in sorted(self.pdb_map):
                if chain == get_chain(num) and check_id(num):
                    chain_only[num]=self.pdb_map[num]
            lines = sorted(chain_only)
            res_start = get_residue_num(lines[0])
            
            renumber_from_one(chain_only, lines[0])
                
            #For now, we only renumber from one.
            #else:
                #chain_only = renumber_from_insert(chain_only, lines[0])      
                    
                
            
            
        
    def change_occupancy(self):
        """
        Changes ALL occupancies in a PDB dictionary to 1.00
        Returns PDB Dictionary.
        """
        
        check = 0
        for key in self.pdb_map:
            if self.pdb_map[key]["occupancy"].rfind("0.00")!=-1:
                print "Changing occupancy of residue " + self.pdb_map[key]["residue_number"] + "To 1.00"
                check =1
            self.pdb_map[key]["occupancy"] = "  1.00"
        if check ==1:
            print "Occupancy Column OK for PyRosetta..."
    
    def save_PDB(self, filename=False):
        """
        Uses a the pdb_map to save the data as a PDB file.
        """

        if filename==False:
            filename = tkFileDialog.asksaveasfilename(mode = 'w', initialdir=global_variables.current_directory,title='Save As...')
            if not filename:return
            global_variables.current_directory = os.path.dirname(filename)

        FILE = open(filename, 'w')
        for line_number in self.pdb_map:
            line = self.morph_line_in_pdb_map_to_pdb_line(line_number)
            FILE.write(line+"\n")
        FILE.close()
        print "PDB File Written..."
        return filename
    
    def morph_line_in_pdb_map_to_pdb_line(self, line_num):
        """
        Oh What fun. ;)
        Magic Numbers?: (6,5,4,3,1,4,8,8,8,4,5);
        """
        
        
        #Here we fix the formating of atom name. If we stripped the atom name.
        """
        atom_name = self.pdb_map[line_num]['atom_name']
        if len(atom_name)==1:
            atom_name=' '+atom_name+'  '
        elif len(atom_name)==2:
            #Note that 2 letter elements like CA (calcium) differ from CA (C-Alpha)
            #If calcium, would go @column 13.  if C-Alpha, column 14.
            atom_name=' '+atom_name+' '
        elif len(atom_name)==3:
            atom_name=' '+atom_name
        elif len(atom_name)==4:
            atom_name=atom_name
        else:
            print "Atom Name missing.  Inserting spaces."
            atom_name = '    '
        """
        
        #Create the PDB line.
        line = (self.pdb_map[line_num]['id']).ljust(6)+             (self.pdb_map[line_num]['atom_number']).rjust(5)+" "+ self.pdb_map[line_num]['atom_name']+ \
               (self.pdb_map[line_num]['alternate_location'])+      ((self.pdb_map[line_num]['three_letter_code']).rjust(3)).ljust(4)+  (self.pdb_map[line_num]['chain'])+             \
               (self.pdb_map[line_num]['residue_number']).rjust(4)+ (self.pdb_map[line_num]['i_code']) +                                              \
               (self.pdb_map[line_num]['x']).rjust(11)+             (self.pdb_map[line_num]['y']).rjust(8)+                  (self.pdb_map[line_num]['z']).rjust(8) +   \
               (self.pdb_map[line_num]['occupancy']).rjust(6)+      (self.pdb_map[line_num]['b_factor']).rjust(6)

        #Note three letter code is wonky due to DA residues.  ljust(4) was not working.
        return line
    
    def pdb_alias(self, pairs, element):
        """
        Replaces ALL occurances of old element with new from pair.
        pair is a dictionary. In C++ it would be an array of pairs.  [string old]:[string new]
        For Specific functions, please see below.
        """
        for num in self.pdb_map:
            for old in pairs:
                if self.pdb_map[num][element]==old:
                    self.pdb_map[num][element]=pairs[old]
                    
    def pdb_atom_alias(self, line_num, pair):
        """
        Replaces atom_names with ones Rosetta is happy with.
        pair is a dictionary. In C++ it would be an array of pairs.  [string MD atom_name]:[string rosetta atom_name]
        """
        for start in pair:
            if self.pdb_map[line_num]["atom_name"]==start:
                print self.pdb_map[line_num]["three_letter_code"]+":"+self.pdb_map[line_num]["atom_name"]+":"+ pair[start]
                self.pdb_map[line_num]["atom_name"]=pair[start]
    
    def pdb_residue_alias(self, pairs):
        """
        Replaces ALL occurances of old residue with new residue.
        pair is a dictionary. In C++ it would be an array of pairs.  [string old residue_name]:[string new residue_name]
        """
        for num in self.pdb_map:
            for old in pairs:
                if self.pdb_map[num]["residue_name"]==old:
                    self.pdb_map[num]["residue_name"]=pairs[old]
    
    def pdb_chain_alias(self, pairs):
        """
        Replaces ALL occurances of old chain with new chain.
        pair is a dictionary. In C++ it would be an array of pairs.  [string old chain]:[string new chain]
        """
        for num in self.pdb_map:
            for old in pairs:
                if self.pdb_map[num]["chain"]==old:
                    self.pdb_map[num]["chain"]=pairs[old]
            
    
    def remove_chain(self, chain):
        """
        Removes chain from pdb_map
        """
        temp_pdb_map = self.pdb_map.copy()
        for num in temp_pdb_map:
            if self.pdb_map[num]["chain"]==chain:
                del self.pdb_map[num]
    
    def read_file_and_replace_b_factors(self, deliminator, filename="", resnum_column=1, chain_column=2, data_column=3, atomname_column=False):
        """
        This function reads a deliminated file with data and inserts the data into the BFactor column.  Used to visualize arbitrary data.
        Use function options to control which column the data is in as well as where your resnums and chains are located.
        If atomname column is given, will insert by atom instead of by residue
        """
        
        if not filename:
            filename = tkFileDialog.askopenfilename(title="Data file", initialdir=global_variables.current_directory)
            if not filename:return
            global_variables.current_directory = os.path.dirname(filename)
        
        INFILE = open(filename, 'r')
        for line in INFILE:
            if line[0] == "#":continue
            line = line.strip()
            lineSP = line.split(deliminator)
            if len(lineSP)<3:
                print "Could not read line.  Must have resnum, chain, and data columns"
                continue
            if not atomname_column:
                self.replace_residue_b_factor(lineSP[resnum_column-1], lineSP[chain_column-1], lineSP[data_column-1])
            else:
                if len(lineSP)<4:
                    print "Could not read line.  Must have resnum, chain, atomname, and data columns"
                    continue
                self.replace_atom_b_factor(lineSP[resnum_column-1], lineSP[chain_column-1], lineSP[atomname_column-1], lineSP[data_column-1])
        INFILE.close()
        
    def replace_residue_b_factor(self, resnum, chain, data):
        """
        Replaces the b factor of each atom in the residue with data.
        Can be all string representations or not.
        """
        
        if type(resnum)!=str:
            resnum = str(resnum)
        if type(data)!=float:
            data=float(data); #In case data is an integer.
        
        #Need to make sure Bfactor column is adjusted correctly.
        
        for line in self.pdb_map:
            if ((self.pdb_map[line]['residue_number']==resnum) and (self.pdb_map[line]['chain']==chain)):
                self.pdb_map[line]['b_factor']="%.2f"%data
            else:
                continue
            
        
    
    def replace_atom_b_factor(self, resnum, chain, atomname, data):
        """
        Replaces the b factor of an atom.
        Can be all string representations or not.
        """
        
        if type(resnum)!=str:
            resnum = str(resnum)
        if type(data)!=float:
            data=float(data);
        
        #Need to make sure Bfactor column is adjusted correctly.
        
        for line in self.pdb_map:
            if ((self.pdb_map[line]['residue_number']==resnum) and (self.pdb_map[line]['chain']==chain) and (self.pdb_map[line]["atom_name"]==atomname)):
                self.pdb_map[line]['b_factor']="%.2f"%data
            else:
                continue
    
    def clean_PDB(self):
        """
        Removes HSD, Waters: Tries to fix atom and residue name inconsistencies.
        HAS worked for changing a single MD pdb (NAMD) frame to Rosetta file.
        PLEASE Expand if possible to alias all residues for Rosetta compatability.
        NOT gaurenteed, but SHOULD work ok.
        """
        
        self.RESIDUES_aliased = False; self.WATER_aliased=False; self.IONS_aliased=False; self.DNA_aliased = False
        
        waters = []; #List of keys that have waters
        print "Attempting to change residue names, atom names, and water"
        for key in self.pdb_map:
            #print self.pdb_map[key]["three_letter_code"]
            def alias_dna():
                if self.pdb_map[key]["three_letter_code"]=="DA":
                    self.DNA_aliased=True
                    self.pdb_map[key]["three_letter_code"]="A"
                    
                elif self.pdb_map[key]["three_letter_code"]=="DT":
                    self.DNA_aliased=True
                    self.pdb_map[key]["three_letter_code"]="T"
                    
                elif self.pdb_map[key]["three_letter_code"]=="DC":
                    self.DNA_aliased=True
                    self.pdb_map[key]["three_letter_code"]="C"
                    
                elif self.pdb_map[key]["three_letter_code"]=="DG":
                    self.DNA_aliased=True
                    self.pdb_map[key]["three_letter_code"]="G"
                    
                else:
                    return
                
            def alias_water():
                if self.pdb_map[key]["three_letter_code"] in ["HOH", "TIP3", "WAT", "TIP5"]:
                    self.WATER_aliased=True
                    self.pdb_map[key]["three_letter_code"]="TP3"; #IO_STRING for TP3 is WAT...Buy still reads TP#?
                    self.pdb_map[key]["id"]="HETATM"
                    waters.append(key)
                    
            #def alias_ions():
                #if self.pdb_map[key]["chain"]=="I":
                    #IONS_aliased= True
                    #self.pdb_map[key]["id"]="HETATM"
                    
            def alias_residues():
                if self.pdb_map[key]["three_letter_code"] == "HSD":
                    self.RESIDUES_aliased = True
                    self.pdb_map[key]["three_letter_code"]="HIS"
                    
            def alias_atoms():
                if self.pdb_map[key]["three_letter_code"]== "SER ":
                    atom_pairs = {"  HG1":"  HG "}
                    
                elif self.pdb_map[key]["three_letter_code"]=="ILE ":
                    atom_pairs = {"  CD ":"  CD1"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="LEU ":
                    atom_pairs = {"  OT1":"  O  ", "  OT2":"  OXT"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="VAL ":
                    atom_pairs = {"  OT1":"  O  ", "  OT2":"  OXT"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="LYS ":
                    atom_pairs = {"  HZ1":"  1HZ", "  HZ2":"  2HZ", "  HZ3":"  3HZ"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="ARG ":
                    atom_pairs = {" HH11":" 1HH1", " HH12":" 2HH1", " HH21":" 1HH2", " HH22":" 2HH2"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="ASN ":
                    atom_pairs = {"HD21":"1HD2", "HD22":"2HD2"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)
                    
                elif self.pdb_map[key]["three_letter_code"]=="PRO ":
                    atom_pairs = {"  OT1":"  O  ", "  OT2":"  OXT", "  HD1":"  1HD", "  HD2":"  2HD", "  HB1":"  1HB", "  HG1":"  1HG", "  HG2":"  2HG"}
                    self.pdb_map = self.pdb_atom_alias(self.pdb_map, key, atom_pairs)    
            
            
            #Unnessessary, but organized.
            alias_water()
            #alias_ions()
            #alias_residues()
            alias_atoms()
            alias_dna()
            
        #Removes Waters. Keeps Ions.
        #for key in waters:
            #self.pdb_map.pop(key)
            
        #Outputs what was found:
        if self.RESIDUES_aliased:
            print "Residues Changed"
            
        if self.WATER_aliased:
            print "Water found...changed to TP3. Remove to decrease calculation time."
            
        if self.IONS_aliased:
            print "Ions found.  Most are able to be read into Rosetta"
        
        if self.DNA_aliased:
            print "DNA found, changed to single letter code."
        
        