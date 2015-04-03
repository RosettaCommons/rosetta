#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/definitions/restype_definitions.py
## @brief  Container class for restypes
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#This object holds data used to define certain things.  CDR's, SASA, Residue types and names, etc.  Hopefully, this will streamline working within python.


#This seems like an ok way to do this.  Keep.
#All the dictionaries for residue Types in design

class definitions():
    def __init__(self):
        self.restypes = ("All", "Charged", "Positive", "Negative", "Non-Polar", "Polar-Uncharged", "Polar", "Aromatic", "Hydroxyl", "Conserved", "etc")
        self.restype_info= dict()
        self.resinfo = dict()
        self.set_mutation_info()
        self.set_residue_info()

    def get_residue_info(self):
        return self.restype_info
    
    def get_mutation_info(self):
        return self.resinfo
    
    def get_one_letter_from_three(self, three_letter_code):
        three_letter_code = three_letter_code.split("_")[0].upper(); #Fix for Rosetta designated chain endings.
        if three_letter_code=="CYD": three_letter_code="CYS"; #Fix for Disulfide
        for triplet in self.restype_info["All"]:
            tripletSP = triplet.split(":")
            if tripletSP[1]==three_letter_code:
                return tripletSP[2]
    
    def get_three_letter_from_one(self, one_letter_code):
        one_letter_code = one_letter_code.upper()
        for triplet in self.restype_info["All"]:
            tripletSP = triplet.split(":")
            if tripletSP[2]==one_letter_code:
                return tripletSP[1]
    
    def get_all_one_letter_codes(self):
        one_letter_codes = []
        for triplet in self.restype_info["All"]:
            tripletSP = triplet.split(":")
            one_letter_codes.append(tripletSP[2])
        return sorted(one_letter_codes)
    
    def set_residue_info(self):
        self.restype_info["Charged"]=("Lysine:LYS:K",
                             "Arginine:ARG:R"
                             "Glutamate:GLU:E",
                             "Aspartate:ASP:D",
                             "Histidine:HIS:H")
        
        self.restype_info["Positive"]=("Lysine:LYS:K",
                              "Arginine:ARG:R",
                              "Histidine:HIS:H")
        
        self.restype_info["Negative"]=("Glutamate:GLU:E",
                              "Aspartate:ASP:D")
        
        self.restype_info["Non-Polar"]=("Glycine:GLY:G",
                               "Alanine:ALA:A",
                               "Valine:VAL:V",
                               "Leucine:LEU:L",
                               "Isoleucine:ILE:I",
                               "Methionine:MET:M",
                               "Proline:PRO:P")
        
        self.restype_info["Polar-Uncharged"]=("Serine:SER:S",
                                     "Threonine:THR:T",
                                     "Cysteine:CYS:C",
                                     "Glutamine:GLN:Q",
                                     "Asparagine:ASN:N")
        
        self.restype_info["Polar"] = ("Serine:SER:S",
                             "Threonine:THR:T",
                             "Cysteine:CYS:C",
                             "Glutamine:GLN:Q",
                             "Asparagine:ASN:N",
                             "Lysine:LYS:K",
                             "Arginine:ARG:R",
                             "Glutamate:GLU:E",
                             "Aspartate:ASP:D",
                             "Histidine:HIS:H")
        
        self.restype_info["Aromatic"] = ("Phenylalanine:PHE:F",
                                "Tyrosine:TYR:Y",
                                "Tryptophan:TRP:W")
        
        self.restype_info["Hydroxyl"] = ("Serine:SER:S",
                                "Threonine:THR:T",
                                "Tyrosine:TYR:Y")
        
        self.restype_info["All"] = ("Alanine:ALA:A",
                           "Arginine:ARG:R",
                           "Asparagine:ASN:N",
                           "Aspartate:ASP:D",
                           "Cysteine:CYS:C",
                           "Glutamate:GLU:E",
                           "Glutamine:GLN:Q",
                           "Glycine:GLY:G",
                           "Histidine:HIS:H",
                           "Leucine:LEU:L",
                           "Isoleucine:ILE:I",
                           "Lysine:LYS:K",
                           "Methionine:MET:M",
                           "Phenylalanine:PHE:F",
                           "Proline:PRO:P",
                           "Serine:SER:S",
                           "Threonine:THR:T",
                           "Tryptophan:TRP:W",
                           "Tyrosine:TYR:Y",
                           "Valine:VAL:V")
        
        self.restype_info["Conserved"] = ("All Conserved Mutations",
                                 "All Conserved Mutations+Self")
        
        self.restype_info["Conserved:ALA"]=("Serine:SER:S",
                                   "Glycine:GLY:G",
                                   "Threonine:THR:T",
                                   "Proline:PRO:P")
        
        self.restype_info["Conserved:ARG"]=("Lysine:LYS:K",
                                   "Histidine:HIS:H",
                                   "Tryptophan:TRP:W",
                                   "Glutamine:GLN:Q")
        self.restype_info["Conserved:ASN"]=("Apartate:ASP:D",
                                   "Serine:SER:S",
                                   "Histidine:HIS:H",
                                   "Lysine:LYS:K",
                                   "Glutamine:GLN:Q",
                                   "Glutamate:GLU:E")
        self.restype_info["Conserved:ASP"]=("Asparagine:ASN:N",
                                   "Glutamate:GLU:E",
                                   "Glutamine:GLN:Q",
                                   "Histidine:HIS:H",
                                   "Glycine:GLY:G")
        self.restype_info["Conserved:CYS"]=("Serine:SER:S")
        self.restype_info["Conserved:GLY"]=("Alanine:ALA:A",
                                   "Serine:SER:S")
        self.restype_info["Conserved:GLN"]=("Glutamate:GLU:E",
                                   "Asparagine:ASN:N",
                                   "Glutamine:GLN:Q",
                                   "Histidine:HIS:H",
                                   "Glycine:GLY:G")
        self.restype_info["Conserved:GLU"]=("Asparagine:ASN:N",
                                   "Glutamine:GLN:Q",
                                   "Asparagine:ASN:N",
                                   "Histidine:HIS:H")
        self.restype_info["Conserved:HIS"]=("Asparagine:ASN:N",
                                   "Glutamine:GLN:Q",
                                   "Aspartate:ASP:D",
                                   "Glutamate:GLU:E",
                                   "Arginine:ARG:R")
        self.restype_info["Conserved:ILE"]=("Valine:VAL:V",
                                   "Leucine:LEU:L",
                                   "Methionine:MET:M")
        self.restype_info["Conserved:LEU"]=("Methionine:MET:M",
                                   "Isoleucine:ILE:I",
                                   "Valine:VAL:V",
                                   "Phenylalanine:PHE:F")
        self.restype_info["Conserved:LYS"]=("Arginine:ARG:R",
                                   "Glutamine:GLN:Q",
                                   "Asparagine:ASN:N")
        self.restype_info["Conserved:MET"]=("Leucine:LEU:L",
                                   "Isoleucine:ILE:I",
                                   "Valine:VAL:V",
                                   "Phenylalanine:PHE:F")
        self.restype_info["Conserved:PHE"]=("Tyrosine:TYR:Y",
                                   "Leucine:LEU:L",
                                   "Isoleucine:ILE:I")
        self.restype_info["Conserved:PRO"]=("Alanine:ALA:A",
                                   "Serine:SER:S")
        self.restype_info["Conserved:SER"]=("Alanine:ALA:A",
                                   "Threonine:THR:T",
                                   "Asparagine:ASN:N",
                                   "Glycine:GLY:G",
                                   "Proline:PRO:P")
        self.restype_info["Conserved:THR"]=("Serine:SER:S",
                                   "Alanine:ALA:A",
                                   "Valine:VAL:V")
        self.restype_info["Conserved:TYR"]=("Phenylalanine:PHE:F",
                                   "Histidine:HIS:H",
                                   "Tryptophan:TRP:W")
        self.restype_info["Conserved:TRP"]=("Phenylalanine:PHE:F",
                                   "Tyrosine:TYR:Y")
        self.restype_info["Conserved:VAL"]=("Isoleucine:ILE:I",
                                   "Leucine:LEU:L",
                                   "Methionine:MET:M")
        self.restype_info["etc"] = ("NATRO", "NATAA", "ALLAA")
        
    def set_mutation_info(self):
        #self.resinfo = SA : Rel Mut : Surf Prob
        self.resinfo = dict()
        self.resinfo["ALA"]=("115", "100", "62")
        self.resinfo["ARG"]=("225", "65", "99")
        self.resinfo["ASN"]=("160", "134", "88")
        self.resinfo["ASP"]=("150", "106", "85")
        self.resinfo["CYS"]=("150", "106", "85")
        self.resinfo["ASP"]=("135", "20", "55")
        self.resinfo["GLU"]=("190", "102", "82")
        self.resinfo["GLN"]=("180", "93", "93")
        self.resinfo["GLY"]=("75", "49", "64")
        self.resinfo["HIS"]=("195", "66", "83")
        self.resinfo["ILE"]=("175", "96", "40")
        self.resinfo["LEU"]=("170", "40", "55")
        self.resinfo["LYS"]=("200", "56", "97")
        self.resinfo["MET"]=("185", "94", "60")
        self.resinfo["PHE"]=("210", "41", "50")
        self.resinfo["PRO"]=("145", "56", "82")
        self.resinfo["SER"]=("115", "120", "78")
        self.resinfo["THR"]=("140", "97", "77")
        self.resinfo["TRP"]=("225", "18", "73")
        self.resinfo["TYR"]=("230", "41", "85")
        self.resinfo["VAL"]=("155", "74", "46")
        self.resinfo["HIS_D"]=self.resinfo["HIS"]
    
