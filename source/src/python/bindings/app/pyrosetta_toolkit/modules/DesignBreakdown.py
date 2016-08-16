#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/DesignBreakdown.py
## @brief  Class for analyzing design results from a fasta of sequences
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import re
import os
import sqlite3
import sys
from optparse import OptionParser, IndentedHelpFormatter

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkSimpleDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.prettytable.prettytable import *
from app.pyrosetta_toolkit.modules.definitions.restype_definitions import definitions

class DesignBreakdown:
    """
    This class functions in organizing results for a Rosetta design run.  Perhaps you can do this using features.  Probably.  Anyhow, this will work for the GUI.
    It can output a text file, a database, and run an R script to graph the results.
    It also has some GUI functions, but can be run as an independant script.
    FASTA should have > format: pdbID region pdbpath OR pdbID pdbpath for whole structure comparisons. (region designated as start:end:chain)
    LIMITATIONS:
       1) Currently, it does not deal with extensions or deletions in the designed poses.
       2) Currently, No DNA, polymers, or other non-cannonicals.
       3) Currently, Does not print or output decoys with mutations of x, y, z at positions a, b, c
          However, it does output a raw_data table in the SQLITE3 database, which will allow advanced querys to get this answer.
    """

    def __init__(self, fasta_path, reference_path, output_directory=False, region=False):

        self.output_directory = output_directory
        self.sequences = []; # List of SequenceInfo objects
        self.reference_path = reference_path


        self.reference_sequence = None
        self.reference_pose = Pose()
        self.main_region = region
        self.regions = dict()

        self.load_sequences(fasta_path)
        self.results = SequenceResults()
        self.aa_codes = definitions().get_all_one_letter_codes()

        #Calculate - If return number is 0, exit.
        if not self.calculate_results():
            return

        #Data Output
        if not self.output_directory:
            self.output_directory = os.path.dirname(fasta_path)+'/'+os.path.basename(fasta_path)+"_RESULTS"
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)
            print "Outputting results to: "+self.output_directory



    def run_outputs(self):
        self.output_basic_table()
        self.output_prettytable()
        self.output_database()
        self.output_plots()

    def load_sequences(self, fasta_path):
        """
        Opens Fasta file and appends a list of SequenceInfo objects.
        """

        FILE = open(fasta_path, 'r')
        for line in FILE:
            if re.search(">", line):
                info = line.strip()
                info = info[1:]
                infoSP = info.split()
                Seq = SequenceInfo()
                if len(infoSP)<2:continue
                elif len(infoSP)==2:
                    Seq.set_pdbID(infoSP[0])
                    Seq.set_pdbpath(infoSP[1])

                elif len(infoSP)==3:
                    Seq.set_pdbID(infoSP[0])
                    Seq.set_region(infoSP[1])
                    self.regions[infoSP[1]]=""
                    Seq.set_pdbpath(infoSP[2])
                self.sequences.append(Seq)
                continue
            line = line.strip()
            if not line:continue
            print line
            self.sequences[-1].set_sequence(line)
        FILE.close()
        print "Sequences Loaded"

    def calculate_results(self):
        print "Calculating Results"

        pose_from_file(self.reference_pose, self.reference_path)
        print self.reference_pose
        self.reference_sequence = self.reference_pose.sequence()
        print self.reference_sequence
        for i in range(0, self.reference_pose.total_residue()):
            self.results.add_reference_residue(i+1, self.reference_sequence[i])
        if not self.sequences:return 0

        #Check to make sure reference pose matches length of all sequences
        #Check that all sequences are the same length
        if not self.regions:
            if not self.are_sequences_same_length:
                print "Sequences are not same length, and no region was specified in fasta.  Cannot continue"
                return 0
            if (self.reference_pose.total_residue() != self.sequences[0].get_length()):
                print "Sequence length of Fasta does not match sequence length of reference pose."
                region = self.main_region
                if not region:
                    region = tkSimpleDialog.askstring(title="Region", prompt ="Please enter a region: start end chain for PDB numbering or start end for Rosetta numbering")
                if not region:return 0
                regionSP = region.split()

                #Make sure it's entered correctly.  If not, try one more time before returning.
                if 1<=len(regionSP)>3:
                    print "Please re-enter region."
                    region = tkSimpleDialog.askstring(title="Region", prompt ="Please enter a region: start end chain for PDB numbering or start end for Rosetta numbering")
                    if not region:return
                regionSP = region.split()
                if 1<=len(regionSP)>3:print "Region not recognized.  Returning."; return
                if len(regionSP)==3:
                    self.regions[":".join(regionSP)]=""
                    for Seq in self.sequences:
                        Seq.set_region(":".join(regionSP))
                elif len(regionSP)==2:
                    if self.reference_pose.pdb_info().pose2pdb(int(regionSP[0])).split()[1]== self.reference_pose.pdb_info().pose2pdb(int(regionSP[1])).split()[1]:
                        print "One chain specified."
                        chain = self.reference_pose.pdb_info().pose2pdb(int(regionSP[0])).split()[1]
                        pdb_start = self.reference_pose.pdb_info().pose2pdb(int(regionSP[0])).split()[0]
                        pdb_end = self.reference_pose.pdb_info().pose2pdb(int(regionSP[0])).split()[0]
                        self.regions[":".join(regionSP)]=""
                        for Seq in self.sequences:
                            Seq.set_region(":".join(regionSP))
                    else:
                        print "Multiple chains in region found. Splitting sequences to match regions."
                        self.split_region_and_fix_sequences(int(regionSP[0]), int(regionSP[1]))




        #Calculate Results

        if not self.regions:
            l = len(self.sequences[0].get_sequence())
            for Seq in self.sequences:
                for i in range(1, l+2):
                    residue = Seq.get_residue(i)
                    self.results.add_residue(i, residue, Seq.get_pdbpath())
        else:
            if not self.are_sequences_for_regions_same_length():
                return 0

            for region in self.regions:
                print region
                regionSP = region.split(":")
                for Seq in self.sequences:
                    #print Seq.get_sequence()
                    if Seq.get_region()==region:
                        #This is so that if there are missing numbers in the PDB between start:end in the region:
                        start = self.reference_pose.pdb_info().pdb2pose(Seq.get_chain(), Seq.get_start_residue())
                        for i in range(0, Seq.get_length()):
                            print i
                            self.results.add_residue(start+i, Seq.get_residue(Seq.get_start_residue()+i), Seq.get_pdbpath())
        return 1

    def output_basic_table(self):
        """
        Outputs a basic table of all data for importing into Excel, R, or other script.  Tab delimited.
        """
        resnums = self.results.get_all_residue_numbers()
        reference_line = "#\t"
        resnum_line = "\t"
        conserved_line = "#\t"
        OUTFILE = open(self.output_directory+"/RAW_DESIGN_TABLE.txt", 'w')
        OUTFILE.write("# TOTAL_SEQUENCES "+repr(len(self.sequences))+"\n")
        for num in resnums:
            pdb_num = self.reference_pose.pdb_info().pose2pdb(num)
            SP = pdb_num.split()
            pdb_num = SP[0]+SP[1]; #Now it it resnumchain like 10A 11B etc.
            resnum_line = resnum_line+pdb_num+"\t"
            if self.reference_sequence:
                reference_line = reference_line+self.results.get_reference_residue(num)+"\t"
                conserved_line = conserved_line+self.results.get_percent_string(num, self.results.get_reference_residue(num))+"\t"
        if self.reference_sequence:
            OUTFILE.write(reference_line+"\n")
            OUTFILE.write(conserved_line+"\n")
        OUTFILE.write(resnum_line+"\n")
        for aa in self.aa_codes:
            line = aa+"\t"
            for num in resnums:
                line=line+self.results.get_percent_string(num, aa)+"\t"
            OUTFILE.write(line+"\n")
        print "Raw file written to RAW_DESIGN_TABLE.txt"
        OUTFILE.close()

    def output_prettytable(self):
        OUTFILE = open(self.output_directory+"/PRETTY_DESIGN_TABLE.txt", 'w')
        OUTFILE.write("# TOTAL_SEQUENCES "+repr(len(self.sequences))+"\n")
        resnums = self.results.get_all_residue_numbers()
        main_row = ["residue"]
        conserved_row = ["conserved"]
        reference_row = ["reference"]
        if not self.regions:
            for num in resnums:
                main_row.append(self.get_correct_pdb_number_string(num))
                if self.reference_sequence:
                    reference_row.append(self.results.get_reference_residue(num))
                    conserved_row.append(self.results.get_percent_string(num, self.results.get_reference_residue(num)))
            table = PrettyTable(main_row)
            if self.reference_sequence:
                table.add_row(reference_row)
                table.add_row(conserved_row)
            for aa in self.aa_codes:
                row = [aa]
                for num in resnums:
                    row.append(self.results.get_percent_string(num, aa))
                table.add_row(row)
            out_string = table.get_string()
            OUTFILE.write(out_string)
        else:
            for region in self.regions:
                OUTFILE.write('# REGION '+region+"\n")
                for num in resnums:

                    if not self.check_if_rosetta_resnum_is_part_of_region(num, region):
                        continue

                    main_row.append(self.get_correct_pdb_number_string(num))
                    if self.reference_sequence:
                        reference_row.append(self.results.get_reference_residue(num))
                        conserved_row.append(self.results.get_percent_string(num, self.results.get_reference_residue(num)))
                table = PrettyTable(main_row)
                if self.reference_sequence:
                    table.add_row(reference_row)
                    table.add_row(conserved_row)
                for aa in self.aa_codes:
                    row = [aa]
                    for num in resnums:
                        if not self.check_if_rosetta_resnum_is_part_of_region(num, region):
                            continue
                        row.append(self.results.get_percent_string(num, aa))
                    table.add_row(row)

                out_string = table.get_string()
                OUTFILE.write(out_string)
        print "PrettyTable file written to PRETTY_DESIGN_TABLE.txt"
        OUTFILE.close()

    def output_database(self):
        self.db_path = self.output_directory+"/SQL_DESIGN_TABLE.db"
        db = sqlite3.connect(self.db_path)
        cur = db.cursor()
        resnums = self.results.get_all_residue_numbers()

        with db:
            #Hard to look at, easy to execute queries on data:
            #Like this: select * from design_data where prob>=.3 and type='design'.  Awesomeness.
            cur.execute("CREATE TABLE IF NOT EXISTS design_data(id integer PRIMARY KEY, region TEXT, type TEXT, pdb_position TEXT, rosetta_position INT, aa TEXT, prob REAL, freq INT, total_sequences INT, ref_name TEXT, decoys TEXT)")

            if not self.regions:
                i=0
                for num in resnums:
                    main_row.append(self.get_correct_pdb_number_string(num))
                    i+=1
                    if self.reference_sequence:
                        cur.execute("INSERT INTO design_data VALUES(NULL, ?,?,?,?,?,?,?,?,?,?)", \
                                ("full", "reference", self.get_correct_pdb_number_string(num), num, self.results.get_reference_residue(num), 1.00, 1, self.get_total_sequences(), self.reference_pose.pdb_info().name(), self.reference_pose.pdb_info().name()))

                    for aa in self.aa_codes:
                        i+=1
                        cur.execute("INSERT INTO design_data VALUES(NULL, ?,?,?,?,?,?,?,?,?,?)", \
                                ("full", "design", self.get_correct_pdb_number_string(num), num, aa, self.results.get_percent(num, aa),self.results.get_freq(num, aa), \
                                 self.get_total_sequences(), self.reference_pose.pdb_info().name(), ":".join(self.results.get_decoys_with_aa(num, aa))))


            else:
                i = 0
                for region in self.regions:
                    for num in resnums:
                        i+=1
                        if not self.check_if_rosetta_resnum_is_part_of_region(num, region):
                            continue

                        cur.execute("INSERT INTO design_data VALUES(NULL, ?,?,?,?,?,?,?,?,?,?)", \
                                (region, "reference", self.get_correct_pdb_number_string(num), num, self.results.get_reference_residue(num), 1.00, 1, self.get_total_sequences(region), self.reference_pose.pdb_info().name(), self.reference_pose.pdb_info().name()))

                        for aa in self.aa_codes:
                            i+=1
                            cur.execute("INSERT INTO design_data VALUES(NULL, ?,?,?,?,?,?,?,?,?,?)", \
                                (region, "design", self.get_correct_pdb_number_string(num), num, aa, self.results.get_percent(num, aa),self.results.get_freq(num, aa), \
                                 self.get_total_sequences(), self.reference_pose.pdb_info().name(),":".join(self.results.get_decoys_with_aa(num, aa))))

        #Raw data table
        #So you can query for combinations and get decoys with specific mutations at positions and above a probablity.
            cur.execute("create table if not exists raw_data(id integer PRIMARY KEY, pdb_position TEXT, rosetta_position INT, aa TEXT, decoy TEXT)")
            if not self.regions:
                l = len(self.sequences[0].get_sequence())
                x=0
                for Seq in self.sequences:
                    for i in range(1, l+2):
                        x+=1
                        residue = Seq.get_residue(i)
                        cur.execute("INSERT INTO raw_data VALUES(NULL, ?,?,?,?)", \
                            (self.get_correct_pdb_number_string(i), i, residue, Seq.get_pdbpath()))

            else:
                x=0
                for region in self.regions:
                    regionSP = region.split(":")
                    for Seq in self.sequences:
                        if Seq.get_region()==region:
                            #This is so that if there are missing numbers in the PDB between start:end in the region:
                            start = self.reference_pose.pdb_info().pdb2pose(Seq.get_chain(), Seq.get_start_residue())
                            for i in range(0, Seq.get_length()):
                                x+=1
                                num = start+i
                                cur.execute("INSERT INTO raw_data VALUES(NULL, ?,?,?,?)", \
                                    (self.get_correct_pdb_number_string(num), num, Seq.get_residue(Seq.get_start_residue()+i), Seq.get_pdbpath()))

        print "Database written to SQL_DESIGN_TABLE.db"

    def output_plots(self):
        script = self.location()+"/R_Scripts/DesignBreakdown.R"
        os.system("Rscript "+script+' '+self.db_path+' '+self.output_directory)
        print "Plots written to PLOTS.pdb for each region."

### Helper Functions ###
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported/ran from anywhere.
        """

        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP[0]

    def are_sequences_same_length(self):
        """
        Determine if all items of the sequence list are the same number of residues
        """
        return all(x.get_length() == self.sequences[0].get_length() for x in self.sequences)

    def are_sequences_for_regions_same_length(self):
        """
        Assertion that sequences are the same length given their region.
        """
        #Setup dictionary for checking
        region_Seq_map = dict()
        for region in self.regions:
            if not region_Seq_map.has_key(region):
                region_Seq_map[region]=[]
        for Seq in self.sequences:
            region_Seq_map[Seq.get_region()].append(Seq)

        #Check length for each region in dictionary.
        same_length=True
        for region in self.regions:
            #for x in region_Seq_map[region]:
                #print x.get_sequence()
            if not all(x.get_length() == region_Seq_map[region][0].get_length() for x in region_Seq_map[region]):
                print "Sequences for region "+region+" are not the same length."
                same_length=False

        return same_length

    def check_if_rosetta_resnum_is_part_of_region(self, resnum, region):
        region_start = region.split(":")[0]
        region_end = region.split(":")[1]
        region_chain = region.split(":")[2]
        pdb_num = self.reference_pose.pdb_info().pose2pdb(resnum)
        SP = pdb_num.split()
        pdb_num = SP[0]; pdb_chain = SP[1]

        if (region_start <=pdb_num<=region_end) and pdb_chain==region_chain:
            return True
        else:
            return False

    def get_correct_pdb_number_string(self, resnum):
        """
        Gets PDB numbering from pose numbering and switches order of chain and num.  chain_num->num_chain
        """
        pdb_num = self.reference_pose.pdb_info().pose2pdb(resnum)
        SP = pdb_num.split()
        pdb_num_string = SP[0]+SP[1]; #Now it it resnumchain like 10A 11B etc.
        return pdb_num_string

    def split_region_and_fix_sequences(self, start, end):
        """
        Splits a Rosetta numbered region into PDB regions for each chain.  Adds regions to self.regions, splits Seq in self.sequences.
        """
        pass

    def get_total_sequences(self, region=False):
        if not region:
            return len(self.sequences)
        else:
            l = 0
            for Seq in self.sequences:
                if Seq.get_region==region:
                    l+=1
            return l

class SequenceResults:
    """
    Simple class for holding, calculating, + accessing result data
    Residue Numbers are in Rosetta numbering.
    """
    def __init__(self):
        self.data = dict()
        self.reference = dict()

    def add_residue(self, resnum, one_letter_code, decoy):
        if not self.data.has_key(resnum):
            self.data[resnum]=dict()
            self.data[resnum][one_letter_code]=dict()
            self.data[resnum][one_letter_code]['freq']=1
            self.data[resnum][one_letter_code]['decoys']=[]
            self.data[resnum][one_letter_code]['decoys'].append(decoy);#This is to keep track of which decoys have which mutation.

        else:
            if not self.data[resnum].has_key(one_letter_code):
                self.data[resnum][one_letter_code]=dict()
                self.data[resnum][one_letter_code]['freq']=0
                self.data[resnum][one_letter_code]['decoys']=[]
            self.data[resnum][one_letter_code]['freq']+=1
            self.data[resnum][one_letter_code]['decoys'].append(decoy)

    def add_reference_residue(self, resnum, one_letter_code):
        self.reference[resnum]=one_letter_code

    def get_freq(self, resnum, one_letter_code):
        try:
            x = self.data[resnum][one_letter_code]['freq']
            return x
        except KeyError:
            return 0

    def get_total(self, resnum):
        total = 0
        for code in self.data[resnum]:
            freq = self.get_freq(resnum, code)
            total = total +freq
        return total

    def get_percent(self, resnum, one_letter_code):
        total = self.get_total(resnum)
        freq = self.get_freq(resnum, one_letter_code)

        percent = float(freq)/float(total)
        return percent

    def get_percent_string(self, resnum, one_letter_code):
        return "%.2f"%self.get_percent(resnum, one_letter_code)

    def get_reference_residue(self, resnum):
        return self.reference[resnum]

    def get_all_residues_observed(self, resnum):
        return sorted(self.data[resnum].keys())

    def get_all_residue_numbers(self):
        return sorted(self.data.keys())

    def get_decoys_with_aa(self, resnum, one_letter_code):
        """
        Returns all decoys with a specific mutation at a position.
        """
        try:
            return self.data[resnum][one_letter_code]['decoys']
        except KeyError:
            return []

    def get_decoys_with_joint_aa(self, resnum_one_letter_code_pair):
        """
        Will output decoys that have x, y, z mutations at positions a, b, c
        """
        pass

    ### reference Comparison Functions ###
    def get_all_mutated_positions(self):
        mutated_positions = []
        for resnum in self.data:
            if not self.reference.has_key(resnum):
                print "Position in data does not match position in reference"

            if self.get_percent(resnum, self.reference[resnum])==1.0:
                pass
            else:
                mutated_positions.append(resnum)

        if mutated_positions:return mutated_positions
        else:print "No mutations found"

    def get_all_reference_percent_observed(self):
        """
        Returns array of tripplets of [postion, one_letter_code, percent] of reference amino acid found.
        """
        tripplet_array = []
        for resnum in self.reference:
            if not self.data.has_key(resnum):
                print "Position in reference does not match any position in data"

            percent = self.get_percent(resnum, self.reference[resnum])
            tripplet = [resnum, self.reference[resnum], percent]
            tripplet_array.append(tripplet)
        return tripplet_array

class SequenceInfo:
    """
    Simple class for holding + accessing sequence metadata
    """

    def __init__(self):
        self.start = ""

    def get_sequence(self):
        return self.sequence
    def get_length(self):
        return len(self.sequence)
    def get_pdbID(self):
        return self.pdbID
    def get_pdbpath(self):
        return self.pdbpath
    def get_region(self):
        return self.region
    def get_start_residue(self):
        return self.start
    def get_end_residue(self):
        return self.end
    def get_chain(self):
        return self.chain

    def get_residue(self, resnum):
        """
        If region is given, resnum is residue number of PDB
        If not, resnum in Rosetta resnum
        """

        if self.start:
            index_num = resnum-self.start
            print index_num
            one_letter_code = self.sequence[index_num]
            return one_letter_code
        else:
            #Rosetta Resnum
            one_letter_code = self.sequence[int(resnum)-1]
            return one_letter_code

    def set_sequence(self, sequence):
        print sequence
        self.sequence = sequence
    def set_pdbID(self, pdbID):
        self.pdbID = pdbID
    def set_pdbpath(self, pdbpath):
        self.pdbpath = pdbpath
    def set_region(self, region):
        self.region = region
        rSP = region.split(":")
        self.start = int(rSP[0])
        self.end = int(rSP[1])
        self.chain = rSP[2]

if __name__ == '__main__':
    """
    For testing and use outside of GUI.

    """
    Tk()
    rosetta.init()

    parser = OptionParser()
    args = sys.argv
    parser.add_option("--fasta", "-f",
        help = "Path to FASTA file for design comparisons"
    )
    parser.add_option("--outpath","-o",
        default="/RESULTS",
        help = "Full output directory path.  Default is fasta file path /Results"
    )

    parser.add_option("--reference", "-r",
        default=None,
        help = "Reference pose for numbering and comparison (Required)"
    )

    parser.add_option("--region", "-g",
        default=None,
        help = "Region  - if none is given in Fasta + Not whole structure used for comparison (region designated as start:end:chain)"
    )

    (options, args) = parser.parse_args(args=args[1:])
    if not options.fasta or not os.path.exists(options.fasta):
        sys.exit("Please specify a FASTA file to use for calculations.")
    if options.outpath == "/RESULTS":
        options.outpath = os.path.dirname(options.fasta)+"/RESULTS"

    if not options.reference:
        sys.exit("Reference pdb required.")


    breakdown = DesignBreakdown(options.fasta, options.reference, options.outpath, options.region)
    breakdown.run_outputs()
