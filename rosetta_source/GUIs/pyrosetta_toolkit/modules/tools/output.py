#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/output.py
## @brief  general output functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os
import re

#Tkinter Imports
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox

#Toolkit Imports
import sequence
import loops
from modules.PDB import *
from modules.tools import loops as loop_tools
from window_main import global_variables
from modules.definitions import restype_definitions

def dumpPDB(p, file, score):
    """
    Dumps the pose using the Py Job Distributor
    """
    jd=PyJobDistributor(file, 100000, score); #This number is high so that it outputs a pose even if one with the name already exists...
    #native_pose = Pose()
    #pose_from_pdb(native_pose, infile)
    #jd.native_pose=native_pose
    jd.output_decoy(p)
    print "Pose written to directory..."
    
def showPose(p, observer):
    
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    #self.obs.add_observer(p)
    score = create_score_function_ws_patch('standard', 'score12')
    print score(p)
    observer.apply(p)
    #self.obs.send_energy(p)
    #self.obs.remove_observer(p)
    #obs.pymol.apply(p)
    return p
    
def save_loop_file(p, loops_as_strings, ask_info=True):
    """
    Saves a Rosetta Loop file.  Also asks to discard residues or not.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    if not loops_as_strings:
        print "\n No loops to save...\n"
        return

    
    if ask_info:
        
        ask_cut_points = tkMessageBox.askyesno(message="Define Cutpoints?", default=tkMessageBox.NO)

        discard_loops =  tkMessageBox.askyesno(message = "Discard Phi/Psi and build using ideal bond lengths and angles?")
    else:
        ask_cut_points = False
        discard_loops = True
    outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output loop file to...")
    if not outfilename: return
    global_variables.current_directory=os.path.dirname(outfilename)
    
    FILE = open(outfilename, 'w')
    for loop_string in loops_as_strings:
        loop_stringSP = loop_string.split(":")
        start = p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[0]))
        end = p.pdb_info().pdb2pose(loop_stringSP[2], int(loop_stringSP[1]))

        #This asks the user for each loop what he/she would like to do in terms of cutpoints.  Default= somewhere near the middle.
        if ask_cut_points:
            cutpoint_known = False
            
            while not cutpoint_known:
                cut = tkSimpleDialog.askstring(title="cutpoint", prompt="Cutpoint for Loop (#, default, 0 for random) " +loop_string, initialvalue="default")
                if cut== "default":
                    cut = (end - start)/2
                    cut = start+cut
                    cutpoint_known=True
                elif cut =="0":
                    cut=0
                    cutpoint_known=True
                else:
                    if ((int(cut) < start) | (int(cut) > end)):
                        tkMessageBox.showerror(message="Invalid CutPoint!")
                        cut = int(cut)
                        cutpoint_known=False
        else:
            cut = (end - start)/2
            cut = start+cut
        FILE.write("LOOP"+" "+repr(start)+" "+repr(end)+" "+repr(cut)+" 0 "+repr(int(discard_loops))+"\n")
        
    FILE.close()
    print "\nLoop File written...\n"
    return

def clean_whitespace(obj):
    if isinstance(obj, basestring):
        return obj.strip()
    elif isinstance(obj, list):
        return [clean_whitespace(o) for o in obj]
    elif isinstance(obj, tuple):
        return tuple(clean_whitespace(o) for o in obj)
    elif isinstance(obj, dict):
        return dict((k, clean_whitespace(v)) for (k,v) in obj.items())
    else:
        return obj

def save_basic_resfile(p):
    """
    Saves an empty resfile, numbered by PDB with NATRO designation
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    ResDic = dict()
    outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output resfile to...")
    if not outfilename: return
    global_variables.current_directory=os.path.dirname(outfilename)
    
    save_resfile_w_designdic(p, ResDic, outfilename)
    
def save_resfile_w_designdic(p, ResDic, filename):
    """
    Saves a resifile, readable by PyRosetta and Rosetta.
    ResDic can be empty.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    tot = p.total_residue()
    FILE = open(filename, 'w')
    FILE.write(" start\n")
    for i in range(1, tot+1):
        chain = p.pdb_info().chain(i)
        chainStr = chain.rjust(3)
        pdbNum = p.pdb_info().number(i)
        pdbStr = str(pdbNum).rjust(4)
        res = repr(pdbNum) + ":" +chain
        if res in ResDic:
            x = ""
            print res
            for residues in ResDic[res]:
                print residues
                if residues == "NATRO":
                    line = pdbStr + chainStr + "  NATRO"+"\n"                      
                elif residues == "NATAA":
                    line = pdbStr + chainStr + "  NATAA"+ "\n"
                elif residues == "ALLAA":
                    line = pdbStr + chainStr + "  ALLAA" + "\n"
                else:
                    residuesAll = residues.split(":")
                    x = x + residuesAll[2]
                    line = pdbStr + chainStr + "  PIKAA  " + x + "\n"
            FILE.write(line)
        else:
            line = pdbStr + chainStr + "  NATRO" + "\n"
            FILE.write(line)
    FILE.close()
    print "\nRes File written....\n"  

def createSeqFile(p, newList):
    """
    Used in conversion to output scwrl sequence file.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    seq = p.sequence()
    print seq
    seq = seq.lower()
    
    #This is to seperate the string so that we can change the elements of each part.
    seqList = []
    for x in seq:
        seqList.append(x)
    
    for list in newList:
        for res in list:
            seqList[res-1]=seqList[res-1].upper()
    
    seq = ""
    for res in seqList:
        seq = seq+res
    
    print seq
    return seq
    

def saveSeqFile(p, fileout, loops_as_strings):
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    newList = loop_tools.loopArea(p, loops_as_strings)
    seq = createSeqFile(p, newList)
    FILE = open(fileout, 'w')
    FILE.write(seq+"\n")
    FILE.close()
    return



def save_basic_blueprint(p, output=True):
    """
    Saves a basic blueprint file to be manually edited.
    If output is false, returns a string of the file for manipulation.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    out_string = ""
    define = restype_definitions.definitions()
    for i in range(1, p.total_residue()+1):
        pdb_num = p.pdb_info().number(i)
        single_letter_code = define.get_one_letter_from_three(p.residue(i).name())
        out_string = out_string+repr(pdb_num)+" "+single_letter_code+" . NATRO\n"
    
    if output:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory)
        if not outfilename:return
        global_variables.current_directory=os.path.dirname(outfilename)
        
        FILE = open(outfilename, 'w')
        FILE.write(out_string)
        FILE.close()
        print "\nBlueprint Saved...\n"
    else:
        return out_string

############PDBLIST TOOLS##########################
def make_PDBLIST(directory=""):
    """
    Makes a list of PDB's from a directory.  Does not walk directory.  This can be an option later.
    Later realize could have used find command...
    """
    if directory=="":
        directory = tkFileDialog.askdirectory(initialdir = global_variables.current_directory)
        if not directory: return
        global_variables.current_directory=directory
    
    contains = tkSimpleDialog.askstring(title="Contains...", prompt="Separate mutliple matches by a coma...", initialvalue=".pdb,")
    containsSP = contains.split(",")
    FILES = os.listdir(directory)
    NEWFILE = open(directory+"/PDBLIST.txt", 'w')
    filenum=1
    for name in FILES:
        match = True; #Assumes true.  If 
        for pattern in containsSP:
            pattern = pattern.strip()
            if re.search(pattern, name) and (re.search("\._", name)== None) and (re.search("~", name)== None):#This shows how stupid python/linux can be...:
                continue
            else:
                match=False
                continue
            
        if match:
            print "File "+repr(filenum)+":"+name
            p = os.path.join(directory, name)
            filenum+=1
            NEWFILE.write(p+"\n")
                
    print "File saved as 'PDBLIST.txt' in directory specified."
    NEWFILE.close()
    return directory+"/PDBLIST.txt"

def make_PDBLIST_recursively(directory=""):
    if directory=="":
        directory = tkFileDialog.askdirectory(initialdir = global_variables.current_directory)
        if not directory: return
        global_variables.current_directory=directory
    
    contains = tkSimpleDialog.askstring(title="Contains...", prompt="Separate mutliple matches by a coma...", initialvalue=".pdb,")
    NEWFILE = open(directory+"/PDBLIST_RECURSIVE.txt", 'w')
    filenum = 1
    containsSP = contains.split(",")
    for root, dirs, files in os.walk(directory, topdown=True):
        #print "Root" + root
        for f in files:
            match = True; #Assumes true.  If 
            for pattern in containsSP:
                pattern = pattern.strip()
                if re.search(pattern, f) and (re.search("\._", f)== None) and (re.search("~", f)== None):#This shows how stupid python/linux can be...:
                    continue
                else:
                    match=False
                    continue
            
            if match:
                print "File "+repr(filenum)+":"+f
                p = os.path.join(root, f)
                filenum+=1
                NEWFILE.write(p+"\n")
    NEWFILE.close()
    print "File saved as 'PDBLIST.txt' in directory specified."
    return directory+"/PDBLIST_RECURSIVE.txt"

def return_rosetta_numbering(loops_as_strings):
    for string in loops_as_strings:
        pass
    
def convert_PDBLIST_to_sqlite3db(pdblist_path):
    """
    Adds each PDB info excluding header information into an SQLITE3 Database.  The module for this is modules/PDB.py.
    Needs to be basename as querying with '/' in a string doesn't seem to work.
    """
    
    if not pdblist_path:
        print "Please choose PDBList..."
        return
    structID = tkSimpleDialog.askstring(title="structID", prompt="These entries will have a structID of...", initialvalue="na")
    PDBLIST = open(pdblist_path, 'r')
    dbname = os.path.dirname(pdblist_path)+"/DATABASE.db"
    print dbname
    DB = PDB("", "", "", False, dbname)
    i = 1
    for filepath in PDBLIST:
        pdbID = os.path.basename(filepath).split(".")[0]
        pdbID = filepath
        DB.set_basic_options(pdbID, i, structID)
        filepath = filepath.strip()
        DB.read_pdb_into_database_flat(filepath, False, False)
        i+=1
    DB.db.close()
    print "Database written to PDBLIST directory"

def extract_pdb_from_sqlite3db():
    dbfilename = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title = "Database filename")
    if not dbfilename: return
    global_variables.current_directory = os.path.dirname(dbfilename)
    
    pdbID = tkSimpleDialog.askstring(title = "pdbID", prompt="Please enter the pdbID/filepath you wish to extract")
    if not pdbID: return
    
    pdbdb = PDB_database(dbfilename)
    print "Database Opened"
    table = pdbdb.scrub("pdb")
    pdbdb.query_pdbID(table, pdbID)
    outname = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory)
    if not outname: return
    global_variables.current_directory = os.path.dirname(outname)
    
    pdbdb.set_output_DIR(os.path.dirname(outname))
    pdbdb.save_cur_as_pdb(os.path.basename(outname))

def extract_pdbs_from_sqlite3db(pdblist_path):
    dbfilename = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title = "Database filename")
    if not dbfilename: return
    global_variables.current_directory = os.path.dirname(dbfilename)
    pdbdb = PDB_database(dbfilename)
    print "Database Opened"
    outdir = tkFileDialog.askdirectory(initialdir = global_variables.current_directory, title= "Choose output directory")
    global_variables.current_directory = outdir
    
    keep_original_filename = tkMessageBox.askyesno(title="Keep Original Filename?", message="Add numerical designation to filename to keep from overwriting same PDBs in list?", default=tkMessageBox.NO)
    strucID = tkSimpleDialog.askstring(title="strucID", prompt="Extract entries with structID of...", initialvalue="na")
    PDBLIST = open(pdblist_path, 'r')
    pdbdb.set_output_DIR(outdir)
    filenum = 1
    table = pdbdb.scrub("pdb")
    for pdbID in PDBLIST:
        pdbdb.query_pdbID_and_strucID(table, os.path.basename(pdbID), strucID)
        newname = os.path.basename(pdbID)
        if not keep_original_filename: newname = newname+"_"+repr(filenum)
        print "Saving "+ newname
        pdbdb.save_cur_as_pdb(newname)
        pdbdb._reset_cursor()
        filenum+=1
        
    print "Finished..."
    PDBLIST.close()
 
def score_PDBLIST(pdblist_path, score):
    """
    Outputs a simple pdb vs score for simple analysis.
    if pdblist_path=False, a dialog box opens.
    """
    
    if not pdblist_path:
        pdblist_path = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title = "PDBLIST")
        if not pdblist_path: return
        global_variables.current_directory = os.path.dirname(pdblist_path)
        PDBLIST = open(pdblist_path, 'r')
    else:
        PDBLIST = open(pdblist_path, 'r')
    SCORED_PDBLIST = open(os.path.dirname(pdblist_path)+"/SCORED_PDBLIST.txt", 'w')
    for path in PDBLIST:
        path = path.strip()
        print path
        p = Pose()
        try:
            pose_from_pdb(p, path)
        except PyRosettaException:
            print "Cannot Load "+path+ " Try using -ignore_unrecognized_residues in options window..."
            continue
        e = score(p)
        SCORED_PDBLIST.write(path+"\t%.3f\n"%e)
    print "\nComplete. File written to SCORED_PDBLIST.txt\n"
    PDBLIST.close()
    SCORED_PDBLIST.close()
        
def convert_PDBLIST_to_rosetta_db(current_directory):
    pass


#### FASTA OUTPUT ####

def save_FASTA(pose, base_name, outfilename = False, loops_as_strings = False ):
    """
    If outfilename is False, will ask for a directory using current_directory.
    If loops_as_strings is given, output FASTA of loops.
    Uses Pyrosetta...
    """
    if not pose.total_residue():
        print "\n No pose loaded...\n"
        return
    
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output FASTA to...")
        if not outfilename: return
        global_variables.current_directory = os.path.dirname(outfilename)
    OUTFILE = open(outfilename, 'w')
    if loops_as_strings:
        for loop_string in loops_as_strings:
            seq = sequence.get_sequence(pose, loop_string)
            header = ">"+base_name+" "+loop_string
            OUTFILE.write(header+"\n")
            OUTFILE.write(seq+"\n")
    else:
        seq = pose.sequence()
        OUTFILE.write(">"+base_name+"\n")
        OUTFILE.write(seq+"\n")
    OUTFILE.close()
    return

def save_FASTA_PDBLIST(pdblist_path, outfilename=False, loops_as_strings=False):
    """
    If outfilename is False, will ask for a filename
    Goes through each member of PDBLIST
    Uses pyrosetta, much slower...Could use my PDB class...
    """
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output FASTA to...")
        if not outfilename:return
        global_variables.current_directory = os.path.dirname(outfilename)
    
    OUTFILE = open(outfilename, 'w')
    PDBLIST = open(pdblist_path, 'r')
    for pdbpath in PDBLIST:
        pdbpath = pdbpath.strip()
        pdb = os.path.basename(pdbpath)
        pdbID = pdb.split(".")[0]
        print pdbID
        pose = Pose()
        try:
            pose_from_pdb(pose, pdbpath)
        except PyRosettaException:
            print "Could not load.. "+pdbID+"..Continueing.."
            continue
        if loops_as_strings:
            for loop_string in loops_as_strings:
                seq = sequence.get_sequence(pose, loop_string)
                header = ">"+pdbID+" "+loop_string+" "+pdbpath
                OUTFILE.write(header+"\n")
                OUTFILE.write(seq+"\n")
        else:
            seq=pose.sequence()
            OUTFILE.write(">"+pdbID+" "+pdbpath+"\n")
            OUTFILE.write(seq+"\n")
    PDBLIST.close()
    OUTFILE.close()
    print "Fasta written..."
    return

    
def exportPDBSCORE(self, score):
    """
    Exports a list of scores.
    """
    
    PDBLIST = tkFileDialog.askopenfile(title = "PDBLIST", initialdir = global_variables.current_directory)
    OUTFILE = tkFileDialog.asksaveasfile(title = "Save As...", initialdir = global_variables.current_directory)
    
    if PDBLIST == None or OUTFILE==None:
        return
    
    for PDBPath in PDBLIST:
        PDBPath = PDBPath.strip()
        
        p = Pose()
        pose_from_pdb(p, PDBPath)
        SCORE = score(p)
        print SCORE
        OUTFILE.write(PDBPath+":%.3f"%SCORE+"\n")
    print "\nComplete\n"
    PDBLIST.close()
    OUTFILE.close()