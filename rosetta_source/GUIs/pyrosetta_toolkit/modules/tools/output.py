#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/output.py
## @brief  general output functions for the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox
from modules.PDB import *
from modules.tools import loops as loop_tools
import sequence
import os
import loops
import re




def dumpPDB(p, file, score):
    """
    Dumps the pose using the Job Distributor
    """
    jd=PyJobDistributor(file, 100000, score); #This number is high so that it outputs a pose even if one with the name already exists...
    #native_pose = Pose()
    #pose_from_pdb(native_pose, infile)
    #jd.native_pose=native_pose
    jd.output_decoy(p)
    print "Pose written to directory..."
    
def showPose(p, observer):
    
    #self.obs.add_observer(p)
    score = create_score_function_ws_patch('standard', 'score12')
    print score(p)
    observer.apply(p)
    #self.obs.send_energy(p)
    #self.obs.remove_observer(p)
    #obs.pymol.apply(p)
    return p
    
def savLoop(p, out, loops_as_strings, ask_info=True):
    """
    Saves a Rosetta Loop file.  Also asks to discard residues or not. Should be rewritten.
    """
    newList = loop_tools.loopArea(p, loops_as_strings)
    if ask_info:
        
        defcut = tkMessageBox.askquestion(message="Define Cutpoints?", default=tkMessageBox.NO)
        discard = tkMessageBox.askquestion(message = "Discard Phi/Psi and build using ideal bond lengths and angles?")
    else:
        defcut = True
        discard = 0
        
    if defcut =="yes": defcut = True
    if discard=='yes': discard= 1
    else: discard=0
    
    FILE = open(out, 'w')
    l = 0;
    for loop in newList:
        start = loop[0]
        end = loop[-1]
        print "Got start and end."
        cut = (end - start)/2
        cut = start+cut
        #This asks the user for each loop what he/she would like to do in terms of cutpoints.  Default= somewhere near the middle.
        if defcut ==True:
            cut = tkSimpleDialog.askstring(title="cutpoint", prompt="Cutpoint for Loop (#, default, 0 for random)" + loops_as_strings[l], initialvalue="default")
            if cut== "default":
                cut = (end - start)/2
                cut = start+cut
            elif cut =="0":
                cut=0
            else:
                LisSp = loops_as_strings[l].split(":")
                startSP = LisSp[0]
                cut = int(cut) - int(startSP)
                cut = start + cut
                if ((cut <= start) | (cut >= end)):
                    #Future - makes sure cut point exists, and then have them point it back in.
                    tkMessageBox.showerror(message="Invalid CutPoint!")
                    savLoop(p, out, loops_as_strings)
                    return
                
        FILE.write("LOOP"+" "+repr(start)+" "+repr(end)+" "+repr(cut)+" 0 "+repr(discard)+"\n")
        l +=1
        
    FILE.close()
    print "Loop File written..."
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
        
def save_resfile_w_designdic(p, ResDic, filename):
    """
    Saves a Design Residue file, readable by PyRosetta and Rosetta.
    """
    
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
    print "Res File written...."  

def createSeqFile(p, newList):
    """
    Used in conversion to output scwrl sequence file.
    """
    
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
    

def saveSeqFile(p, fileout, loopsLis):
    newList = loop_tools.loopArea(p, loopsLis)
    seq = createSeqFile(p, newList)
    FILE = open(fileout, 'w')
    FILE.write(seq+"\n")
    FILE.close()
    return

def make_PDBLIST(current_directory, directory=""):
    """
    Makes a list of PDB's from a directory.  Does not walk directory.  This can be an option later.
    Later realize could have used find command...
    """
    if directory=="":
        directory = tkFileDialog.askdirectory(initialdir = current_directory)
    if directory ==None:
        return
    FILES = os.listdir(directory)
    NEWFILE = open(directory+"/PDBLIST.txt", 'w')
    for names in FILES:
        if re.search(".pdb", names) and (re.search("\._", names)== None) and (re.search("~", names)== None):#This shows how stupid python/linux can be...:
            print names
            p = os.path.join(directory, names)
            NEWFILE.write(p+"\n")
    print "Written..."
    print "File saved as 'PDBLIST.txt' in directory specified."
    NEWFILE.close()
    return directory+"/PDBLIST.txt"

def make_PDBLIST_recursively(current_directory, directory=""):
    if directory=="":
        directory = tkFileDialog.askdirectory(initialdir = current_directory)
    if directory ==None:
        return
    contains = ".pdb"
    NEWFILE = open(directory+"/PDBLIST_RECURSIVE.txt", 'w')
    filenum = 1
    for root, dirs, files in os.walk(directory, topdown=True):
        #print "Root" + root
        for f in files:
            if re.search(".pdb", f) and (re.search("\._", f)== None) and (re.search("~", f)== None):#This shows how stupid python/linux can be...:
                print "File_"+repr(filenum)+"_"+f
                p = os.path.join(root, f)
                filenum+=1
                NEWFILE.write(p+"\n")
                filenum+=1
    NEWFILE.close()
    print "File saved as 'PDBLIST.txt' in directory specified."
    return directory+"/PDBLIST_RECURSIVE.txt"

def return_rosetta_numbering(loops_as_strings):
    for string in loops_as_strings:
        pass
def convert_PDBLIST_to_sqlite3db(current_directory, filename=""):
    
    
    if filename=="":
        filename = tkFileDialog.askopenfilename(initialdir = current_directory)
    if filename ==None:
        return
    
    PDBLIST = open(filename, 'r')
    dbname = os.path.dirname(filename)+"/DATABASE.db"
    print dbname
    DB = PDB("", "", "", False, dbname)
    i = 1
    for filepath in PDBLIST:
        pdbID = os.path.basename(filepath).split(".")[0]
        DB.set_basic_options(pdbID, i, "x")
        filepath = filepath.strip()
        DB.read_pdb_into_database_flat(filepath, False, False)
        i+=1
    DB.db.close()
    print "Database written to PDBLIST directory"

"""
def extract_pdbs_from_sqlite3db(current_directory, list="", dbfile = ""):
    if list=="":
        filename = tkFileDialog.askopenfilename(initialdir = current_directory, message="PDBLIST")
        
    if filename ==None:
        return
    dir = os.path.split(PDBLIST)[0]
    
    PDBLIST = open(filename, 'r')
    
    if dbfile=="":
        dbfilename = tkFileDialog.askopenfilename(initialdir = dir, message="dbfilename")
        dir = os.path.split(PDBLIST)[0]
    if dbfilename ==None:
        return
    
    pdbdb = PDB_database()
"""

def score_PDBLIST(pdblist_path, score):
    """
    Outputs a simple pdb vs score for simple analysis.
    if pdblist_path=False, a dialog box opens.
    """
    
    if not pdblist_path:
        pdblist_path = tkFileDialog.askopenfilename(title = "PDBLIST")
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
    print "Complete"
    PDBLIST.close()
    SCORED_PDBLIST.close()
        
def convert_PDBLIST_to_rosetta_db(current_directory):
    pass


#### FASTA OUTPUT ####

def save_FASTA(pose, base_name, outfilename = False, current_directory=False, loops_as_strings = False ):
    """
    If outfilename is False, will ask for a directory using current_directory.
    If loops_as_strings is given, output FASTA of loops.
    Uses Pyrosetta...
    """
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = current_directory, message="Output FASTA to...")
    if not outfilename:
        return
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
def save_FASTA_PDBLIST(pdblist_path, outfilename=False, current_directory=False, loops_as_strings=False):
    """
    If outfilename is False, will ask for a filename
    Goes through each member of PDBLIST
    Uses pyrosetta, much slower...Could use my PDB class...
    """
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = current_directory, message="Output FASTA to...")
    if not outfilename:
        return
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

    
def exportPDBSCORE(self, cwd, score):
    """
    Exports a list of scores.
    """
    PDBLIST = tkFileDialog.askopenfile(title = "PDBLIST", initialdir = cwd)
    OUTFILE = tkFileDialog.asksaveasfile(title = "Save As...", initialdir = cwd)
    
    if PDBLIST == None or OUTFILE==None:
        return
    
    for PDBPath in PDBLIST:
        PDBPath = PDBPath.strip()
        
        p = Pose()
        pose_from_pdb(p, PDBPath)
        SCORE = score(p)
        print SCORE
        OUTFILE.write(PDBPath+":%.3f"%SCORE+"\n")
    print "Complete"
    PDBLIST.close()
    OUTFILE.close()