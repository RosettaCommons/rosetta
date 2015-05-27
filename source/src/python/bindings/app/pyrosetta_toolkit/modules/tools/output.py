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
import tarfile
import multiprocessing
from multiprocessing import Process
import time

#Tkinter Imports
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox

#Toolkit Imports
import sequence
import loops
from app.pyrosetta_toolkit.modules.SQLPDB import *
from app.pyrosetta_toolkit.modules.tools import loops as loop_tools
from app.pyrosetta_toolkit.window_main import global_variables
from app.pyrosetta_toolkit.modules.definitions import restype_definitions
from collections import defaultdict

def dumpPDB(p, native_pose, filepath, score, overwrite=False):
    """
    Dumps the pose using the Py Job Distributor
    """
    if not overwrite:
        jd=PyJobDistributor(filepath, 100000, score); #This number is high so that it outputs a pose even if one with the name already exists...
        #native_pose = Pose()
        #pose_from_pdb(native_pose, infile)
        jd.native_pose=native_pose
        jd.output_decoy(p)
        os.remove(jd.current_name+".in_progress")
    else:
        print "Removing .fasc + .pdb with same output name."
        
        if os.path.exists(filepath+".fasc"):
            os.remove(filepath+".fasc")
        if os.path.exists(filepath+"_1.pdb"):
            os.remove(filepath+"_1.pdb")
        native_pose.dump_pdb(filepath+"_1.pdb")
        output_scorefile(p, filepath, filepath+"_1.pdb", filepath+".fasc", score, 1, native_pose)
        
    print "Pose written to "+os.path.dirname(filepath)
    
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
    
def save_loop_file(p, regions, ask_info=True, outfilename=False):
    """
    Saves a Rosetta Loop file.  Also asks to discard residues or not.
    Migrate to use Regions.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    if not regions:
        print "\n No loops to save...\n"
        return

    
    if ask_info:
        
        ask_cut_points = tkMessageBox.askyesno(message="Define Cutpoints?", default=tkMessageBox.NO)

        discard_loops =  tkMessageBox.askyesno(message = "Discard Phi/Psi and build using ideal bond lengths and angles?")
    else:
        ask_cut_points = False
        discard_loops = True
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output loop file to...")
    if not outfilename: return
    global_variables.current_directory=os.path.dirname(outfilename)
    
    FILE = open(outfilename, 'w')
    for region in regions:
        
        start = region.get_rosetta_start(p)
        end = region.get_rosetta_end(p)
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

def save_basic_resfile(p, outname=None):
    """
    Saves an empty resfile, numbered by PDB with NATRO designation
    If outname is not given, will ask where to save.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    ResDic = dict()
    if not outname:
        outname = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output resfile to...")
        if not outname: return
    global_variables.current_directory=os.path.dirname(outname)
    
    save_resfile_w_designdic(p, ResDic, outname)
    
def save_resfile_w_designdic(p, ResDic, filename):
    """
    Saves a resifile, readable by PyRosetta and Rosetta.
    ResDic can be empty.
    ResDic is [string pdbNum:pdbChain]:[array name:three_letter:one_letter string]
    If NC - ResDic should be [string pdbNum:pdbChain]:[array 'NC':residue string]
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
        #res - string pdbNum:chain
        if res in ResDic:
            x = ""
            print res
            for residue_string in ResDic[res]:
                print residue_string
                if residue_string == "NATRO":
                    line = pdbStr + chainStr + "  NATRO"+"\n"                      
                elif residue_string == "NATAA":
                    line = pdbStr + chainStr + "  NATAA"+ "\n"
                elif residue_string == "ALLAA":
                    line = pdbStr + chainStr + "  ALLAA" + "\n"
                elif residue_string.split(":")[0]=="NC":
                    type = residue_string.split(":")[1]+" "
                    x = x+" NC "+type
                    #Looks like for each NC you need NC designation.
                    line = pdbStr + chainStr + x + "\n";
                else:
                    residuesAll = residue_string.split(":")
                    x = x + residuesAll[2]
                    line = pdbStr + chainStr + "  PIKAA  " + x + "\n"; #Just gets rewritten if more residues to add
            FILE.write(line)
        else:
            #If residue not found in Resdic, give NATRO
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
    

def saveSeqFile(p, fileout=None, loops_as_strings=None):
    """
    Saves a scwrl seq file.  Needs to be more robust to account for all Seths changes ( some NCAA, carbohydrates) when it is eventually released.
    Needs to migrate to use region class.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    if not fileout:
        fileout = tkFileDialog.asksaveasfilename(initialdir=global_variables.current_directory)
    if not fileout:return
    global_variables.current_directory = os.path.dirname(fileout)
    newList = loop_tools.loopArea(p, loops_as_strings)
    seq = createSeqFile(p, newList)
    FILE = open(fileout, 'w')
    FILE.write(seq+"\n")
    FILE.close()
    print "Seq file saved.."
    return


def save_basic_blueprint(p, output=True, outfilename=None):
    """
    Saves a basic blueprint file to be manually edited.
    If output is false, returns a string of the file for manipulation.
    If outfilename is not given, will ask where to save.
    """
    if not p.total_residue():
        print "\n No pose loaded...\n"
        return
    
    out_string = ""
    define = restype_definitions.definitions()
    seq = p.sequence()
    for i in range(1, p.total_residue()+1):
        pdb_num = p.pdb_info().number(i)
        single_letter_code = seq[i-1]
        out_string = out_string+repr(pdb_num)+" "+single_letter_code+" . NATRO\n"
    
    if output:
        if not outfilename:
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
    Makes a list of PDB's from a directory.  Does not walk directory. 
    Later realize could have used find command...
    """
    if directory=="":
        directory = tkFileDialog.askdirectory(title = "Choose directory with PDB files", initialdir = global_variables.current_directory)
        if not directory: return
        global_variables.current_directory=directory
    
    contains = tkSimpleDialog.askstring(title="Contains...", prompt="Separate mutliple match criteria by a coma...", initialvalue=".pdb")
    containsSP = contains.split(",")
    FILES = os.listdir(directory)
    
    filenum=1
    if len(FILES)<=1:
        print "No PDBs found.  Returning."
        return
    
    matches = []
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
            print "File "+repr(filenum)+": "+name
            p = os.path.join(directory, name)
            filenum+=1
            matches.append(p)
    
    if matches:
        NEWFILE = open(directory+"/PDBLIST.txt", 'w')
        for match in matches:
            NEWFILE.write(match+"\n")
        NEWFILE.close()
        print "File saved as 'PDBLIST.txt' in directory specified."
        return directory+"/PDBLIST.txt"
    else:
        print "No matches found.."
        return None

def make_PDBLIST_recursively(directory=""):
    if directory=="":
        directory = tkFileDialog.askdirectory(initialdir = global_variables.current_directory)
        if not directory: return
        global_variables.current_directory=directory
    
    contains = tkSimpleDialog.askstring(title="Contains...", prompt="Separate mutliple match criteria by a coma...", initialvalue=".pdb,")
    NEWFILE = open(directory+"/PDBLIST_RECURSIVE.txt", 'w')
    filenum = 1
    containsSP = contains.split(",")
    matches = []
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
                matches.append(p)
                filenum+=1
    
    if matches:
        NEWFILE = open(directory+"/PDBLIST.txt", 'w')
        for match in matches:
            NEWFILE.write(match+"\n")
        NEWFILE.close()
        print "File saved as 'PDBLIST.txt' in directory specified."
        return directory+"/PDBLIST.txt"
    else:
        print "No matches found.."
        return None

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
    DB = SQLPDB("", "", "", False, dbname)
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

def rescore_single_pdb(path, scorefunction, manager_dict):
    """
    Used by score_PDBList for multiprocessing rescoring of decoy structures.
    Only speed up if PDBs are large.
    """
    p = Pose()
    try:
        pose_from_pdb(p, path)

    except PyRosettaException:
        print "Cannot Load "+path+ " Try using -ignore_unrecognized_residues in options window..."
        return
    
    #Extra measure of protection to make sure we get to the end of the function.
    if p.total_residue()>0:
        print "Loaded PDB"
        score = scorefunction(p)
        manager_dict[path]=score
    else:
        return
    
def score_PDBLIST(pdblist_path, score, processors, output_class = None):
    """
    Outputs a simple pdb vs score for simple analysis.
    if pdblist_path=False, a dialog box opens.
    will grab the number of processors from output_class (processor StringVar variable) and attempt multiprocessing rescoring of all PDBs.
    """
    
    if not pdblist_path:
        pdblist_path = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title = "PDBLIST")
        if not pdblist_path: return
        global_variables.current_directory = os.path.dirname(pdblist_path)
        PDBLIST = open(pdblist_path, 'r')
    else:
        PDBLIST = open(pdblist_path, 'r')
    SCORED_PDBLIST = open(os.path.dirname(pdblist_path)+"/SCORED_PDBLIST.txt", 'w')

    if processors==1:
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
        
    #Multiprocessing - Slower for extremely small PDBs, great for large ones.
    else:
        if output_class:
            output_class.terminal_output.set(1)
        manager = multiprocessing.Manager()
        result_map = manager.dict(); #[path]:[score]
        workers = []
        i=1
        for line in PDBLIST:
            line = line.strip()
            if not os.path.exists(line):
                print "Could not find "+line
                print "Skipping"
                continue
            result_map[line]="NA"
            worker = Process(name = line, target=rescore_single_pdb, args=(line, score, result_map))
            workers.append(worker)
            i+=1
        total_allowed_jobs = processors
        print "Total allowed jobs: "+repr(total_allowed_jobs)
        total_running_jobs = 0
        job_complete=False
              
        #Run the protocol
        while not job_complete:

            
            time.sleep(1)
            for worker in workers:
                if worker.is_alive():
                    pass
                    #There is no code for worker hasn't started yet that I can figure out.  So, something else needs to check it!
                elif result_map[worker.name]!="NA":
                    if worker.exitcode!=0:
                        print "%s.exitcode = %s" %(worker.name, worker.exitcode)
                    
                    workers.pop(workers.index(worker)); #If the job is done, pop it.
                    total_running_jobs-=1
                    print "Total running jobs: "+repr(total_running_jobs)
                    print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                    
            if len(workers)==0:
                job_complete=True
                break
            
            if total_running_jobs<total_allowed_jobs:
                for worker in workers:
                    if not worker.is_alive():
                        print "Starting Worker"
                        try:
                            worker.start()
                        except AssertionError:
                            continue
                        print "Total running jobs: "+repr(total_running_jobs)
                        print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                        total_running_jobs+=1
                        if total_running_jobs>=total_allowed_jobs: break
        
            if total_running_jobs==0:
                job_complete=True

        d = defaultdict()

        #Convert Multiprocessing dictionary to regular one so it can be sorted
        for path in result_map.keys():
            d[path] = result_map[path]

        for path in sorted(d, key=d.get):
            e = d[path]
            print path+"\t%.3f\n"%e
            SCORED_PDBLIST.write(path+"\t%.3f\n"%e)
            
        print "\nComplete. File written to "+os.path.dirname(pdblist_path)+"/SCORED_PDBLIST.txt\n"
        if output_class:
            output_class.terminal_output.set(0)
        SCORED_PDBLIST.close()

    return os.path.dirname(pdblist_path)+"/SCORED_PDBLIST.txt"

#### FASTA OUTPUT ####

def save_FASTA(pose, base_name, outfilename = None, regions = None ):
    """
    If outfilename is False, will ask for a directory using current_directory.
    If loops_as_strings is given, output FASTA of loops. Base_name is used as label >base_name (region) for fasta
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
    if regions:
        region_array = regions.get_regions()
        for region in region_array:
            if not region.region_exists(pose):continue
            else:
                seq = region.get_sequence(pose)
            header = ">"+base_name+" "+region.get_region_string_with_all_residues(pose)
            OUTFILE.write(header+"\n")
            OUTFILE.write(seq+"\n")
    else:
        seq = pose.sequence()
        OUTFILE.write(">"+base_name+"\n")
        OUTFILE.write(seq+"\n")
    OUTFILE.close()
    print "FASTA written."
    return

def save_FASTA_PDBLIST(pdblist_path, outfilename=None, regions=None):
    """
    If outfilename is False, will ask for a filename
    Goes through each member of PDBLIST
    Uses pyrosetta to get sequence.
    """
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output FASTA to...")
        if not outfilename:return
        global_variables.current_directory = os.path.dirname(outfilename)
    
    OUTFILE = open(outfilename, 'w')
    PDBLIST = open(pdblist_path, 'r')
    i = 1
    for pdbpath in PDBLIST:
        pdbpath = pdbpath.strip()
        if not pdbpath: continue

        pdb = os.path.basename(pdbpath)
        pdbID = pdb.split(".")[0]
        print pdbID
        pose = Pose()
        region_array = regions.get_regions()
        try:
            pose_from_pdb(pose, pdbpath)
        except PyRosettaException:
            print "Could not load.. "+pdbID+"..continueing.."
            continue
        if regions:
            for region in region_array:
                #if not region.region_exists(pose):continue
                seq = region.get_sequence(pose)
                

                header = ">"+repr(i)+"_"+os.path.basename(pdbpath)+" "+region.get_region_string_with_all_residues(pose)
                print header
                print seq
                OUTFILE.write(header+"\n")
                OUTFILE.write(seq+"\n\n")
        else:
            seq=pose.sequence()
            OUTFILE.write(">"+pdbID+" "+pdbpath+"\n")
            OUTFILE.write(seq+"\n")
        i+1

    PDBLIST.close()
    OUTFILE.close()
    print "Fasta written..."
    return

    
def exportPDBSCORE(score):
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
    
def output_molfile_to_params():
    """
    Uses molfile to params in pyrosetta bindings to convert.
    Maybe should be converted to a window for more options.
    """
    
    print "Using molfile_to_params.py script located in pyrosetta/toolbox/molfile2params written by Ian W Davis.  For more options, please use script."
    script_path = os.environ["PYROSETTA"]+"/toolbox/molfile2params/molfile_to_params.py"
    
    if not os.path.exists(script_path):
        print "Untarring script"
        extract_path = os.environ["PYROSETTA"]+"/toolbox"
        tar_path = extract_path+"/molfile2params.tar.gz"
        
        if tarfile.is_tarfile(tar_path):
            tfile = tarfile.open(tar_path)
            tfile.extractall(extract_path)
        else:
            print "Could not extract tar file."
            return
    
    mdl_file = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title = "Open MDL, MOL, MOL2, or SDF file")
    if not mdl_file:return
    global_variables.current_directory=os.path.dirname(mdl_file)
    
    output_kinemage = tkMessageBox.askyesno(title = "kinemage", message="Output kinemage file for ligand visualization?")
    
    options = " "+mdl_file+" "
    if output_kinemage:
        options = options+"-k "
    
    outdir = os.path.dirname(mdl_file)+"/"+os.path.basename(mdl_file).split(".")[0]
    if not os.path.exists(outdir): os.mkdir(outdir)
    
    prefix = outdir+"/"+os.path.basename(mdl_file).split(".")[0]
    
    options = options +"-c --clobber "+"-p "+prefix
    
    print "Running molfile_to_params with these options: "+options
    os.system("python "+script_path+options)
    print "Parameters generated. Output directed to: "+outdir
    
    return
    
def save_param_path_list(array_of_paths, outfilename=False):
    """
    Saves a file of paths.
    """
    if not array_of_paths:
        print "No extra params enabled."
        return
    if not outfilename:
        outfilename = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title="Output Parm pathList to...")
    if not outfilename:return
    global_variables.current_directory = os.path.dirname(outfilename)
    FILE = open(outfilename, 'w')
    d = dict()
    #Uniqify the output
    for path in array_of_paths:d[path]=0
    
    for path in d:
        FILE.write(path+"\n")
    FILE.close()
