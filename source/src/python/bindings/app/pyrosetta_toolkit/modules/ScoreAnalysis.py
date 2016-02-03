#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/ScoreAnalysis
## @brief  Class for analyzing a scored PDB list/FASC/SC.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os
import re
from multiprocessing import Process

#Tkinter Imports
from Tkinter import *
import tkMessageBox
import tkSimpleDialog
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools.analysis import rmsd
from app.pyrosetta_toolkit.window_main import global_variables
from collections import defaultdict

class ScoreAnalysis:
    """
    This class reads fasc, sc, or a ScoredPDBList.  It only cares about total score- for now.
    Will get top scores, write them to a txt file or copy them to a new directory.
    Can do rmsd vs energy - but for now, we only use 1 processor to do this
    """
    def __init__(self):
        
        self.score_map = defaultdict(); #[float score]:[path/name].
        self.score_pairs = []; #[[float score, string path], [float score, string path]] - For Redundancy
        
        self.score_rmsd_triplet = []; #[[float score, float RMSD, string path],[x, x, x]]
        self.score_rmsd_triplet_loops = defaultdict(); #[loop_string]:score_rmsd_triplet
        self.top_score_map = defaultdict()
        self.top_scoring_by_percent_map = defaultdict()
        self.top_scoring_by_number_map = defaultdict()
        self.filepath=""
        
        #self.read_scores()
        
    def set_filepath(self, filepath):
        """
        Sets filepaths and loads data
        """
        if not filepath:return
        self.filepath = filepath
        self.read_scores()
        
    def read_scores(self, filepath = None):
        """
        Reads scores - Either .fasc, .fc, or .txt (.txt only for now)
        """
        if filepath:
            self.filepath = filepath
            
        if re.search(".txt", self.filepath):
            self.read_PDBList(self.filepath)
        elif re.search(".fasc", self.filepath):
            self.read_FASC(self.filepath)
        elif re.search(".sc", self.filepath):
            self.read_SC(self.filepath)
        else:
            print "Could not read filetype"
            return
        print "\nScores read. \n"
        
    def get_score_map(self):
        return self.score_map
    
    def get_score_pairs(self):
        return self.score_pairs
    
    def get_top_score_map(self):
        return self.top_score_map
    
    def get_top_scoring(self, ask_copy_results= True):
        """
        Gets top score - sets self.top_score_map.  Runs copy_results
        """
        if not self.score_map:print "Please load scores.";return
        scores = sorted(self.score_map)
        top_score = scores[0]
        
        print "Top Score: %.2f"%top_score+" "+self.score_map[top_score]
        fullpath = self.score_map[top_score]
        self.top_score_map[fullpath]=top_score
        if ask_copy_results:
            self.copy_results(self.top_score_map, "TOP_SCORE")
        else:
            return fullpath
        #return top_score, fullpath
    
    def get_top_scoring_by_percent(self, percent=False, ask_copy_results=True):
        """
        Gets top percent of poses - sets self.top_scoring_by_percent.  Runs copy_results
        """
        if not self.score_pairs:print "Please load scores.";return
        if not percent:

            percent = tkSimpleDialog.askfloat(title="Cutoff", prompt="Please enter a percent for top cutoff", initialvalue=.05)
            if not percent:return
        if percent>1.0:
            percent = percent/float(100)
        total_num = round(percent*len(self.score_pairs))
        print repr(percent)
        for i in range(0, int(total_num)):
            score = self.score_pairs[i][0]; fullpath = self.score_pairs[i][1]
            print "%.2f"%score+" "+fullpath
            self.top_scoring_by_percent_map[fullpath]=score
        
        if ask_copy_results:   
            self.copy_results(self.top_scoring_by_percent_map, "TOP_%.2f"%percent+"_PERCENT")
        else:
            return self.top_scoring_by_percent_map
        
    def get_top_scoring_by_number(self, top_number=False, ask_copy_results = True):
        """
        Gets top scoring poses - Sets self.top_scoring_by_number_map.  Runs copy_results
        """
        if not self.score_pairs:print "Please load scores.";return
        if not top_number:
            top_number = tkSimpleDialog.askinteger(title="Cutoff", prompt="Please enter how many top scoring poses you would like", initialvalue=10)
            if not top_number:return
        
        print repr(top_number)
        for i in range(0, top_number):
            score = self.score_pairs[i][0]; fullpath = self.score_pairs[i][1]
            print "%.2f"%score+" "+fullpath
            self.top_scoring_by_number_map[fullpath]=score
        
        if ask_copy_results:
            self.copy_results(self.top_scoring_by_number_map, "TOP_"+str(top_number))
        else:
            return self.top_scoring_by_number_map

    def get_score_vs_rmsd(self, native_pose, loops_as_strings):
        if not self.score_pairs: print "Please load scores.";return
        for pair in self.score_pairs:
            
            fullpath = pair[1]
            score = pair[0]
            p = Pose()
            pose_from_file(p, fullpath)
            rms, loop_rms_map = rmsd(native_pose, p, loops_as_strings)
            triplet = [score, rms, fullpath]
            self.score_rmsd_triplet.append(triplet)
            if loops_as_strings:
                for loop_string in loop_rms_map:
                    if not self.score_rmsd_triplet_loops.has_key(loop_string):
                        self.score_rmsd_triplet_loops[loop_string]=[]
                        self.score_rmsd_triplet_loops[loop_string].append(triplet)
                    else:
                        self.score_rmsd_triplet_loops[loop_string].append(triplet)
        
        print "Score vs RMSD data now loaded into GUI.."
        #Output Results
        outpath = tkFileDialog.asksaveasfilename(title="Save as..", initialdir=global_variables.current_directory)
        if not outpath:
            return
            #outpath = os.path.dirname(self.filepath)
        global_variables.current_directory = os.path.dirname(outpath)
        
        
        #result = tkMessageBox.askyesno(title="Save as Db?", message="Save to SQLite3db in addition to a txt file?")
        self.write_score_vs_rmsd_to_txt(outpath)
        #if result:
            #self.write_score_vs_rmsd_to_db(outpath)
        
        
    def plot_score_vs_rmsd(self, native_pose, loops_as_strings):
        if not self.score_rmsd_triplet:
            self.get_score_vs_rmsd(native_pose, loops_as_strings)
        pass
    
    
    def write_score_vs_rmsd_to_db(self, outpath):
        if not self.score_rmsd_triplet:
            return
        
        try:
            import sqlite3
        except ImportError:
            print "SQL functions unavailable"
            print "Saving txt file instead."
            self.write_score_vs_rmsd_to_txt(outpath)
            return
            
            
        pass
    
    def write_score_vs_rmsd_to_txt(self, outpath):
        if not self.score_rmsd_triplet:return
        OUTFILE = open(outpath, 'w')
        header = "region score rmsd fullpath\n"
        OUTFILE.write(header)
        
        for triplet in sorted(self.score_rmsd_triplet, reverse=True):
            line = "full %.3f"%triplet[0]+" %.3f"%triplet[1]+' '+triplet[2]+"\n"
            OUTFILE.write(line)
        if self.score_rmsd_triplet_loops:
            for loop_string in sorted(self.score_rmsd_triplet_loops):
                line=loop_string+" "
                for triplet in sorted(self.score_rmsd_triplet_loops[loop_string], reverse=True):
                    line = line+"%.3f"%triplet[0]+" %.3f"%triplet[1]+' '+triplet[2]+"\n"
                    OUTFILE.write(line)
        print "Results written to : "+outpath
        OUTFILE.close()
        
    
    def read_PDBList(self, filepath):
        """
        Reads scored PDBList (path score)
        """
        
        FILE = open(filepath, 'r')
        for line in FILE:
            line = line.strip()
            lineSP = line.split()
            score = float(lineSP[1]); fullpath = lineSP[0]
            pair = [score, fullpath]
            self.score_map[score]=fullpath
            self.score_pairs.append(pair)
        self.score_pairs = sorted(self.score_pairs)    
        FILE.close()
        self.score_pairs = sorted(self.score_pairs) 
        
    def read_FASC(self, filepath):
        FILE = open(filepath, 'r')
        for line in FILE:
            lineSP = line.strip().split()
            if lineSP[0]=="filename:":
                fullpath = lineSP[1]
                total_score = float(lineSP[3])
                self.score_map[total_score]=fullpath
                pair = [total_score, fullpath]
                self.score_pairs.append(pair)
        self.score_pairs = sorted(self.score_pairs)
        FILE.close()
        
    def read_SC(self, filepath):
        FILE = open(filepath, 'r')
        print "Assuming PDB's are within the same path as the score file."
        for line in FILE:
            lineSP = line.strip().split()
            if lineSP[1]=="score":
                continue
            score = float(lineSP[1])
            filename = lineSP[-1]; #Does not give full path. - So, we assume .sc is within the same directory.
            fullpath = os.path.dirname(filepath)+"/"+filename+"/.pdb"
            pair = [score, fullpath]
            self.score_pairs.append(pair)
            self.score_map[score]=fullpath
        self.score_pairs = sorted(self.score_pairs)
        FILE.close()
        
    def copy_results(self, path_dict, type, ask = True, dir_suffix = ""):
        """
        Should run after top x is found.  Asks user to copy results.
        path_dict is [path]:[float score]
        type is a string for type of result - top score, top scoring, top percent, top_by_number.  Used to name files and directories
        """
        if ask:
            result = tkMessageBox.askyesno(title="Copy", message="Copy resultant structures?")
            if not result:return
        
        new_dir = os.path.dirname(self.filepath)+dir_suffix+"/TOP"
        if not os.path.exists(new_dir): os.mkdir(new_dir)
        
        new_dir2 = new_dir+"/"+type.upper()
        if not os.path.exists(new_dir2): os.mkdir(new_dir2)
        print "Results are in "+new_dir +" and /"+type.upper()
        
        file_name = new_dir+'/'+os.path.basename(self.filepath).split(".")[0]+'_'+type.lower()+'.txt'
        FILE = open(file_name, 'w')
        
        for path in sorted(path_dict, key = path_dict.get):
            print path;
            FILE.write(path+"\t%.3f"%path_dict[path]+"\n")
            os.system("cp "+path+' '+new_dir2+'/'+os.path.basename(path))
        FILE.close()
        
    def copy_results2(self, path_dict, txt_out_dir, pdb_out_dir, top_pdblist_out_name):
        """
        Same as copy results, but we specify where we put top models and the top models text.
        """
        if not os.path.exists(txt_out_dir): os.mkdir(txt_out_dir)

        if not os.path.exists(pdb_out_dir): os.mkdir(pdb_out_dir)
        print "Results are in "+txt_out_dir +" and "+pdb_out_dir

        file_name = txt_out_dir+'/'+top_pdblist_out_name+'.txt'
        FILE = open(file_name, 'w')
        RAW_LIST = open(pdb_out_dir+"/ORDERED_PDBLIST.txt")
        for path in sorted(path_dict, key = path_dict.get):
            #print path;
            FILE.write(path+"\t%.3f"%path_dict[path]+"\n")
            RAW_LIST.write(path+"\n")
            os.system("cp "+path+' '+pdb_out_dir+'/'+os.path.basename(path))
        FILE.close()
        RAW_LIST.close()
        
    
if __name__ == '__main__':
    """
    For testing.
    """
    
    ScoredPDBList = ""
    Tk()
    analyzer = ScoreAnalysis(ScoredPDBList)
    analyzer.get_top_scoring()
    analyzer.get_top_scoring_by_percent(5.0)
    analyzer.get_top_scoring_by_number(3)
    
    
