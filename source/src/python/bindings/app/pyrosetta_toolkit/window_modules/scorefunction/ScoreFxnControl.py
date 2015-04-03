#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/scorefunction/ScoreFxnControl.py
## @brief  Scorefunction window.  Holds current scorefunction accross GUI.  Use to set scorefunction permanently
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.scoring.methods import *

#Python Imports
import glob
from os import environ as env
from os import path
import re

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.window_main import global_variables
#from app.pyrosetta_toolkit.modules.RegionalScoring import *


class ScoreFxnControl():
    """
    This is the class where you can pick from any of rosetta's energy functions and classes, as well as edit them on the fly, save them, and load
    custom built ones.  Useful for testing new scorefunctions.
    Use this class to set the scorefunction globally.
    Needs a bit of refactoring.
    """
    
    def __init__(self):
        self.weight_dir = rosetta.rosetta_database_from_env() +"/scoring/weights/"
        self.wd = self.location()[0]
        self.readDefaults()
        self.score = create_score_function_ws_patch(self.ScoreType.get(), self.ScorePatch.get()); #This should be a score function.  We are going to have this class hold all score information.  
        
        #Fix to have Hbonds in scoring and total score.  Only 3 years too late.  Bargh.
        emopts = EnergyMethodOptions(self.score.energy_method_options())
        emopts.hbond_options().decompose_bb_hb_into_pair_energies( True)
        self.score.set_energy_method_options(emopts)
        
        self.pose = Pose(); #Should this be here?
        
        self.ScoreOptionsdef = StringVar(); self.ScoreOptionsdef.set("Score Pose")
        self.OPTIONS = [
            "Set ScoreFunction",
            "Score Pose",
        ]
      
        self.common_scorefunctions = [
            "talaris2013",
            "score12prime",
            "score12",
            "orbitals",
            "soft_rep",
            "soft_rep_design",
            "mm_std",
            "ligand",
            "interchain_cen",
            "docking",
            "cen_std",
        ]
        
        self.show_all_scorefunctions = IntVar(); self.show_all_scorefunctions.set(False)
        
    def makeWindow(self, main, p, r=0, c=0):
        """
        Creates the window to edit, analyze, and explore score functions and score terms.
        """
        
        #Listboxes with Labels
        
        self.pose = p
        self.main = main
        self.main.title("ScoreFunction Control")
        
        #Menu:
        self.MenBar = Menu(self.main)
        self.export = Menu(self.MenBar, tearoff=0)
        self.export.add_command(label = "Save to WD", command = lambda: self.saveToWD())
        self.export.add_command(label = "Save as..", command = lambda: self.saveAS())
        self.export.add_separator()
        #self.export.add_command(label = "Score list of PDB's - OUTPUT SCORE", command = lambda: self.exportSCORE())
        self.export.add_command(label = "Score list of PDB's - OUTPUT PDB_PATH:SCORE", command = lambda: output_tools.score_PDBLIST(False, self.score))
        self.MenBar.add_cascade(label = "Export", menu=self.export)
        
        self.options = Menu(self.MenBar, tearoff=0)
        #self.options.add_command(label = "Set Working Directory (WD)", command = lambda: self.setWD())
        self.options.add_checkbutton(label = "Show ALL ScoreFunctions", variable=self.show_all_scorefunctions)
        self.options.add_command(label = "Set Default Score Function (Does not save edits)", command = lambda: self.saveDefaults())
        self.MenBar.add_cascade(label = "Options", menu = self.options)
        
        self.main.config(menu=self.MenBar)
        
        self.TypeLabel = Label(self.main, text = "Score Function")
        self.ScoreTypeListbox = Listbox(self.main, width=30, height = 18)
        self.ScoreTypeScroll = Scrollbar(self.main)
        self.ScoreTypeListbox.config(yscrollcommand=self.ScoreTypeScroll.set); self.ScoreTypeScroll.config(command=self.ScoreTypeListbox.yview)
        
        self.PatchLabel = Label(self.main, text = "Patches")
        self.ScorePatchListbox = Listbox(self.main, width=30, height = 18)
        self.ScorePatchScroll = Scrollbar(self.main)
        self.ScorePatchListbox.config(yscrollcommand=self.ScorePatchScroll.set); self.ScorePatchScroll.config(command = self.ScorePatchListbox.yview)
        
        self.NonzeroLabel = Label(self.main, text = "NonZero Terms (Double click to change)")
        self.NonzeroListbox = Listbox(self.main, width=30, height = 18)
        self.NonzeroScroll = Scrollbar(self.main)
        self.NonzeroListbox.config(yscrollcommand=self.NonzeroScroll.set); self.NonzeroScroll.config(command = self.NonzeroListbox.yview)
        
        self.ZeroLabel = Label(self.main, text = "Zero Terms (Double click to add)")
        self.ZeroListbox = Listbox(self.main, width=30, height = 18)
        self.ZeroScroll = Scrollbar(self.main)
        self.ZeroListbox.config(yscrollcommand=self.ZeroScroll.set); self.ZeroScroll.config(command = self.ZeroListbox.yview)
        
        #Showing Score Selected:
        self.ShoTypeLabel = Label(self.main, textvariable = self.ScoreType)
        self.ShoPatchLabel = Label(self.main, textvariable = self.ScorePatch)
        
        
        #First Actions after selected:
        self.change_scorefunction_button = Button(self.main, text = "Change ScoreFunction", command = lambda: self.updateScore())
        self.change_patch_button = Button(self.main, text = "Change/Add Patch", command = lambda: self.updatePatch())
        self.ScoreButton = Button(self.main, text = "Score Pose", command = lambda: self.scoreOption("Score Pose"))
        self.ScoreOptions = OptionMenu(self.main, self.ScoreOptionsdef, *self.OPTIONS)
        self.RemovePatchButton = Button(self.main, text = "Remove Patch", command = lambda: self.removePatch(), justify=CENTER)
        self.ScoreTypeListbox.bind("<Double-Button-1>", lambda event: self.updateScore())
        self.ScorePatchListbox.bind("<Double-Button-1>", lambda event: self.updatePatch())
        self.NonzeroListbox.bind("<Double-Button-1>", lambda event: self.changeWeight(self.NonzeroListbox))
        self.ZeroListbox.bind("<Double-Button-1>", lambda event: self.changeWeight(self.ZeroListbox))
        
        
        #ShoTk
        self.TypeLabel.grid(row = r+0, column=c+0); self.PatchLabel.grid(row=r+0, column=c+2)
        self.ScoreTypeListbox.grid(row=r+1, column=c+0, rowspan=6); self.ScoreTypeScroll.grid(row=r+1, column=c+1, rowspan=6, sticky=E+N+S)
        self.ScorePatchListbox.grid(row=r+1, column=c+2, rowspan=6); self.ScorePatchScroll.grid(row=r+1, column=c+3, rowspan=6, sticky=E+N+S)
        
        self.ShoTypeLabel.grid(row=r+7, column=c+0); self.ShoPatchLabel.grid(row=r+7, column=c+2)
        self.change_scorefunction_button.grid(row=r+8, column=c+0, columnspan=2, sticky=W+E); self.change_patch_button.grid(row=r+8, column=c+2, sticky=W+E, columnspan=2)
        #elf.ScoreButton.grid(row=r+9, column=c+0, columnspan = 4, sticky=W+E);
        #self.ScoreOptions.grid(row=r+8, column=c+1, sticky = W+E)
        self.RemovePatchButton.grid(row=r+10, column=c+2, columnspan = 2, sticky=W+E)
        
        self.NonzeroLabel.grid(row=r+11, column=c+0);self.ZeroLabel.grid(row=r+11, column=c+2)
        self.NonzeroListbox.grid(row = r+12, column = c+0); self.NonzeroScroll.grid(row=r+12, column=c+1, sticky=E+N+S)
        self.ZeroListbox.grid(row = r+12, column = c+2); self.ZeroScroll.grid(row=r+12, column=c+3, sticky=E+N+S)
        
        #This Gets all of the score types from PyRosetta Database
        self.populateScoreList()
        self.updateScoreandTerms()
        self.show_all_scorefunctions.trace_variable('w', self.update_scorelist_callback)
    
    
    def update_scorelist_callback(self, name, index, mode):
        self.populateScoreList()
        
    ###Public Methods###
    def set_scorefunction(self, name, patch=None):
        self.ScoreType.set(name)
        self.ScorePatch.set(patch)
        self.scoreOption("Set ScoreFunction")
        
        #Check if window has been run and is running somewhere
        try:
            self.updateScoreandTerms()
        except AttributeError:
            return
        
    def readDefaults(self):
        """
        Opens text file with default settings.
        """
        
        self.ScoreType = StringVar();
        self.ScorePatch = StringVar();
        DEFAULT = open(self.wd+"/defaults.txt", 'r')
        self.defaultDic = dict()
        for line in DEFAULT:
            lineSP = line.strip().split("=")
            self.defaultDic[lineSP[0]]=lineSP[1]
            
        self.setDefaults()
        DEFAULT.close()
        print self.ScoreType.get()
        print self.ScorePatch.get()
        
    def setDefaults(self):
        self.ScoreType.set(self.defaultDic["SCORE"])
        if self.defaultDic.has_key("PATCH"):
            self.ScorePatch.set(self.defaultDic["PATCH"])
        if self.defaultDic.has_key("WD"):
            self.wd.set(self.defaultDic["WD"])
            
    def saveDefaults(self):
        """
        Saves text with default settings.
        """
        DEFAULT = open(self.wd+"/defaults.txt", 'w')
        if self.ScorePatch.get() ==None:
            DEFAULT.write("SCORE="+self.ScoreType.get())
        else:
            DEFAULT.write("SCORE="+self.ScoreType.get()+"\n")
            DEFAULT.write("PATCH="+self.ScorePatch.get())
        DEFAULT.close()
        
    def populateScoreList(self):
        """
        Populates the ScoreFunction and Patch listboxes with files or common types.
        """
        
        self.ScoreTypeListbox.delete(0, END)
        print repr(rosetta.rosetta_database_from_env())
        if self.show_all_scorefunctions.get():
            
            
            Scorefiles = glob.glob(rosetta.rosetta_database_from_env() +"/scoring/weights/*.wts")
            
            for scorefile in sorted(Scorefiles):
                scorefile = os.path.basename(scorefile).split(".")[0]
                self.ScoreTypeListbox.insert(END, scorefile)
        else:
            for scorefile in self.common_scorefunctions:
                self.ScoreTypeListbox.insert(END, scorefile)
                
        self.ScorePatchListbox.insert(END, "None")
        
        
        Patchfiles = glob.glob(rosetta.rosetta_database_from_env()+"/scoring/weights/*.wts_patch")
        for patchfile in sorted(Patchfiles):
            patchfile = os.path.basename(patchfile).split(".")[0]
            self.ScorePatchListbox.insert(END, patchfile)
    
    def populateETerms(self, nonzeroList, zeroList):
        self.NonzeroListbox.delete(0, END)
        self.ZeroListbox.delete(0, END)
        for term in sorted(nonzeroList):
            self.NonzeroListbox.insert(END, term)
        
        for term in sorted(zeroList):
            self.ZeroListbox.insert(END, term)
            
    def scoreOption(self, Option):
        if Option=="Set ScoreFunction":
            if self.ScorePatch.get() == "None":
                self.score = create_score_function(self.ScoreType.get())
                print "Score Function Set to: "+self.ScoreType.get()
            else:

                self.score = create_score_function_ws_patch(self.ScoreType.get(), self.ScorePatch.get())
                
                #Fix to have Hbonds in scoring and total score.
                emopts = EnergyMethodOptions(self.score.energy_method_options())
                emopts.hbond_options().decompose_bb_hb_into_pair_energies( True)
                self.score.set_energy_method_options(emopts)

            return
        elif Option=="Set Temp Score":
            if self.ScorePatch.get() == "None":
                scoretemp = create_score_function(self.ScoreType.get())
            else:
                scoretemp = create_score_function_ws_patch(self.ScoreType.get(), self.ScorePatch.get())
            return scoretemp
        elif Option=="Score Pose":
            scoretemp = self.scoreOption("Set Temp Score")
            
            if self.pose.total_residue()==0:
                print "No Pose loaded..."
                return
            
            print scoretemp(self.pose)
            self.score.show(self.pose)


            return
        elif Option=="Breakdown ScoreFxn":
            #self.scoreOption("Set ScoreFunction")
            y = str(self.score.weights())
            ylis = list(y)
            ylis2 = ylis[2:(len(ylis)-2)]
            y = "".join(ylis2)
            ysp = y.split(") (")
            zeroscores = []
            nonzeroscores = []
            for ele in ysp:
                eleSP = ele.split("; ")
                if eleSP[1] == "0":
                    zeroscores.append(eleSP[0].strip())
                else:
                    nonzeroscores.append(ele.strip())
            
            #print "\nZero Weights:"
            #for ele in sorted(zeroscores):
                #print ele
            #print "\nNonZero Weights"
            #for ele in sorted(nonzeroscores):
                #print ele
            return zeroscores, nonzeroscores
        
        elif Option=="Graph ScoreFxn":
            (zeroscores, nonzeroscores) = self.scoreOption("Breakdown ScoreFxn")
            try:
                from Modules import plotting
            except ImportError:
                return
            
            #Setup of the Data:
            x = []; y=[]
            data = (x, y)
            Settings = dict()
            Settings['width']=.5
            Settings['size']='x-small'
            Settings['color']='b'
            Settings['rot']='vertical'
            
            for ele in nonzeroscores:
                eleSP = ele.split(";")
                xPiece = eleSP[0].strip(); yPiece = float(eleSP[1].strip())
                data[0].append(xPiece); data[1].append(yPiece)
            plotting.generalPlot().BarPlot(data, "Function", "Weight", "Score Function Breakdown: "+self.ScoreType.get()+" "+self.ScorePatch.get(), Settings)
            return
        else:
            return
        
    def updateScore(self):
        self.ScoreType.set(self.ScoreTypeListbox.get(self.ScoreTypeListbox.curselection()))
        self.scoreOption("Set ScoreFunction")
        self.updateScoreandTerms()
    def updatePatch(self):
        self.ScorePatch.set(self.ScorePatchListbox.get(self.ScorePatchListbox.curselection()))
        self.scoreOption("Set ScoreFunction")
        self.updateScoreandTerms()
    def removePatch(self):
        self.ScorePatch.set("None")
        self.scoreOption("Set ScoreFunction")
        self.updateScoreandTerms()
    def changeWeight(self, box):
        """
        Changes the weights on the term selected.  Updates the Listbox.  Updates the energy term.  If now NonZero, adds it to the NonZeroList and updates.  ViceVersa.
        """
        term = box.get(box.curselection())
        termSP = term.split(";")
        if len(termSP)>=2:
            newWeight = tkSimpleDialog.askfloat(title="New", prompt = "Enter New Weight", initialvalue=float(termSP[1]))
            
            if  newWeight == 0:
                print "Zeroing E term"
                self.score.set_weight(eval(termSP[0]), 0)
                self.zero, self.nonzero = self.scoreOption("Breakdown ScoreFxn")
                self.populateETerms(self.nonzero, self.zero)
            elif  not newWeight: return
            else:
                print "Changing E term"
                self.score.set_weight(eval(termSP[0]), newWeight)
                self.zero, self.nonzero = self.scoreOption("Breakdown ScoreFxn")
                self.populateETerms(self.nonzero, self.zero)
        else:
            newWeight = tkSimpleDialog.askfloat(title="weight", prompt = "Enter New Weight")
            if newWeight == None or newWeight == 0:
                return
            else:
                self.score.set_weight(eval(termSP[0]), newWeight)
                self.zero, self.nonzero = self.scoreOption("Breakdown ScoreFxn")
                self.populateETerms(self.nonzero, self.zero)
          
    def ret_Refs(self):
        self.getRefs()
        return self.refLine
    
    def saveAS(self, outfilepath=False):
        """
        Saves New Energy Weights
        """
        
        if not outfilepath:
            outfilepath= tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory)
        if not outfilepath:
            return
        self.zero, self.nonzero = self.scoreOption("Breakdown ScoreFxn")
        
        global_variables.current_directory = path.dirname(outfilepath)
        FILE = open(outfilepath, 'w')
        FILE.write(self.ret_Refs())
        for term in sorted(self.nonzero):
            termSP = term.split("; ")
            FILE.write(termSP[0]+" "+termSP[1]+"\n")
        FILE.close()
        
    def saveToWD(self):
        name = tkSimpleDialog.askstring(title = "Save", prompt = "Name")
        if not name:
            return
        FILE = open(global_variables.current_directory + "/"+name, 'w'); FILE.write(self.ret_Refs())
        for term in sorted(self.nonzero):
            termSP = term.split("; ")
            FILE.write(termSP[0]+" "+termSP[1]+"\n")
        FILE.close()
        
    def setWD(self):
        self.wd = tkFileDialog.askdirectory(title = "Set WD", initialdir = self.wd)
        
    def updateScoreandTerms(self):
        """
        Updates Term listboxes.
        """
        #self.scoreOption("Set ScoreFunction")
        self.zero, self.nonzero = self.scoreOption("Breakdown ScoreFxn")
        self.populateETerms(self.nonzero, self.zero)
    
    def getRefs(self):
        """
        Loads the name of the ScoreFunction and parses the reference energies as a 1 line string.
        We then write this when outputing the ScoreFunction.
        """
        
        SCOREFILE = open(rosetta.rosetta_database_from_env() +"/scoring/weights/"+self.ScoreType.get()+".wts")
        for line in SCOREFILE:
            if re.search("METHOD_WEIGHTS", line):
                self.refLine = line
                break
        SCOREFILE.close()
        
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported/ran from anywhere.
        """
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP
