#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/menu.py
## @brief  Handles main menu component of the GUI 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from Tkinter import *
from window_modules import *
#from window_modules.interactive_terminal import IPython
from modules.tools import output as output_tools
from modules.tools import general_tools
from modules.tools import sequence
from modules import help as help_tools
from modules.tools import input as input_tools
from modules import calibur
import webbrowser
class Menus():
    def __init__(self, main, toolkit):
	self.main = main
	self.toolkit = toolkit
	self.MenBar=Menu(self.main)
	self.options_class = Option_System_Manager(self.toolkit.current_directory.get())
    def setTk(self):

	"""
    OPTION - Enable Antibodies

    #### Antibody Menu ####
	self.Antibody = Menu(self.MenBar, tearoff=0)
	self.Antibody.add_command(label="General Tools", command=lambda: self.shoCDRAnalysis())
	self.Antibody.add_command(label="Antibody Designer")
	self.Antibody.add_command(label="Antibody Design Analysis")

	self.MenBar.add_cascade(label = "Antibodies", menu=self.Antibody)
	"""



    #### Protein Design Menu ####
	self.MenDesign=Menu(self.MenBar, tearoff=0)
	self.MenDesign.add_command(label="Design File ToolBox", command=lambda: self.shoDesign1())
	#self.MenDesign.add_command(label="Design Ab Initio", command=lambda: self.shoProAb())
	#self.MenDesign.add_command(label="Loop Ab Initio")
	#self.MenDesign.add_command(label="Constraint File ToolBox")
	#self.MenDesign.add_command(label="Full Control ToolBox", command=lambda: self.shoFull())
	self.MenBar.add_cascade(label="Protein Design", menu=self.MenDesign)

    #### Advanced Control Menu ####
	self.FineControl=Menu(self.MenBar, tearoff=0)
	self.FineControl.add_command(label="Advanced PyMOL Visualization", command=lambda: self.toolkit.PyMOLObject.makeWindow(0, 0, Toplevel(self.main), self.toolkit.ScoreObject))
	self.FineControl.add_command(label="Full Control Toolbox", command=lambda: self.toolkit.FullControlObject.makeWindow(Toplevel(self.main)))
	self.FineControl.add_command(label="ScoreFxn Control", command =lambda: self.toolkit.ScoreObject.makeWindow(Toplevel(self.main), self.toolkit.pose))
	self.FineControl.add_separator()
	self.FineControl.add_command(label="Configure Option System",command = lambda: self.shoOptions())
	self.FineControl.add_command(label="Interactive Terminal", foreground='red',command = lambda: self.shoIPythonWindow())
	self.FineControl.add_command(label="Jump into Session", foreground='red', command = lambda: embed())
	#self.FineControl.add_command(label="Constraints")


	self.MenBar.add_cascade(label = "Advanced", menu = self.FineControl)

    #### Rosetta Tool Menu ####
	self.MenRosetta = Menu(self.MenBar, tearoff=0)
	#self.MenRosetta.add_command(label = "Return Rosetta Numbering", command = lambda: output_tools.return_rosetta_numbering())
	self.MenRosetta.add_command(label ="Setup PDB for Rosetta", command=lambda: FixPDB().runfixPDBWindow(self.main, 0, 0))
	self.MenRosetta.add_separator()
	self.MenRosetta.add_command(label ="Enable Constraints", foreground='red')
	self.MenRosetta.add_command(label ="Enable Symmetry", foreground='red')
	self.MenRosetta.add_command(label ="Enable Non-Standard Residues", foreground='red')
	self.MenRosetta.add_separator()
	##self.MenRosetta.add_command(label= "Rosetta Script Creator")
	self.MenRosetta.add_command(label= "Rosetta Command-Line Creator", command = lambda: self.shoRosettaProtocolSetup())
	##self.MenRosetta.add_command(label = "Run PsiPred")
	self.MenBar.add_cascade(label="Rosetta Tools", menu=self.MenRosetta)
    #### Molecular Dynamics Menu ####
	"""'
	self.MenMD=Menu(self.MenBar, tearoff=0)
	self.MenMD.add_command(label="Setup Fixed Atoms")
	self.MenMD.add_command(label="Setup Harmonic Constraints")
	self.MenMD.add_command(label="Setup User-Defined Forces")
	self.MenMD.add_separator()
	self.MenMD.add_command(label="Run NAMD")
	self.MenBar.add_cascade(label="Dynamics", menu=self.MenMD)
	"""''

    #### PDB Menu #### NEEDS WORK
	"""
	self.MenPDB=Menu(self.MenBar, tearoff=0)
	#self.MenPDB.add_command(label="Load PDB into Python", command = lambda: setPDBdic(tools.pdbs.pdbTools().parsePDB(filename1.get(), pdbDic)))
	self.MenPDB.add_command(label= "Setup PDB for Rosetta", command=lambda: FixPDB.runfixPDBWindow(self.main, 0, 0))
	#self.MenPDB.add_separator()
	#self.MenPDB.add_command(label = "Meiler label_'s Clean PDB", command = lambda: self.printPDB())
	#self.MenPDB.add_command(label = "Clean MD File", command = lambda: setPDBdic(tools.pdbs.pdbTools().cleanPDB(pdbDic)))
	#self.MenPDB.add_separator()
	#self.MenPDB.add_command(label= "Change Occupancy to 1.0", command= lambda: tools.pdbs.pdbTools().chaOcc(pdbDic))
	#self.MenPDB.add_command(label= "Remove Alternate Residues", command = lambda: tools.pdbs.pdbTools().RemAlt(pdbDic))
	#self.MenPDB.add_command(label= "Remove Element Column", command = lambda: tools.pdbs.pdbTools().RemEle(pdbDic))
	#self.MenPDB.add_separator()
	#self.MenPDB.add_command(label = "Save PDB File", command = lambda: setPDBdic(tools.pdbs.pdbTools().savePDB(pdbDic)))
	#self.MenBar.add_cascade(label = "PDB Tools", menu = self.MenPDB)
	self.MenRel=Menu(self.MenBar, tearoff=0)
	self.RelClose=Menu(self.MenBar, tearoff=0)

    #### Analysis Menu #### NEEDS WORK
	self.MenAnal=Menu(self.MenBar, tearoff=0)
	self.MenCluster = Menu(self.MenBar, tearoff=0)
	##self.MenCluster.add_command(label = "Cluster PDBList using Rosetta")
	##self.MenCluster.add_command(label = "Cluster PDBList using Calibur (Recommended)")
	##self.MenAnal.add_cascade(label = "Cluster", menu = self.MenCluster)
	##self.MenAnal.add_command(label="Analysis Menu", command=lambda: self.shoAnalyze())
	#self.MenAnal.add_command(label="Pose Management")
	self.MenAnal.add_separator()
	self.Menddg=Menu(self.MenBar, tearoff=0)
	##self.Menddg.add_checkbutton(label="Save To File")
	self.Menddg.add_command(label="Calculate", command=lambda: tools.analysis.interface().analyzeDDG(self.toolkit.pose))
	##self.MenAnal.add_cascade(label="Calculate ddg", menu=self.Menddg)
	self.MenBar.add_cascade(label="Analysis", menu=self.MenAnal)
	"""

    #### Import Menu #### REMOVE?
	self.MenFrag=Menu(self.MenBar, tearoff=0)
	self.MenFrag.add_command(label="Loop File", command = lambda: self.load_loop())
	##self.MenFrag.add_command(label="Design Resfile")
	#self.MenFrag.add_separator()
	#self.MenBar.add_cascade(label="Import", menu=self.MenFrag)

    #### Export Menu ####
	self.MenExport = Menu(self.MenBar, tearoff=0)
	self.MenExport.add_command(label="Save SCWRL seq File", command=lambda: self.savSeq())
	self.MenExport.add_separator()
	self.MenExport.add_command(label="Save Rosetta Loop File", foreground = 'red', command = lambda: self.savLoop())
	##self.MenExport.add_command(label="Save Rosetta ResFile",command=lambda: Design1.saveDesign())
	##self.MenExport.add_command(label="Save Rosetta Blueprint File")
	self.MenExport.add_separator()
	self.MenExport.add_command(label = "Save FASTA (Pose)", command=lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False, self.toolkit.current_directory.get() ))
	self.MenExport.add_command(label = "Save FASTA (Loops)", command = lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False, self.toolkit.current_directory.get(), self.toolkit.loops_as_strings))
	self.MenExport.add_separator()
	#self.MenExport.add_command(label="Save Database")
	#self.MenExport.add_command(label="Save loops as new PDBs")
	self.MenBar.add_cascade(label="Export", menu=self.MenExport)
	

    #### Tools Menu ####
	self.MenTools = Menu(self.MenBar, tearoff=0)

	self.MenTools.add_command(label = "Create PDBList", command = lambda: self.toolkit.PDBLIST.set(output_tools.make_PDBLIST(self.toolkit.current_directory.get())))
	self.MenTools.add_command(label = "Create PDBList Recursively", command = lambda: self.toolkit.PDBLIST.set(output_tools.make_PDBLIST_recursively(self.toolkit.current_directory.get())))
	##self.MenTools.add_command(label = "Convert PDBList to Silent File")
	self.MenTools.add_separator()
	self.MenTools.add_command(label = "Convert PDBList to SQLite3 DB", command = lambda: output_tools.convert_PDBLIST_to_sqlite3db(self.toolkit.current_directory.get(), self.toolkit.PDBLIST.get()))
	self.MenTools.add_command(label = "Convert PDBList to Rosetta DB", foreground='red')
	self.MenTools.add_separator()
	self.MenTools.add_command(label = "Cluster PDBList using Rosetta", foreground='red')
	self.MenTools.add_command(label = "Cluster PDBList using Calibur", foreground='red')
	self.MenTools.add_separator()
	self.MenTools.add_command(label = "Rename PDBs + Copy to Outpath", command = lambda: output_tools.rename_and_save(self.toolkit.PDBLIST.get()))
	self.MenTools.add_separator()

	self.MenScore = Menu(self.MenBar, tearoff=0)
	""" There are scripts in RosettaTools that do this welll...
	self.MenScore.add_command(label = "Rescore PDBList", command = lambda: output_tools.score_PDBLIST(self.toolkit.PDBLIST.get(), self.toolkit.ScoreObject))
	#self.MenScore.add_command(label = "Load Scores")
	#self.MenScore.add_command(label = "Filter  PDBList")
	#self.MenScore.add_command(label = "Filter  PDBList by Loop Energy")
	#self.MenScore.add_command(label = "Filter  PDBList by E Component")

	self.MenTools.add_cascade(label = "Score", menu=self.MenScore)
	"""

	self.MenSequence = Menu(self.MenBar, tearoff=0)
	self.MenSequence.add_command(label = "Output FASTA for Each PDB", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.PDBLIST.get(), False, self.toolkit.current_directory.get()))
	self.MenSequence.add_command(label = "Output FASTA for Each Region", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.PDBLIST.get(), False, self.toolkit.current_directory.get(), self.toolkit.loops_as_strings))
	#If you need this, you know how to program: self.MenSequence.add_command(label = "Output LOOP file for each PDB", command = lambda: output_tools.save_LOOP_PDBLIST(self.toolkit.PDBLIST.get()))
	self.MenSequence.add_command(label = "Output Design Breakdown for each Loop", foreground='red')
	self.MenTools.add_cascade(label = "Sequence", menu=self.MenSequence)
	##self.MenTools.add_command(label = "Rename Files Recursively")
	self.MenBar.add_cascade(label = "PDBList Tools", menu=self.MenTools)
	

    #### Help Menu ####
	self.MenHelp=Menu(self.MenBar, tearoff=0)
	self.MenHelp.add_command(label="About", command=lambda: help_tools.about())
	self.MenHelp.add_command(label = "License", command = lambda: help_tools.show_license())
	self.MenHelp.add_command(label="PyRosetta Tutorials", command = lambda: webbrowser.open("http://www.pyrosetta.org/tutorials"))
	self.MenHelp.add_command(label="Rosetta Manual", command = lambda: webbrowser.open("http://www.rosettacommons.org/manual_guide"))
	self.MenDevel = Menu(self.MenBar, tearoff=0)
	self.MenDevel.add_command(label = "Wiki Page", command = lambda: webbrowser.open("https://wiki.rosettacommons.org/index.php/PyRosetta_Toolkit"))
	self.MenHelp.add_cascade(label = "Developers: ", menu = self.MenDevel)
	self.MenHelDes=Menu(self.MenBar, tearoff=0)
	self.MenHelDesMut =Menu(self.MenBar, tearoff=0)
	self.MenHelDesMut.add_command(label= "Accessible Surface Area", command=lambda: help_tools.mutSA())
	self.MenHelDesMut.add_command(label = "Relative Mutability", command=lambda: help_tools.mutRM())
	self.MenHelDesMut.add_command(label = "Surface Probability", command=lambda: help_tools.mutSP())
	self.MenHelp.add_cascade(label="Mutability Data", menu=self.MenHelDesMut)

	self.MenHelp.add_separator()
	self.MenBar.add_cascade(label="Help", menu=self.MenHelp)
    
    #### Window CheckButtons ######
	self.checkbutton_fullcontrol = Checkbutton(self.main, text="FullControl")
	self.checkbutton_pymol = Checkbutton(self.main, text = "PyMOL Visualization")

    def shoTk(self):
	self.main.config(menu=self.MenBar)




##### MENU FUNCTIONS #######

    def printPDB(self):
	print pdbDic
    def savSeq(self):
	"""
	Saves SCWRL sequence file
	"""

	if not self.toolkit.input_class.loops_as_strings:
	    tkMessageBox.showerror(title="Loop", message="Please Specify what you would like to remodel...")
	    return
	else:
	    out = tkFileDialog.asksaveasfilename(initialdir = self.toolkit.input_class.current_directory.get(), title ="Save As...")
	    output_tools.saveSeqFile(self.toolkit.pose, out, self.toolkit.input_class.loops_as_strings)
	outSP = out.split("/")
	length = len(outSP)
	folder = outSP[length-2]
	if folder == "FragSets_Designs":
	    prot.LisProt1Frag.insert(END, outSP[length-1])

    def load_loop(self):
	if not self.toolkit.pose:
	    print "Please load a pose..."
	    return
	
	f = tkFileDialog.askopenfilename(initialdir = self.toolkit.current_directory.get(), title="Loop file")
	if not f:
	    return

	loops_as_strings = input_tools.load_loop_file(self.toolkit.pose, f)
	for string in loops_as_strings:
	    self.toolkit.loops_as_strings.append(string)
	    self.toolkit.input_class.loops_listbox.insert(END, string)

    def savLoop(self):
	"""
	Saves a Loop File for Rosetta
	"""
	if not self.toolkit.input_class.loops_as_strings:
	    tkMessageBox.showerror(title="Loop", message="Please Specify what you would like to remodel...")
	    return
	else:
	    out = tkFileDialog.asksaveasfilename(initialdir = self.toolkit.input_class.current_directory.get(), title="Save As...")
	    if not out:
		return
	    output_tools.savLoop(self.toolkit.pose, out, self.toolkit.loops_as_strings)
	return
    
#### WINDOWS ##### (ADD NEW WINDOWS TO THIS) #######
    def shoOptions(self):
	"""
	Main Design window interacting with options system
	"""
	#Broken?  What the fuck?
	WinOptions = Toplevel(self.main)
	self.options_class.setTk(WinOptions)
	self.options_class.shoTk()
    def shoDesign1(self):
	"""
	Main Design window for creating a ResFile
	"""
	#Broken?  What the fuck?
	WinDesign1 = Toplevel(self.main)
	design1 = Design(WinDesign1, self.toolkit.DesignDic, self.toolkit.pose)
	design1.setTk()
	design1.shoTk()
	design1.setTypes()
    def shoAnalyze(self):
	"""
	Main Analysis window
	"""

	WinAnalyze = Toplevel(self.main); Analyze1 = AnalyzePyRosetta(WinAnalyze)
	Analyze1.setTk(); Analyze1.shoTk(); Analyze1.setLis()
    def shoPdbTools(self):
	"""
	PDB Tools for cleaning pdb files, hopfully more stuff in the future..
	"""

	WinPDB = Toplevel(self.main);
	PDB1 = PdbTools(WinPDB);
	PDB1.setTk(); PDB1.shoTk()
    def shoCDRAnalysis(self):
	"""
	CDR Analysis window.  Needs to be optional as 300 mb of pdb is a lot.
	"""

	WinCDRAnal = Toplevel(self.main)
	CDRWindow = CDR_Analysis.CDRAnalysis(WinCDRAnal, 0, 0)
	seq = StringVar()
	if VarLoopSeq.get() =="":
	    seq.set("Sequence")
	else:
	    print ":"+VarLoopSeq.get()+":"
	    seq.set(VarLoopSeq.get())
	seq2 = StringVar()
	seq2.set("New")
	CDRWindow.setTk(seq, seq2)
	CDRWindow.grid()

    def shoRosettaProtocolSetup(self):
	"""
	Rosetta Protocols - Used to make commands from lists of possible options.
	"""

	RosProtSet = Toplevel(self.main)
	ProtSet = RosettaProtocols.ProtocolSetup(RosProtSet)
	ProtSet.setTk()
	ProtSet.shoTk(0, 0)
	ProtSet.setMenu(RosProtSet)
    
    def shoIPythonWindow(self):
	"""
	IPython Interactive Window.  Isolated from variables for now.
	"""
	term = IPythonView(Toplevel(self.main))
	term.pack()
	
