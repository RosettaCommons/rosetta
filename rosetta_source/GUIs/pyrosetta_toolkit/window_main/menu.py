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
import webbrowser

#from window_modules.interactive_terminal import IPython
from modules.tools import output as output_tools
from modules.tools import general_tools
from modules.tools import sequence
from modules import help as help_tools
from modules.tools import input as input_tools
from modules import calibur
from modules import ImportExport

#Windows
from window_modules.options_system.OptionSystemManager import OptionSystemManager
from window_modules.clean_pdb.FixPDBWindow import FixPDBWindow
from window_modules.rosetta_tools.RosettaProtocolBuilder import RosettaProtocolBuilder
from window_modules.design.ResfileDesignWindow import ResfileDesignWindow

class Menus():
    def __init__(self, main, toolkit):
	self.main = main
	self.toolkit = toolkit 
	self.MenBar=Menu(self.main)
	self.options_class = OptionSystemManager(self.toolkit.current_directory.get())
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


	
      #### Import ####
	self.MenImport = Menu(self.MenBar, tearoff=0)
	self.MenImport.add_command(label="Loop File", command = lambda: self.load_loop())
	
	
      #### Export ####
	self.MenExport = Menu(self.MenBar, tearoff=0)
	self.MenExport.add_command(label="SCWRL seq File", command=lambda: self.save_SCWRL_sequence_file())
	self.MenExport.add_separator()
	self.MenExport.add_command(label="Rosetta Loop File", foreground = 'red', command = lambda: self.save_loop_file())
	self.MenExport.add_command(label="Rosetta ResFile", foreground='red')
	self.MenExport.add_command(label="Rosetta Blueprint File", foreground='red')
	self.MenExport.add_separator()
	self.MenExport.add_command(label = "FASTA (Pose)", command=lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False, self.toolkit.current_directory.get() ))
	self.MenExport.add_command(label = "FASTA (Loops)", command = lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False, self.toolkit.current_directory.get(), self.toolkit.loops_as_strings))
	#self.MenExport.add_command(label="Save Database")
	#self.MenExport.add_command(label="Save loops as new PDBs")
	
    #### File Menu ####
	self.MenFile=Menu(self.MenBar, tearoff=0)
	input_object = ImportExport.InputFiles(self, self.toolkit)
	self.MenFile.add_command(label="Load PDB", command=lambda: input_object.choose_load_pose())
	self.MenFile.add_command(label="Load PDB list", command=lambda: input_object.set_PDBLIST())
	self.MenFile.add_cascade(label="Import", menu=self.MenImport)
	self.MenFile.add_cascade(label="Export", menu=self.MenExport)
	self.MenFile.add_checkbutton(label="Set Pymol Observer", variable=self.toolkit.PyMOLObject.auto_send) #this option should be set only once.
	self.MenFile.add_command(label="Configure Option System",command = lambda: self.show_OptionsSystemManager())
	self.MenFile.add_command(label ="Setup PDB for Rosetta", command=lambda: FixPDBWindow().runfixPDBWindow(self.main, 0, 0))
	self.MenFile.add_command(label ="Enable Constraints", foreground='red')
	self.MenFile.add_command(label ="Enable Symmetry", foreground='red')
	self.MenFile.add_command(label ="Enable Non-Standard Residues", foreground='red')
	self.MenBar.add_cascade(label="File", menu=self.MenFile)
	self.MenFile.add_separator()
	self.MenFile.add_command(label= "Rosetta Command-Line Creator", command = lambda: self.show_RosettaProtocolBuilder())
	
    #### Protein Design Menu ####
	self.MenDesign=Menu(self.MenBar, tearoff=0)
	self.MenDesign.add_command(label="Setup Resfile for PDB", command=lambda: self.show_ResfileDesignWindow())
	self.MenDesign.add_command(label="Setup Blueprint for PDB", foreground='red')
	self.MenDesign.add_command(label="Structure Editing and Grafting", foreground='red')
	self.MenBar.add_cascade(label="Protein Design", menu=self.MenDesign)

    #### Advanced Control Menu ####
	self.FineControl=Menu(self.MenBar, tearoff=0)
	self.FineControl.add_command(label="PyMOL Visualization", command=lambda: self.toolkit.PyMOLObject.makeWindow(0, 0, Toplevel(self.main), self.toolkit.ScoreObject))
	self.FineControl.add_command(label="Full Control Toolbox", command=lambda: self.toolkit.FullControlObject.makeWindow(Toplevel(self.main)))
	self.FineControl.add_command(label="ScoreFxn Control + Creation", command =lambda: self.toolkit.ScoreObject.makeWindow(Toplevel(self.main), self.toolkit.pose))
	self.FineControl.add_separator()
	self.FineControl.add_command(label="Interactive Terminal", foreground='red',command = lambda: self.show_IpythonWindow())
	self.FineControl.add_command(label="Jump into Session", foreground='red', command = lambda: embed())
	self.MenBar.add_cascade(label = "Advanced", menu = self.FineControl)

	
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
	self.MenHelp.add_command(label="Rosetta Glossary", command=lambda: help_tools.print_glossary())
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

    def shoTk(self):
	self.main.config(menu=self.MenBar)




##### MENU FUNCTIONS #######

    def save_SCWRL_sequence_file(self):
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

    def save_loop_file(self):
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
    
    
    
#### WINDOWS ##### (ADD NEW WINDOWS TO THIS THAT NEED TO BE SET UP) #######

    def show_OptionsSystemManager(self):
	"""
	Main Design window interacting with options system
	"""
	
	top_level_tk = Toplevel(self.main)
	self.options_class.setTk(top_level_tk)
	self.options_class.shoTk()
    def show_ResfileDesignWindow(self):
	"""
	Main Design window for creating a ResFile
	"""
	
	top_level_tk = Toplevel(self.main)
	resfile_design_window = ResfileDesignWindow(top_level_tk, self.toolkit.DesignDic, self.toolkit.pose)
	resfile_design_window.setTk()
	resfile_design_window.shoTk()
	resfile_design_window.setTypes()

    def show_RosettaProtocolBuilder(self):
	"""
	Rosetta Protocols - Used to make commands from lists of possible options.
	"""

	top_level_tk = Toplevel(self.main)
	rosetta_protocol_builder = RosettaProtocolBuilder(top_level_tk)
	rosetta_protocol_builder.setTk()
	rosetta_protocol_builder.shoTk(0, 0)
	rosetta_protocol_builder.setMenu(top_level_tk)
    
    def show_IpythonWindow(self):
	"""
	IPython Interactive Window.  Isolated from variables for now.
	"""
	term = IPythonView(Toplevel(self.main))
	term.pack()
	
