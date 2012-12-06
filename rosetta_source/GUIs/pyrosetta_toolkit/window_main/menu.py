#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/menu.py
## @brief  Handles main menu component of the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Tk and Python Imports
from Tkinter import *
import webbrowser
import tkSimpleDialog

#Module Imports
from modules.tools import output as output_tools
from modules.tools import analysis as analysis_tools
from modules.tools import general_tools
from modules.tools import sequence

from modules import help as help_tools
from modules.tools import input as input_tools
from modules import calibur

from modules.protocols.DockingProtocols import DockingProtocols
from modules.protocols.DesignProtocols import DesignProtocols
from modules.protocols.LowResLoopModelingProtocols import LowResLoopModelingProtocols
from modules.protocols.HighResLoopModelingProtocols import HighResLoopModelingProtocols

#Window Imports
from window_modules.options_system.OptionSystemManager import OptionSystemManager
from window_modules.clean_pdb.FixPDBWindow import FixPDBWindow
from window_modules.rosetta_tools.RosettaProtocolBuilder import RosettaProtocolBuilder
from window_modules.design.ResfileDesignWindow import ResfileDesignWindow
#from window_modules.interactive_terminal import IPython
from window_modules.ligand_ncaa_ptm_manager.ligand_ncaa_ptm_manager import ligand_ncaa_ptm_manager
import global_variables

class Menus():
    def __init__(self, main, toolkit):
	self.main = main
	self.toolkit = toolkit
	self.main_menu=Menu(self.main)
	


    def setTk(self):
	"""
	Sets up the Menu.  Order is important here.
	"""


	self._set_file_menu()
	self._set_advanced_menu()
	self._set_protocols_menu()
	self._set_protein_design_menu()
	self._set_pdblist_menu()
	self._set_help_menu()

    def shoTk(self):
	"""
	Shows the menu via main.
	"""

	self.main.config(menu=self.main_menu)
	
    def _set_file_menu(self):
	"""
	Sets import, export, and main File menu
	"""

      #### Import ####
	self.import_menu = Menu(self.main_menu, tearoff=0)
	self.import_menu.add_command(label="Param PathList File", command = lambda: self.toolkit.input_class.load_param_list())
	self.import_menu.add_command(label="Rosetta Loop File", command = lambda: self.toolkit.input_class.load_loop())
	#self.import_menu.add_command(label="Rosetta Resfile", command = lambda: input_tools.load_resfile())


      #### Export ####
	self.export_menu = Menu(self.main_menu, tearoff=0)
	self.export_menu.add_command(label="SCWRL seq File", command=lambda: self.save_SCWRL_sequence_file())
	self.export_menu.add_separator()
	self.export_menu.add_command(label="Rosetta Loop File", command = lambda: output_tools.save_loop_file(self.toolkit.pose, self.toolkit.input_class.loops_as_strings))
	self.export_menu.add_command(label="Rosetta Basic ResFile", command = lambda: output_tools.save_basic_resfile(self.toolkit.pose))
	self.export_menu.add_command(label="Rosetta Basic Blueprint File", command = lambda: output_tools.save_basic_blueprint(self.toolkit.pose))
	self.export_menu.add_separator()
	self.export_menu.add_command(label="Param PathList File", command = lambda: output_tools.save_param_path_list(self.toolkit.input_class.param_paths))
	self.export_menu.add_separator()
	self.export_menu.add_command(label = "FASTA (Pose)", command=lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False ))
	self.export_menu.add_command(label = "FASTA (Loops)", command = lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.outname.get(), False, self.toolkit.input_class.loops_as_strings))
	#self.export_menu.add_command(label="Save Database")
	#self.export_menu.add_command(label="Save loops as new PDBs")

    #### File Menu ####
	self.file_menu=Menu(self.main_menu, tearoff=0)
	
	self.file_menu.add_command(label="Load PDB", command=lambda: self.toolkit.input_class.choose_load_pose())
	self.file_menu.add_command(label="Load PDBList", command=lambda: self.toolkit.input_class.set_PDBLIST())
	self.file_menu.add_cascade(label="Import", menu=self.import_menu)
	self.file_menu.add_cascade(label="Export", menu=self.export_menu)
	self.file_menu.add_checkbutton(label="Set Pymol Observer", variable=self.toolkit.pymol_class.auto_send) #this option should be set only once.
	self.file_menu.add_command(label = "Show Pose in PyMOL", command = lambda: self.toolkit.pymol_class.pymover.apply(self.toolkit.pose))
	#self.file_menu.add_checkbutton(label="Set stdout to terminal", variable = self.toolkit.terminal_output)
	self.file_menu.add_command(label="Configure Option System",command = lambda: self.show_OptionsSystemManager())
	self.file_menu.add_command(label ="Setup PDB for Rosetta", command=lambda: FixPDBWindow().runfixPDBWindow(self.main, 0, 0))
	self.main_menu.add_cascade(label="File", menu=self.file_menu)
	self.file_menu.add_separator()
	self.file_menu.add_command(label= "Rosetta Command-Line Creator", command = lambda: self.show_RosettaProtocolBuilder())

    def _set_advanced_menu(self):
	"""
	Sets the advanced control menu
	"""

	self.advanced_menu=Menu(self.main_menu, tearoff=0)
	self.analysis_menu = Menu(self.main_menu, tearoff=0)
	self.analysis_menu.add_command(label = "Interface Analyzer", command=lambda: analysis_tools.analyze_interface(self.toolkit.pose, self.toolkit.score_class.score))
	self.analysis_menu.add_command(label = "Packing Analyzer", command = lambda: analysis_tools.analyze_packing(self.toolkit.pose))
	self.analysis_menu.add_command(label = "Loops Analyzer", command = lambda: analysis_tools.analyze_loops(self.toolkit.pose, self.toolkit.input_class.loops_as_strings))
	self.analysis_menu.add_command(label = "VIP Analyzer", foreground='red', command = lambda: analysis_tools.analyze_vip(self.toolkit.pose, self.toolkit.score_class.score, self.toolkit.pwd))
	self.advanced_menu.add_cascade(label = "Analysis", menu = self.analysis_menu)
	self.advanced_menu.add_separator()
	self.advanced_menu.add_command(label ="Enable Constraints", command = lambda: input_tools.add_constraints_to_pose_and_scorefunction(self.toolkit.pose, self.toolkit.score_class.score))
	#self.advanced_menu.add_command(label ="Enable Symmetry", foreground='red')
	
	#NonStandard
	self.non_standard_menu = Menu(self.main_menu, tearoff=0)
	self.non_standard_menu.add_command(label = "Ligand/NCAA/PTM Manager", command = lambda: self.show_ligand_ncaa_ptm_manager())
	self.non_standard_menu.add_separator()
	self.non_standard_menu.add_command(label = "Download Rosetta NCAA Rotamer Library", command=lambda: webbrowser.open("http://carl.bio.nyu.edu/~renfrew/ncaa/"))
	self.non_standard_menu.add_command(label = "Convert molfile to param files", command = lambda: output_tools.output_molfile_to_params())
	self.non_standard_menu.add_separator()
	self.non_standard_menu.add_command(label = "Documentation", command = lambda: webbrowser.open("http://www.pyrosetta.org/obtaining-and-preparing-ligand-pdb-files"))
	
	self.advanced_menu.add_cascade(label ="Enable NCAA/PTM/Ligands", menu=self.non_standard_menu)
	self.advanced_menu.add_separator()
	self.advanced_menu.add_command(label="PyMOL Visualization", command=lambda: self.toolkit.pymol_class.makeWindow(0, 0, Toplevel(self.main), self.toolkit.score_class))
	self.advanced_menu.add_command(label="ScoreFunction Control and Creation", command =lambda: self.toolkit.score_class.makeWindow(Toplevel(self.main), self.toolkit.pose))
	self.advanced_menu.add_command(label="Per Residue Control and Analysis", command=lambda: self.toolkit.fullcontrol_class.makeWindow(Toplevel(self.main)))
	self.advanced_menu.add_separator()
	#self.advanced_menu.add_command(label="Extract PDB from SQLite3 DB", command = lambda: output_tools.extract_pdb_from_sqlite3db())
	#self.advanced_menu.add_command(label="Interactive Terminal", foreground='red',command = lambda: self.show_IpythonWindow())
	#self.advanced_menu.add_command(label="Jump into Session", foreground='red', command = lambda: embed())
	self.main_menu.add_cascade(label = "Advanced", menu = self.advanced_menu)


	
    def _set_protein_design_menu(self):
	"""
	Sets the Protein Design menu.
	"""

	self.protein_design_menu=Menu(self.main_menu, tearoff=0)
	self.protein_design_menu.add_command(label="Setup Resfile for PDB", command=lambda: self.show_ResfileDesignWindow())
	self.protein_design_menu.add_command(label="Setup Blueprint for PDB", foreground='red')
	self.protein_design_menu.add_command(label="Structure Editing and Grafting", foreground='red')
	self.main_menu.add_cascade(label="Protein Design", menu=self.protein_design_menu)

    def _set_protocols_menu(self):
	"""
	Menu for running protocols through PyRosetta.
	"""
	
	self.protocols_menu = Menu(self.main_menu, tearoff=0)

	self.protocols_menu.add_command(label = "Set processors to use", command = lambda: self.toolkit.output_class.processors.set(tkSimpleDialog.askinteger(title="Processesors", prompt="Please set the number of processess you wish to create for protocol runs.", initialvalue=self.toolkit.output_class.processors.get())))
	self.protocols_menu.add_command(label = "Enable MPI Mode", foreground = 'red')

	self.protocols_menu.add_separator()
	
	#Setup Protocol Classes. Should this go to Main File? They are very light weight classes.
	self.design_class = DesignProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.docking_class = DockingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.low_res_loop_modeling_class = LowResLoopModelingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.high_res_loop_modeling_class = HighResLoopModelingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	
	#Design
	
	self.design_protocols = Menu(self.main_menu, tearoff=0)
	self.design_protocols.add_command(label = "FixedBB", command = lambda: self.design_class.packDesign())
	self.design_protocols.add_command(label = "Grafting", foreground='red')
	self.design_protocols.add_command(label = "Remodel", foreground='red')
	
	#Docking
	
	self.docking_protocols = Menu(self.main_menu, tearoff=0)
	self.docking_protocols.add_command(label = "Low Resolution", command = lambda: self.docking_class.low_res_dock())
	self.docking_protocols.add_command(label = "High Resolution", command = lambda: self.docking_class.high_res_dock())
	self.docking_protocols.add_command(label = "Ligand", foreground='red')
	
	#Loop Modeling

	self.loop_modeling_protocols=Menu(self.main_menu, tearoff=0)
	self.loop_modeling_protocols_low = Menu(self.main_menu,tearoff=0)
	self.loop_modeling_protocols_low.add_command(label = "CCD", command = lambda: self.low_res_loop_modeling_class.default_CCD())
	self.loop_modeling_protocols_low.add_command(label = "KIC", command = lambda: self.low_res_loop_modeling_class.default_KIC())
	self.loop_modeling_protocols_high = Menu(self.main_menu, tearoff=0)
	self.loop_modeling_protocols_high.add_command(label = "CCD", command = lambda: self.high_res_loop_modeling_class.default_CCD())
	self.loop_modeling_protocols_high.add_command(label = "KIC", command = lambda: self.high_res_loop_modeling_class.default_KIC())
	self.loop_modeling_protocols.add_cascade(label = "Low Resolution", menu = self.loop_modeling_protocols_low)
	self.loop_modeling_protocols.add_cascade(label = "High Resolution", menu = self.loop_modeling_protocols_high)
	
	#General Modeling (Homology/Abinitio)
	self.general_protocols = Menu(self.main_menu, tearoff=0)
	self.servers = Menu(self.main_menu, tearoff = 0); #Because why not?
	
	#Servers.  Because there are tons of work in these, and I can't add everything.
	self.servers.add_command(label = "ROSIE Servers", command = lambda: webbrowser.open("http://rosie.rosettacommons.org/"))
	self.servers.add_separator()
	self.servers.add_command(label = "Backrub", command = lambda: webbrowser.open("https://kortemmelab.ucsf.edu/backrub/cgi-bin/rosettaweb.py?query=submit"))
	self.servers.add_command(label = "Fragments", command = lambda: webbrowser.open("http://robetta.bakerlab.org/fragmentsubmit.jsp"))
	self.servers.add_command(label = "Interface Alanine Scan", command = lambda: webbrowser.open("http://robetta.bakerlab.org/alascansubmit.jsp"))
	self.servers.add_command(label = "DNA Interface Scan", command = lambda: webbrowser.open("http://robetta.bakerlab.org/dnainterfacescansubmit.jsp"))
	
	self.general_protocols.add_cascade(label = "Web Servers", menu = self.servers)
	
	self.general_protocols.add_separator()
	self.general_protocols.add_command(label = "Ab Initio", foreground='red')
	self.general_protocols.add_command(label = "Homology Modeling", foreground='red')
	self.protocols_menu.add_cascade(label = "General", menu = self.general_protocols)
	self.protocols_menu.add_separator()
	self.protocols_menu.add_cascade(label = "Loop Modeling", menu = self.loop_modeling_protocols)
	self.protocols_menu.add_cascade(label = "Docking", menu = self.docking_protocols)
	self.protocols_menu.add_cascade(label = "Design", menu = self.design_protocols)
	
	
	self.main_menu.add_cascade(label = "Protocols", menu=self.protocols_menu)

    def _set_pdblist_menu(self):
	"""
	Sets a Menu for interacting with a simple PDBList.  Useful since Rosetta is pretty much all about thousands of models.
	"""

	self.pdblist_tools_menu = Menu(self.main_menu, tearoff=0)
	self.pdblist_analysis_menu = Menu(self.main_menu, tearoff=0)
	#There are scripts in RosettaTools that do this welll...
	self.pdblist_analysis_menu.add_command(label = "Energy vs RMSD", foreground='red')
	self.pdblist_analysis_menu.add_command(label = "Energy vs RMSD (Region(s))", foreground='red')
	self.pdblist_analysis_menu.add_separator()
	self.pdblist_analysis_menu.add_command(label = "Rescore PDBList", command = lambda: output_tools.score_PDBLIST(self.toolkit.input_class.PDBLIST.get(), self.toolkit.score_class.score))
	self.pdblist_analysis_menu.add_command(label = "Load Scores", foreground='red')
	self.pdblist_analysis_menu.add_command(label = "Get top model", foreground='red')
	self.pdblist_analysis_menu.add_command(label = "Get top % models", foreground='red')
	#self.pdblist_analysis_menu.add_command(label = "Filter  PDBList by Loop Energy")
	#self.pdblist_analysis_menu.add_command(label = "Filter  PDBList by E Component")
	self.pdblist_tools_menu.add_cascade(label = "Analysis", menu=self.pdblist_analysis_menu)

	self.sequence_menu = Menu(self.main_menu, tearoff=0)
	self.sequence_menu.add_command(label = "Output FASTA for Each PDB", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.input_class.PDBLIST.get(), False))
	self.sequence_menu.add_command(label = "Output FASTA for Each Region", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.input_class.PDBLIST.get(), False, self.toolkit.input_class.loops_as_strings))
	#If you need this, you know how to program: self.sequence_menu.add_command(label = "Output LOOP file for each PDB", command = lambda: output_tools.save_LOOP_PDBLIST(self.toolkit.input_class.PDBLIST.get()))
	self.sequence_menu.add_command(label = "Output Design Breakdown for each Loop", foreground='red')
	self.pdblist_tools_menu.add_cascade(label = "Sequence", menu=self.sequence_menu)
	self.pdblist_tools_menu.add_separator()
	self.pdblist_tools_menu.add_command(label = "Create PDBList", command = lambda: self.toolkit.input_class.PDBLIST.set(output_tools.make_PDBLIST()))
	self.pdblist_tools_menu.add_command(label = "Create PDBList Recursively", command = lambda: self.toolkit.input_class.PDBLIST.set(output_tools.make_PDBLIST_recursively()))
	self.pdblist_tools_menu.add_separator()
	self.pdblist_tools_menu.add_command(label = "Cluster PDBList using Calibur", foreground='red')
	#self.pdblist_tools_menu.add_command(label = "Convert PDBList to SQLite3 DB", command = lambda: output_tools.convert_PDBLIST_to_sqlite3db(self.toolkit.input_class.PDBLIST.get()))
	#self.pdblist_tools_menu.add_command(label = "Extract PDBList from SQLite3 DB", command = lambda: output_tools.extract_pdbs_from_sqlite3db(self.toolkit.input_class.PDBLIST.get()))
	#self.pdblist_tools_menu.add_separator()
	self.pdblist_tools_menu.add_command(label = "Rename All PDBs Recursively + Copy to Outpath", command = lambda: output_tools.rename_and_save(self.toolkit.input_class.PDBLIST.get()))

	self.main_menu.add_cascade(label = "PDBList Tools", menu=self.pdblist_tools_menu)

    def _set_help_menu(self):
	"""
	Sets the help Menu.  TK kinda sucks for formating dialog windows.  Just FYI.
	"""

	self.help_menu=Menu(self.main_menu, tearoff=0)
	self.help_menu.add_command(label="About", command=lambda: help_tools.about())
	self.help_menu.add_command(label = "License", command = lambda: help_tools.show_license())
	self.help_menu.add_separator()
	self.help_menu.add_command(label="Rosetta Manual", command = lambda: webbrowser.open("http://www.rosettacommons.org/manual_guide"))
	self.help_menu.add_command(label="Rosetta Glossary", command=lambda: help_tools.print_glossary())
	self.help_menu.add_command(label="Rosetta Forums", command = lambda: webbrowser.open("www.rosettacommons.org/forum"))
	self.help_menu.add_command(label="Rosetta BugTracker", command = lambda: webbrowser.open("www.bugs.rosettacommons.org"))
	self.help_menu.add_separator()
	self.help_menu.add_command(label="PyRosetta Tutorials", command = lambda: webbrowser.open("http://www.pyrosetta.org/tutorials"))
	self.help_devel_menu = Menu(self.main_menu, tearoff=0)
	self.help_devel_menu.add_command(label = "Wiki Page", command = lambda: webbrowser.open("https://wiki.rosettacommons.org/index.php/PyRosetta_Toolkit"))
	self.help_menu.add_cascade(label = "Developers: ", menu = self.help_devel_menu)
	self.help_desmut_menu =Menu(self.main_menu, tearoff=0)
	self.help_desmut_menu.add_command(label= "Accessible Surface Area", command=lambda: help_tools.mutSA())
	self.help_desmut_menu.add_command(label = "Relative Mutability", command=lambda: help_tools.mutRM())
	self.help_desmut_menu.add_command(label = "Surface Probability", command=lambda: help_tools.mutSP())
	self.help_menu.add_cascade(label="Mutability Data", menu=self.help_desmut_menu)
	
	self.main_menu.add_cascade(label="Help", menu=self.help_menu)






##### MENU FUNCTIONS #######





#### WINDOWS ##### (ADD NEW WINDOWS TO THIS THAT NEED TO BE SET UP) #######

    def show_ligand_ncaa_ptm_manager(self):
	top_level_tk = Toplevel(self.main)
	ptm = ligand_ncaa_ptm_manager(self.toolkit.input_class, self.toolkit.score_class, self.toolkit.pose)
	ptm.setTk(top_level_tk)
	ptm.shoTk(0, 0)
	
    def show_OptionsSystemManager(self):
	"""
	Main Design window interacting with options system
	"""

	top_level_tk = Toplevel(self.main)
	self.toolkit.options_class.setTk(top_level_tk)
	self.toolkit.options_class.shoTk()
	
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
