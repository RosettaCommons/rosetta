#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/menu.py
## @brief  Handles main menu component of the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import webbrowser
import time

#Tkinter Imports
from Tkinter import *
import tkSimpleDialog
import tkFileDialog
import tkMessageBox

#Module Imports
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.modules.tools import analysis as analysis_tools
from app.pyrosetta_toolkit.modules.tools import general_tools
from app.pyrosetta_toolkit.modules.tools import sequence

from app.pyrosetta_toolkit.modules import help as help_tools
from app.pyrosetta_toolkit.modules.tools import input as input_tools
from app.pyrosetta_toolkit.modules import calibur
from app.pyrosetta_toolkit.modules.ScoreAnalysis import ScoreAnalysis
from app.pyrosetta_toolkit.modules.DesignBreakdown import DesignBreakdown
from app.pyrosetta_toolkit.modules.protocols.DockingProtocols import DockingProtocols
from app.pyrosetta_toolkit.modules.protocols.DesignProtocols import DesignProtocols
from app.pyrosetta_toolkit.modules.protocols.LowResLoopModelingProtocols import LowResLoopModelingProtocols
from app.pyrosetta_toolkit.modules.protocols.HighResLoopModelingProtocols import HighResLoopModelingProtocols
from app.pyrosetta_toolkit.modules.protocols.GraftingProtocols import GraftMoverWindow
from app.pyrosetta_toolkit.modules.protocols.AnalysisProtocols import AnalysisProtocols
from app.pyrosetta_toolkit.modules.protocols.FloppyTailProtocol import FloppyTailProtocol
from app.pyrosetta_toolkit.modules.protocols.MinimizationProtocols import MinimizationProtocols

from app.pyrosetta_toolkit.window_modules.clean_pdb.FixPDBWindow import FixPDBWindow


from app.pyrosetta_toolkit.window_modules.design.ResfileDesignWindow import ResfileDesignWindow
from app.pyrosetta_toolkit.window_modules.ligand_ncaa_ptm_manager.ligand_ncaa_ptm_manager import ligand_ncaa_ptm_manager
from app.pyrosetta_toolkit.window_modules.full_control.FullControlWindow import FullControlWindow
from app.pyrosetta_toolkit.window_modules.insert_bfactor.InsertBFactor import InsertBFactor

import global_variables
from app.pyrosetta_toolkit.window_main.IO import SessionIO

#### Import of other GUIs in rosetta_source.  Nessessary check for eventual inclusion into PyRosetta Binaries.
flag_file_builder_import = True
try:
    from rosetta_flag_file_builder.RosettaFlagFileBuilder import RosettaFlagFileBuilder
except ImportError:
    flag_file_builder_import = False
####


class Menus():
    def __init__(self, main, toolkit):
	"""
	@main: Tk
	@toolkit: pyrosetta_toolkit
	"""
	self.main = main
	self.toolkit = toolkit; #Need to figure out how to tell Python what this 'toolkit' is so IDEs can make sense out of it!
	self.main_menu=Menu(self.main)
	self.score_analyzer = ScoreAnalysis(); #Needs to exist as instance due to holding data from load data.  Eventually these functions can go into it's own window.


    def setTk(self):
	"""
	Sets up the Menu.  Order is important here.
	"""


	self._set_file_menu()
	self._set_options_menu()
	self._set_visualization_menu()
	self._set_advanced_menu()
	self._set_protocols_menu()
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
        self.import_menu.add_command(label="GUI Session", command = lambda:SessionIO.SessionIO(self.toolkit).load_session())
        self.import_menu.add_separator()
	self.import_menu.add_command(label="Param PathList File", command = lambda: self.toolkit.input_class.load_param_list())
	self.import_menu.add_command(label="Rosetta Loop File", command = lambda: self.toolkit.input_frame.load_loop())
	#self.import_menu.add_command(label="Rosetta Resfile", command = lambda: input_tools.load_resfile())


      #### Export ####
	self.export_menu = Menu(self.main_menu, tearoff=0)
        self.export_menu.add_command(label="GUI Session", command = lambda: SessionIO.SessionIO(self.toolkit).save_session())
        self.export_menu.add_separator()
	self.export_menu.add_command(label="SCWRL seq File", command=lambda: output_tools.saveSeqFile(self.toolkit.pose, None, self.toolkit.input_class.loops_as_strings))
	self.export_menu.add_separator()
	self.export_menu.add_command(label="Rosetta Loop File", command = lambda: output_tools.save_loop_file(self.toolkit.pose, self.toolkit.input_class.regions))
	self.export_menu.add_command(label="Rosetta Basic ResFile", command = lambda: output_tools.save_basic_resfile(self.toolkit.pose))
	self.export_menu.add_command(label="Rosetta Basic Blueprint File", command = lambda: output_tools.save_basic_blueprint(self.toolkit.pose))
	self.export_menu.add_separator()
	self.export_menu.add_command(label="Param PathList File", command = lambda: output_tools.save_param_path_list(self.toolkit.input_class.param_paths))
	self.export_menu.add_separator()
	self.export_menu.add_command(label = "FASTA (Pose)", command=lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.output_class.outname.get(), False ))
	self.export_menu.add_command(label = "FASTA (Regions)", command = lambda: output_tools.save_FASTA(self.toolkit.pose, self.toolkit.output_class.outname.get(), False, self.toolkit.input_class.regions))
	#self.export_menu.add_command(label="Save Database")
	#self.export_menu.add_command(label="Save loops as new PDBs")

    #### File Menu ####
        self.file_menu=Menu(self.main_menu, tearoff=0)

        self.file_menu.add_command(label="Load PDB", command=lambda: self.toolkit.input_class.select_pose_then_launch_fixpdb())
        self.file_menu.add_command(label="Fetch PDB", command = lambda: self.toolkit.input_class.fetch_pdb())
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Reload PDB", command = lambda: self.toolkit.input_class.load_pose(self.toolkit.input_class.pdb_path.get()))
        self.file_menu.add_separator()
        self.file_menu.add_cascade(label="Import", menu=self.import_menu)
        self.file_menu.add_cascade(label="Export", menu=self.export_menu)
        self.file_menu.add_separator()
        self.file_menu.add_command(label ="Setup PDB for Rosetta", command=lambda: self.show_fxpdb_window())

        ## Note:  There is no disable for menu items.  Which is why we either show the item or not, instead of disabling it.
        if flag_file_builder_import:
            self.file_menu.add_command(label= "Rosetta Flag File Builder", command = lambda: self.show_RosettaProtocolBuilder())

        self.main_menu.add_cascade(label="File", menu=self.file_menu)
    def _set_options_menu(self):
        self.options_menu = Menu(self.main_menu, tearoff=0)

        self.rounds_options_menu=Menu(self.main_menu, tearoff=0)

        #Could be it's own window in the future.
        self.rounds_options_menu.add_command(label = "Set processors to use", command = lambda: self.toolkit.output_class.processors.set(tkSimpleDialog.askinteger(title="Processesors", prompt= "Please set the number of processess you wish to create for protocol runs.", initialvalue=self.toolkit.output_class.processors.get())))
        self.rounds_options_menu.add_checkbutton(label = "Use Boltzmann criterion each round", variable=self.toolkit.output_class.use_boltzmann)
        self.rounds_options_menu.add_command(label = "Set Temperature", command = lambda: self.toolkit.output_class.set_temperature())
        self.rounds_options_menu.add_checkbutton(label = "Recover lowest energy from all rounds", variable = self.toolkit.output_class.recover_low)
        self.options_menu.add_command(label="Configure Option System",command = lambda: self.show_OptionsSystemManager())
        self.options_menu.add_command(label="Configure Score Function", command =lambda: self.toolkit.score_class.makeWindow(Toplevel(self.main), self.toolkit.pose))
        self.options_menu.add_cascade(label="Set General Protocol Options", menu = self.rounds_options_menu)
        self.options_menu.add_separator()

	self.options_menu.add_checkbutton(label="Set output to terminal", variable = self.toolkit.output_class.terminal_output)


	self.main_menu.add_cascade(label="Options", menu=self.options_menu)

    def _set_visualization_menu(self):
	self.vis_menu = Menu(self.main_menu, tearoff=0)
	self.vis_menu.add_command(label = "Show Pose in PyMOL", command = lambda: self.toolkit.pymol_class.pymover.apply(self.toolkit.pose))
	#self.vis_menu.add_command(label = "Show current regions")
	self.vis_menu.add_command(label = "Show current region", command = lambda: self.toolkit.input_frame.color_region_from_entry())
	self.vis_menu.add_command(label = "Show current residue", command = lambda: self.toolkit.input_frame.color_residue())

	self.vis_observer_menu = Menu(self.main_menu, tearoff=0)
	self.vis_observer_menu.add_checkbutton(label="Set Pymol Observer", variable=self.toolkit.pymol_class.auto_send)
	self.vis_observer_menu.add_separator()
	self.vis_observer_menu.add_checkbutton(label="Send as new states", variable=self.toolkit.pymol_class.keep_history)
	self.vis_observer_menu.add_checkbutton(label="Send energies", variable=self.toolkit.pymol_class.send_energies)



	self.vis_menu.add_separator()

	self.vis_menu.add_checkbutton(label = "Color current region on addition", variable=self.toolkit.pymol_class.auto_send_region_colors)
	self.vis_menu.add_checkbutton(label = "Color current residue on selection", variable = self.toolkit.pymol_class.auto_send_residue_colors)

	self.vis_menu.add_separator()
	self.vis_menu.add_cascade(label = "PyMol Observer", menu=self.vis_observer_menu)
	self.vis_menu.add_command(label="Advanced PyMOL Visualization", command=lambda: self.toolkit.pymol_class.makeWindow(0, 0, Toplevel(self.main), self.toolkit.score_class))

	self.main_menu.add_cascade(label="Visualization", menu=self.vis_menu)

    def _set_advanced_menu(self):
	"""
	Sets the advanced control menu
	"""

	self.advanced_menu=Menu(self.main_menu, tearoff=0)

	###Analysis###
	self.analysis_menu = Menu(self.main_menu, tearoff=0)

	self.analysis_menu.add_command(label = "Interface Analyzer", command=lambda: analysis_tools.analyze_interface(self.toolkit.pose, self.toolkit.score_class.score, self.toolkit))
	self.analysis_menu.add_command(label = "Packing Analyzer", command = lambda: analysis_tools.analyze_packing(self.toolkit.pose))
	self.analysis_menu.add_command(label = "Loops Analyzer", command = lambda: analysis_tools.analyze_loops(self.toolkit.pose, self.toolkit.input_class.loops_as_strings))
	self.analysis_menu.add_command(label = "VIP Analyzer", command = lambda: AnalysisProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class).analyze_vip())
	self.analysis_menu.add_separator()
	self.analysis_menu.add_command(label = "Insert Data into B-Factor", command=lambda:InsertBFactor(self.toolkit.input_class).show_window(self.main))
	self.advanced_menu.add_cascade(label = "Analysis", menu = self.analysis_menu)

	###Design###
	self.design_menu=Menu(self.main_menu, tearoff=0)
	self.design_menu.add_command(label="Setup Resfile for PDB", command=lambda: self.show_ResfileDesignWindow())
	self.design_menu.add_command(label="Setup Blueprint for PDB", foreground='red')
	#self.design_menu.add_command(label="Structure Editing and Grafting", foreground='red')
	self.advanced_menu.add_cascade(label = "Design", menu=self.design_menu)
	self.advanced_menu.add_separator()
	self.advanced_menu.add_command(label ="Add Constraints", command = lambda: self.toolkit.input_class.constraint_file_paths.append(\
				    input_tools.add_constraints_to_pose_and_scorefunction(self.toolkit.pose, self.toolkit.score_class.score)))
	#self.advanced_menu.add_command(label ="Enable Symmetry", foreground='red')

	###NonStandard###
	self.non_standard_menu = Menu(self.main_menu, tearoff=0)
	self.non_standard_menu.add_command(label = "Ligand/NCAA/PTM Manager", command = lambda: self.show_ligand_ncaa_ptm_manager())
	self.non_standard_menu.add_separator()
	self.non_standard_menu.add_command(label = "Download Rosetta NCAA Rotamer Library", command=lambda: webbrowser.open("http://carl.bio.nyu.edu/~renfrew/ncaa/"))
	self.non_standard_menu.add_command(label = "Convert molfile to param files", command = lambda: output_tools.output_molfile_to_params())
	self.non_standard_menu.add_separator()
	self.non_standard_menu.add_command(label = "Documentation", command = lambda: webbrowser.open("http://www.pyrosetta.org/obtaining-and-preparing-ligand-pdb-files"))

	self.advanced_menu.add_cascade(label ="Enable NCAA/PTM/Ligands", menu=self.non_standard_menu)
	self.advanced_menu.add_separator()


	self.advanced_menu.add_command(label="Per Residue Control and Analysis", command=lambda: self.show_fullcontrol_window())
	#self.advanced_menu.add_command(label="Extract PDB from SQLite3 DB", command = lambda: output_tools.extract_pdb_from_sqlite3db())
	#self.advanced_menu.add_command(label="Interactive Terminal", foreground='red',command = lambda: self.show_IpythonWindow())
	#self.advanced_menu.add_command(label="Jump into Session", foreground='red', command = lambda: embed())
	self.main_menu.add_cascade(label = "Advanced", menu = self.advanced_menu)

    def _set_protocols_menu(self):
	"""
	Menu for running protocols through PyRosetta.
	"""

	self.protocols_menu = Menu(self.main_menu, tearoff=0)

	#self.protocols_menu.add_command(label = "Enable MPI Mode", foreground = 'red')

	self.protocols_menu.add_separator()

	#Setup Protocol Classes. Should this go to Main File? They are very light weight classes.
	self.design_class = DesignProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.docking_class = DockingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.low_res_loop_modeling_class = LowResLoopModelingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.high_res_loop_modeling_class = HighResLoopModelingProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.floppy_tail_class = FloppyTailProtocol(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
	self.minimize_class = MinimizationProtocols(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)

	#Design

	self.design_protocols = Menu(self.main_menu, tearoff=0)
	self.design_protocols.add_command(label = "FixedBB", command = lambda: self.design_class.setupPackDesign(self.main))
	self.design_protocols.add_command(label = "Grafting", command = lambda:self.show_graftmover_window())
	#self.design_protocols.add_command(label = "Remodel", foreground='red')

	#Docking

	self.docking_protocols = Menu(self.main_menu, tearoff=0)
	self.docking_protocols.add_command(label = "Low Resolution", command = lambda: self.docking_class.low_res_dock())
	self.docking_protocols.add_command(label = "High Resolution", command = lambda: self.docking_class.high_res_dock())
	#self.docking_protocols.add_command(label = "Ligand", foreground='red')

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

	#Minimization.  Just so it's here too.
	self.minimization_protocols = Menu(self.main_menu, tearoff=0)
	self.minimization_protocols.add_command(label = "RePack Region Rotamers", command = lambda:self.minimize_class.optimize_rotamers())
	self.minimization_protocols.add_command(label = "RePack Region Rotamers (SCWRL)", command = lambda:self.minimize_class.SCWRL())
	self.minimization_protocols.add_separator()
	self.minimization_protocols.add_command(label = "Minimize Backbone+Sidechains", command = lambda:self.minimize_class.minimize())
	self.minimization_protocols.add_command(label = "Minimize Backbone only", command = lambda:self.minimize_class.minimize(False, False, True))
	self.minimization_protocols.add_command(label = "Minimize Sidechains only", command = lambda:self.minimize_class.minimize(False, False, False, True))
	self.minimization_protocols.add_separator()
	self.minimization_protocols.add_command(label = "FastRelax Backbone+Sidechains", command = lambda:self.minimize_class.relax(1))
	self.minimization_protocols.add_command(label = "FastRelax Backbone only", command = lambda:self.minimize_class.relax(1, False, True))
	self.minimization_protocols.add_command(label = "FastRelax Sidechains only", command = lambda:self.minimize_class.relax(1, False, False, True))
	self.minimization_protocols.add_separator()
	self.minimization_protocols.add_command(label = "ClassicRelax Backbone+Sidechains", command = lambda:self.minimize_class.relax(0))
	self.minimization_protocols.add_command(label = "ClassicRelax Backbone only", command = lambda:self.minimize_class.relax(0, False, True))
	self.minimization_protocols.add_command(label = "ClassicRelax Sidechains only", command = lambda:self.minimize_class.relax(0, False, False, True))

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
	self.servers.add_command(label = "Scaffold Select", command = lambda: webbrowser.open("http://rosettadesign.med.unc.edu/scaffold/"))
	self.rna_protocols = Menu(self.main_menu, tearoff=0)

	self.dna_protocols = Menu(self.main_menu, tearoff=0)

	self.general_protocols.add_cascade(label = "Web Servers", menu = self.servers)
	self.general_protocols.add_separator()
	#self.general_protocols.add_cascade(label = "RNA", menu = self.rna_protocols)
	#self.general_protocols.add_cascade(label = "DNA", menu = self.dna_protocols)
	self.general_protocols.add_command(label = "FloppyTail", command=lambda: self.floppy_tail_class.setup_protocol(self.main))
	#self.general_protocols.add_command(label = "Ab Initio", foreground='red')
	#self.general_protocols.add_command(label = "Homology Modeling", foreground='red')
	self.protocols_menu.add_cascade(label = "General", menu = self.general_protocols)
	self.protocols_menu.add_separator()

	self.protocols_menu.add_cascade(label = "Loop Modeling", menu = self.loop_modeling_protocols)
	self.protocols_menu.add_cascade(label = "Docking", menu = self.docking_protocols)
	self.protocols_menu.add_cascade(label = "Design", menu = self.design_protocols)
	#self.protocols_menu.add_separator()
	self.protocols_menu.add_cascade(label = "Minimization", menu=self.minimization_protocols)

	self.main_menu.add_cascade(label = "Protocols", menu=self.protocols_menu)

    def _set_pdblist_menu(self):
	"""
	Sets a Menu for interacting with a simple PDBList.  Useful since Rosetta is pretty much all about thousands of models.
	"""

	self.pdblist_tools_menu = Menu(self.main_menu, tearoff=0)
	#There are scripts in RosettaTools that do this welll...
	self.pdblist_tools_menu.add_command(label = "Load PDBList", command=lambda: self.toolkit.input_class.set_PDBLIST())
	self.pdblist_tools_menu.add_command(label = "Create PDBList + Load", command = lambda: self.toolkit.input_class.PDBLIST.set(output_tools.make_PDBLIST()))
	self.pdblist_tools_menu.add_command(label = "Create PDBList Recursively + Load", command = lambda: self.toolkit.input_class.PDBLIST.set(output_tools.make_PDBLIST_recursively()))
	self.pdblist_tools_menu.add_separator()
	self.score_menu = Menu(self.main_menu, tearoff=0)
	self.score_menu.add_command(label = "Rescore PDBList + Output ScoredPDBList + Load", command = lambda: self.score_analyzer.set_filepath(output_tools.score_PDBLIST(self.toolkit.input_class.PDBLIST.get(), self.toolkit.score_class.score, self.toolkit.output_class.processors.get(), self.toolkit.output_class)))
	self.score_menu.add_command(label = "Load Scores (.sc, .fasc, or ScoredPDBList)", command = lambda: self.load_scores_for_score_analysis())
	self.score_menu.add_separator()
    	self.score_menu.add_command(label = "Get top model", command = lambda: self.score_analyzer.get_top_scoring())
	self.score_menu.add_command(label = "Get top % models", command = lambda: self.score_analyzer.get_top_scoring_by_percent())
	self.score_menu.add_command(label = "Get top # models", command = lambda: self.score_analyzer.get_top_scoring_by_number())
	self.score_menu.add_separator()
	self.score_menu.add_command(label = "Energy vs RMSD", command=lambda:self.score_analyzer.get_score_vs_rmsd(self.toolkit.pose, self.toolkit.input_class.loops_as_strings))
	self.score_menu.add_separator()

	#self.plot_menu = Menu(self.main_menu, tearoff=0)
	#self.plot_menu.add_command(label = "Plot Score KDE", command = lambda:self.score_analyzer.plot_score_kde(), foreground='red')
	#self.plot_menu.add_command(label = "Plot Score vs RMSD", command = lambda: self.score_analyzer.plot_score_vs_rmsd(), foreground='red')
	#self.plot_menu.add_command(label = "Plot Score vs RMSD KDE", command = lambda: self.score_analyzer.plot_score_vs_rmsd_KDE(), foreground='red')
	#self.score_menu.add_cascade(label = "R Plots", menu = self.plot_menu)

	self.pdblist_tools_menu.add_cascade(label = "Score Analysis", menu=self.score_menu)

	self.sequence_menu = Menu(self.main_menu, tearoff=0)
	self.sequence_menu.add_command(label = "Output FASTA for Each PDB", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.input_class.PDBLIST.get(), False))
	self.sequence_menu.add_command(label = "Output FASTA for Each Region", command = lambda: output_tools.save_FASTA_PDBLIST(self.toolkit.input_class.PDBLIST.get(), False, self.toolkit.input_class.regions))
	#If you need this, you know how to program: self.sequence_menu.add_command(label = "Output LOOP file for each PDB", command = lambda: output_tools.save_LOOP_PDBLIST(self.toolkit.input_class.PDBLIST.get()))
	self.sequence_menu.add_command(label = "Use FASTA for design statistics", command = lambda: self.run_design_breakdown())
	self.pdblist_tools_menu.add_cascade(label = "Sequence Analysis", menu=self.sequence_menu)
	self.pdblist_tools_menu.add_separator()

	#self.pdblist_tools_menu.add_command(label = "Cluster PDBList using Calibur", foreground='red')
	#self.pdblist_tools_menu.add_command(label = "Convert PDBList to SQLite3 DB", command = lambda: output_tools.convert_PDBLIST_to_sqlite3db(self.toolkit.input_class.PDBLIST.get()))
	#self.pdblist_tools_menu.add_command(label = "Extract PDBList from SQLite3 DB", command = lambda: output_tools.extract_pdbs_from_sqlite3db(self.toolkit.input_class.PDBLIST.get()))
	#self.pdblist_tools_menu.add_separator()
	self.pdblist_tools_menu.add_command(label = "Rename All PDBs Recursively + Copy to Outpath", command = lambda: output_tools.rename_and_save(self.toolkit.input_class.PDBLIST.get()))

	self.main_menu.add_cascade(label = "PDBLists", menu=self.pdblist_tools_menu)

    def _set_help_menu(self):
	"""
	Sets the help Menu.  TK kinda sucks for formating dialog windows.  Just FYI.
	"""

	self.help_menu=Menu(self.main_menu, tearoff=0)
	self.help_menu.add_command(label="About", command=lambda: help_tools.about())
	self.help_menu.add_command(label = "License", command = lambda: help_tools.show_license())
	self.help_menu.add_separator()
	self.help_menu.add_command(label = "Region Selection", command = lambda: help_tools.region_selection())
	self.help_menu.add_command(label = "PyMOL Setup", command = lambda:webbrowser.open("http://www.pyrosetta.org/pymol_mover-tutorial"))
	self.help_menu.add_separator()
	self.help_menu.add_command(label="Rosetta Manual", command = lambda: webbrowser.open("http://www.rosettacommons.org/manual_guide"))
	self.help_menu.add_command(label="Rosetta Glossary", command=lambda: help_tools.print_glossary())
	self.help_menu.add_command(label="Rosetta Forums", command = lambda: webbrowser.open("http://www.rosettacommons.org/forum"))
	self.help_menu.add_command(label="Rosetta BugTracker", command = lambda: webbrowser.open("http://bugs.rosettacommons.org"))
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
	self.help_menu.add_separator()
	self.help_menu_scwrl = Menu(self.main_menu, tearoff=0)
	self.help_menu_scwrl.add_command(label = "Download", command = lambda: webbrowser.open("http://dunbrack.fccc.edu/scwrl4/SCWRL4.php"))
	self.help_menu_scwrl.add_command(label = "Install", command = lambda: webbrowser.open("http://dunbrack.fccc.edu/scwrl4/SCWRL4.php#installation"))
	self.help_menu_scwrl.add_command(label = "GUI Install", command = lambda: help_tools.scwrl())
	self.help_menu.add_cascade(label = "SCWRL", menu=self.help_menu_scwrl)
	self.main_menu.add_cascade(label="Help", menu=self.help_menu)


##### MENU FUNCTIONS #######

    def load_scores_for_score_analysis(self):
        filename = tkFileDialog.askopenfilename(title="Open ScoredPDBList, .fasc, or .sc file", initialdir=global_variables.current_directory)
        if not filename:return
        global_variables.current_directory = os.path.dirname(filename)
        self.score_analyzer.set_filepath(filename)

    #These are bullshit.
    def get_path(self, string):
        path = tkFileDialog.askopenfilename(title=string,initialdir = global_variables.current_directory)
        return path

    def run_design_breakdown(self):

    #Super huge bug in Tkinter, which is not allowing multiple dialog boxes to be called in succession. ~jadolfbr
    #Jan 2013 MacBookPro~2010 OS 10.6, Python 2.6.

        if self.toolkit.pose.total_residue()==0:
            print "Please load a pose for reference."
            return

        fasta_path = self.get_path("FASTA Path")
        if not fasta_path:return
        global_variables.current_directory = os.path.dirname(fasta_path)
        outpath = os.path.dirname(fasta_path)+"/RESULTS"



        #else:
            #answer = tkMessageBox.askokcancel(title = "Current", message = "Using current pose as reference.  Continue?")
            #if not answer: return
        reference_path=self.toolkit.input_class.pdb_path.get()
        if not reference_path:return
        #native = tkMessageBox.askyesno(title="Native", message="Use pose as Native?")

        breakdown = DesignBreakdown(fasta_path, reference_path,  outpath)
        breakdown.run_outputs()

#### WINDOWS ##### (ADD NEW WINDOWS TO THIS THAT NEED TO BE SET UP) #######
    def show_fullcontrol_window(self):
        self.fullcontrol_class = FullControlWindow(self.toolkit.score_class, self.toolkit.pose, self.toolkit.input_class, self.toolkit.output_class);
        self.fullcontrol_class.show_window(self.main, 0, 0)

    def show_graftmover_window(self):
        grafter = GraftMoverWindow(self.toolkit.pose, self.toolkit.score_class, self.toolkit.input_class, self.toolkit.output_class)
        top_level_tk = Toplevel(self.main)
        grafter.setTk(top_level_tk)
        grafter.shoTk(0,0)

    def show_fxpdb_window(self):
        cleaner = FixPDBWindow(self.toolkit.input_class, self.toolkit.score_class, self.toolkit.pose)
        cleaner.runfixPDBWindow(self.main, 0, 0)

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
        self.toolkit.input_class.options_manager.setTk(top_level_tk)
        self.toolkit.input_class.options_manager.shoTk()

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
        rosetta_protocol_builder = RosettaFlagFileBuilder(top_level_tk)
        if not rosetta_protocol_builder.result:return
        rosetta_protocol_builder.setTk()
        rosetta_protocol_builder.shoTk(0, 0)
        rosetta_protocol_builder.setMenu(top_level_tk)

    def show_IpythonWindow(self):
        """
        IPython Interactive Window.  Isolated from variables for now.
        """
        term = IPythonView(Toplevel(self.main))
        term.pack()
