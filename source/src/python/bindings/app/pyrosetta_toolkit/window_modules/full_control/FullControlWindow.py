#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/full_control/FullControl.py
## @brief  Full Control window
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.pose import remove_variant_type_from_pose_residue

#Python Imports
import glob

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog

#Toolkit Imports
import app.pyrosetta_toolkit.modules.tools.analysis as analysis_tools
from app.pyrosetta_toolkit.modules.definitions.restype_definitions import *
from app.pyrosetta_toolkit.modules.tools import general_tools
from app.pyrosetta_toolkit.modules.tools import protocols as protocol_tools
from app.pyrosetta_toolkit.modules.protocols.DesignProtocols import DesignProtocols
from app.pyrosetta_toolkit.modules.protocols.MinimizationProtocols import MinimizationProtocols
from app.pyrosetta_toolkit.modules.RegionalScoring import *
from app.pyrosetta_toolkit.window_modules.ligand_ncaa_ptm_manager.ligand_ncaa_ptm_manager import ligand_ncaa_ptm_manager
from app.pyrosetta_toolkit.window_modules.scorefunction.ScoreFxnControl import ScoreFxnControl
#from window_main.IO.GUIInput import GUIInput
#from window_main.IO.GUIOutput import GUIOutput

class FullControlWindow():
    """
    Window for full control of protein backbone and sidechains
    """

    def __init__(self, ScoreObject, pose, input_class, output_class):

	#Main classes
	self.pose = pose
	self.input_class = input_class
	self.design_protocols = DesignProtocols(self.pose, ScoreObject, input_class, output_class)
	self.loop_protocols = MinimizationProtocols(self.pose, ScoreObject, input_class, output_class)
	self.residue_definitions = definitions(); #Defines all residue type information.
	self.score_object = ScoreObject
	self.variant_manager = ligand_ncaa_ptm_manager(input_class, ScoreObject, pose)


	#GUI variables
	self.resnum = StringVar()
	self.chain = StringVar()
	
	
	self.PhiVar = StringVar(); #phi string
	self.PsiVar = StringVar(); #psi string
	self.OmeVar = StringVar(); #omega string
	self.restype_var = StringVar(); #restype string
	self.rotamer_probabilty = StringVar()
	self.rotamer_energy = StringVar()
	self.eterm_type = StringVar()
	self.eterm_type.set("Residue Energy: ")
	self.residue_energy_total = StringVar()
	self.residue_energy_of_Eterm = StringVar(); #Energy of the particular residue

	self.variant = StringVar()
	self.variant_map = dict(); # [string variant]:[string names]
        
	
        #Ignore this.  It is for Komodo autocomplete
        if 0:
            self.score_object = ScoreFxnControl()
	    self.input_class = GUIInput()
		
    def __exit__(self):
	exit()

    def shoInfo(self, res, chain):
	if self.pose.total_residue()==0:return
	if not res:return
	if not chain:return
	
	#self.populate_restype_listbox()
	self.populate_energy_term_listbox()
	self.regional_scoring = RegionalScoring(self.pose, self.score_object.score); #Happens every time just in case the user selects a dif scorefunction
	
	#Much exception handling to make sure no errors are thrown for particular residues, numbering inconsistancies, DNA/ligands
	
	#Numbering Check
	try:
	    res = self.pose.pdb_info().pdb2pose(chain, int(res))
            
            if res==0:
                print "Residue does not exist in PDB"
                return
	    self.restype_var.set(self.pose.residue(res).name()); #Exception comes from this! Segfualts on .psi etc.  This needs to be first.
	    self.input_class.residue_rosetta_resnum.set(repr(res))
	except PyRosettaException:
	    print "Residue does not exist in PDB.."
	    return
	
	#BackBone Check
	try:
	    self.PhiVar.set("%.3f"%self.pose.phi(res)); self.PsiVar.set("%.3f"%self.pose.psi(res)); self.OmeVar.set("%.3f"%self.pose.omega(res))
	except PyRosettaException:
	    print "Residue does not have phi, psi or omega"
	    
	
	#Rotamer Check
	try:
	    self.rotamer_energy.set("%.3f REU"%analysis_tools.return_energy(self.pose, res))
	    self.rotamer_probabilty.set("%.3f"%analysis_tools.return_probability(self.pose, res))
	except PyRosettaException:
	    print "Residue does not have a rotamer energy"
	    
	#Calculate Total Energy.
	self.residue_energy_total.set("%.3f REU"%self.regional_scoring.ret_total_weighted_residue_energy(res))
	if not self.eterm_type.get()=="Residue Energy: ":
	    #get_residue_energy function should split so it accepts the energy type directly.
	    self.get_residue_energy(self.last_type_string_selected)

    def nextRes(self, res, chain):
	res = int(res)+1
	self.resnum.set(res)
	self.shoInfo(res, chain)
    def prevRes(self, res, chain):
	res = int(res)-1
	self.resnum.set(res)
	self.shoInfo(res, chain)

    def show_window(self, master, r=0, c=0):
	
	try :
	    print self.pose.pdb_info().name()
	except AttributeError:
	    tkMessageBox.showwarning(message = 'Please Load a Pose...')
	    return

	try:
	    self.regional_scoring = RegionalScoring(self.pose, self.score_object.score)

	except PyRosettaException:
	    print "Please Load a pose"
	    return
	
    
	self.main = Toplevel(master); #This is to prevent Ghosting of windows.
	self.main.title("Full Control")
	self.main.grid_columnconfigure(ALL, weight=1)
	self.entry_Res = Entry(self.main, textvariable=self.resnum, justify=CENTER)
	self.entry_Cha = Entry(self.main, textvariable=self.chain, justify=CENTER)
	self.entry_Nam = Entry(self.main, textvariable=self.restype_var, justify=CENTER)
	self.label_residue = Label(self.main, text = "Residue #")
	self.label_chain = Label(self.main, text = "Chain ID")
	self.label_restype = Label(self.main, text = "Residue Name")
	self.button_Get = Button(self.main, text="Get Info", command=lambda: self.shoInfo(self.resnum.get(), self.chain.get()))
	self.button_previous_residue = Button(self.main, text="Previous Residue", command=lambda: self.prevRes(self.resnum.get(), self.chain.get()))
	self.button_next_residue = Button(self.main, text = "Next Residue", command=lambda: self.nextRes(self.resnum.get(), self.chain.get()))
	self.entry_Phi = Entry(self.main, textvariable=self.PhiVar, justify=CENTER)
	self.entry_Psi = Entry(self.main, textvariable=self.PsiVar, justify=CENTER)
	self.entry_Ome = Entry(self.main, textvariable=self.OmeVar, justify=CENTER)
	self.label_Phi = Label(self.main, text="Phi")
	self.label_Psi = Label(self.main, text="Psi")
	self.label_Ome = Label(self.main, text="Omega")
	self.button_Delta = Button(self.main, text = "Change Dihedrals", command=lambda: self.pose.assign(protocol_tools.bbMover(self.pose, int(self.resnum.get()), self.chain.get(), float(self.PhiVar.get()), float(self.PsiVar.get()), float(self.OmeVar.get()), self.score_object.score)))
	self.label_rotamer_probability = Label(self.main, text = "Aprox Rotamer Probability")
	self.label_rotamer_energy = Label(self.main, text = "Rotamer Energy")
	self.entry_rotamer_probability = Label(self.main, textvariable=self.rotamer_probabilty, justify=CENTER, relief=SUNKEN)
	self.entry_rotamer_energy = Label(self.main, textvariable=self.rotamer_energy, justify=CENTER, relief=SUNKEN)
	self.slider_phi = Scale(self.main, from_=-180, to=180, orient=HORIZONTAL,variable = self.PhiVar)
	self.slider_psi = Scale(self.main, from_=-180, to=180, orient=HORIZONTAL, variable = self.PsiVar)


	#### Quick Protocols ####
	self.button_relax_residue = Button(self.main, text = "Relax Residue", command=lambda: self.relax_residue())
	self.button_backrub_residue= Button(self.main, text = "Backrub Residue", state=DISABLED)
	self.button_relax_residue_and_neighbors=Button(self.main, text= "Relax X-Res-X", command = lambda: self.relax_residue_neighbors(False))
	self.button_relax_residue_and_neighbors_bb=Button(self.main, text= "Relax X-Res-X (BB only)", command = lambda: self.relax_residue_neighbors(True))
	self.button_backrub_residue_and_neighbors=Button(self.main, text = "Backrub X-Res-X", state=DISABLED)
	self.button_Pack = Button(self.main, text = "Pack Rotamer", command=lambda: self.packRes())
	self.button_full_design = Button(self.main, text="Design", command = lambda: self.desRes())
	self.listbox_mutate_restypes = Listbox(self.main)
	self.variant_listbox = Listbox(self.main)
	self.add_variant_button = Button(self.main, text = "Add Variant", command = lambda: self.mutate_to_variant())
	self.remove_variant_button = Button(self.main, text = "Remove Variant", command = lambda: self.remove_variant())
	self.label_Etot = Label(self.main, text = "Total Residue Energy")
	self.entry_Etot = Label(self.main, textvariable = self.residue_energy_total, relief=SUNKEN)
	self.label_Terms = Label(self.main, text = "energy terms")
	self.entry_ResE= Label(self.main, textvariable = self.residue_energy_of_Eterm, relief=SUNKEN)
	self.label_Etype = Label(self.main, textvariable = self.eterm_type)
	self.listbox_energy_terms = Listbox(self.main)



	#### GRID ####
	self.entry_Res.grid(row=r+5, column=c+0); self.entry_Cha.grid(row=r+5, column=c+1); self.entry_Nam.grid(row=r+5, column=c+2)
	self.label_residue.grid(row=r+6, column=c+0); self.label_chain.grid(row=r+6, column=c+1); self.label_restype.grid(row=r+6,column=c+2)
	self.button_previous_residue.grid(row=r+7, column=c+0); self.button_Get.grid(row=r+7, column=c+1, sticky=W+E); self.button_next_residue.grid(row=r+7, column=c+2)
	self.entry_Phi.grid(row=r+8, column=c+0); self.entry_Psi.grid(row=r+8, column=c+1); self.entry_Ome.grid(row=r+8, column=c+2)
	self.slider_phi.grid(row=r+9, column=c+0); self.slider_psi.grid(row=r+9, column=c+1)
	self.label_Phi.grid(row=r+10, column=c+0); self.label_Psi.grid(row=r+10, column=c+1); self.label_Ome.grid(row=r+10, column=c+2)
	self.button_Delta.grid(row=r+11, column=c+1);
	self.entry_rotamer_probability.grid(row=r+12, column=c+0); self.entry_rotamer_energy.grid(row=r+12, column=c+2)
	self.label_rotamer_probability.grid(row=r+13, column=c+0); self.label_rotamer_energy.grid(row=r+13, column=c+2)

	#### Mutagenesis, Repacking, Residue Energies ####
	self.button_full_design.grid(row=r+13, column=c+1)
	self.button_Pack.grid(row=r+15, column=c+1)
	self.listbox_mutate_restypes.grid(row=r+17, column=c+1, rowspan=6)
	self.variant_listbox.grid(row = r+23, column=c+1, rowspan=6)
	self.remove_variant_button.grid(row=r+23, column=c); self.add_variant_button.grid(row=r+23, column=c+2)
	self.entry_Etot.grid(row=r+29, column=c+0); self.entry_ResE.grid(row=r+29, column=c+2)
	self.label_Etot.grid(row=r+30, column=c+0); self.label_Terms.grid(row=r+30, column=c+1); self.label_Etype.grid(row=r+30, column=c+2)
	self.listbox_energy_terms.grid(row=r+31, column = 1, rowspan=6)
	self.listbox_mutate_restypes.bind("<Double-Button-1>", lambda event: self.mutRes())
	self.variant_listbox.bind("<ButtonRelease-1>", lambda event: self.show_variant_info())
	self.variant_listbox.bind("<Double-Button-1>", lambda event: self.mutate_to_variant())
	self.listbox_energy_terms.bind("<Double-Button-1>", lambda event: self.get_residue_energy(self.listbox_energy_terms.get(self.listbox_energy_terms.curselection())))
	self.main.grid_columnconfigure(ALL, weight=1)
	
	#### Quick Min ####
	self.button_relax_residue.grid(row=r+17, column=c+2, sticky= W+E)
	self.button_backrub_residue.grid(row=r+17, column=c, sticky= W+E)
	self.button_relax_residue_and_neighbors.grid(row=r+18, column=c+2, sticky= W+E)
	#self.button_relax_residue_and_neighbors_bb.grid(row=r+18, column=c+2, sticky= W+E)
	self.button_backrub_residue_and_neighbors.grid(row=r+18, column=c, sticky= W+E)

	#Populating the Listboxes:
	self.populate_restype_listbox()
	self.populate_energy_term_listbox()
	self.populate_variant_listbox()
	
	#These are to hook up the sub window with the sequence entry from the main window.
	self.main_resnum = self.input_class.residue_resnum
	self.main_chain = self.input_class.residue_chain
	self.main_resnum.trace_variable('w', self.resnum_callback)
	self.main_chain.trace_variable('w', self.resnum_callback)
	self.resnum.set(self.main_resnum.get())
	self.chain.set(self.main_chain.get())
	self.shoInfo(self.main_resnum.get(), self.main_chain.get())
    
    def resnum_callback(self, name, index, mode):
	"""
	This is to update GUI upon a new resnum selected from main window.
	"""
	self.shoInfo(self.main_resnum.get(), self.main_chain.get())
	self.resnum.set(self.main_resnum.get())
	self.chain.set(self.main_chain.get())
	
    def get_residue_energy(self, type_weight_string):
	self.last_type_string_selected = type_weight_string
	type_string = type_weight_string.split(";")[0]
	weight = float(type_weight_string.split(";")[1])
	res = int(self.resnum.get())
	chain = self.chain.get()
	posenum = self.pose.pdb_info().pdb2pose(chain, res)
	emap = self.regional_scoring.ret_residue_energy(posenum)
	
	#e = weight*emap[self.regional_scoring.fa_terms[type_string]]
	e = weight*emap[rosetta.core.scoring.score_type_from_name(type_string)]
	self.eterm_type.set(type_string)
	self.residue_energy_of_Eterm.set("%.3f REU"%e)

    def populate_restype_listbox(self):
	self.listbox_mutate_restypes.delete(0, END)
	for type in self.residue_definitions.restype_info["All"]:
	    self.listbox_mutate_restypes.insert(END, type)
    
    def populate_energy_term_listbox(self):
	self.listbox_energy_terms.delete(0, END)
	
	(zeroscores, nonzeroscores) = self.score_object.scoreOption("Breakdown ScoreFxn")
	for type in nonzeroscores:
	    self.listbox_energy_terms.insert(END, type)
	    #For now, I will show this.  Later, I will use an object to hold all Residue and Design/CDR info.  Useful in the future.
    
    def populate_variant_listbox(self):
	self.variant_listbox.delete(0, END)

	#Patches first
	for variant in sorted(self.variant_manager.patch_type_map):
	    if not self.variant_map.has_key(variant):
		self.variant_map[variant]=[]
		for type in self.variant_manager.patch_type_map[variant]:
		    self.variant_map[variant].append(type)

	"""
	#Now params - Adding variants from params -  Protonated, deprotonated, doesn't seem to work.  So, we are commmenting this out for now.
	for name in self.variant_manager.param_name_map:
	    p = self.variant_manager.param_name_map[name]

	    if p.variant.get():
		variant = p.variant.get().lower()
		if not self.variant_map.has_key(variant):
		    self.variant_map[variant]=[]
		    self.variant_map[variant].append(p.name.get())
		else:
		    self.variant_map[variant].append(p.name.get())
	"""

	for variant in sorted(self.variant_map):
	    self.variant_listbox.insert(END, variant)

    def show_variant_info(self):
	variant = self.variant_listbox.get(self.variant_listbox.curselection())
	p = ""
	for name in self.variant_map[variant]:
	    p = p+" "+name+","

	print "Types: "+p
	self.variant.set(variant)


    def mutate_to_variant(self):
	variant = self.variant_listbox.get(self.variant_listbox.curselection())

	self.variant_manager.current_variant.set(variant.upper())
	residue = self.pose.pdb_info().pdb2pose(self.chain.get(), int(self.resnum.get()))
	try:
	    self.variant_manager.current_main_selection.set("patch")
	except AttributeError:
	    pass

	self.variant_manager.mutate(residue)
	self.shoInfo(self.resnum.get(), self.chain.get())

    def remove_variant(self):
	resnum = self.pose.pdb_info().pdb2pose(self.chain.get(), int(self.resnum.get()))
	remove_variant_type_from_pose_residue(self.pose, self.variant.get().upper(), resnum)
	self.shoInfo(self.resnum.get(), self.chain.get())

    #### Functions ####
    def relax_residue(self):
	"""
	Relaxes the given residue
	"""

	self.loop_protocols.relax_residue_and_neighbors(1, (self.resnum.get()), self.chain.get(), 0)
	self.shoInfo(self.resnum.get(), self.chain.get())

    def relax_residue_neighbors(self, bbonly=False):
	"""
	Relaxes residues and neighbors
	"""

	if bbonly:
	    self.loop_protocols.relax_residue_and_neighbors(1, (self.resnum.get()), self.chain.get(), 1, True)
	else:
	    self.loop_protocols.relax_residue_and_neighbors(1, (self.resnum.get()), self.chain.get())
	self.shoInfo(self.resnum.get(), self.chain.get())

    def packRes(self):
	self.design_protocols.pack_residue(self.resnum.get(), self.chain.get())
	self.shoInfo(self.resnum.get(), self.chain.get())

    def mutRes(self):
	self.design_protocols.mutateRes(self.resnum.get(), self.chain.get(), self.listbox_mutate_restypes.get(self.listbox_mutate_restypes.curselection()))
	self.shoInfo(self.resnum.get(), self.chain.get())

    def desRes(self):
	self.design_protocols.design_residue(self.resnum.get(), self.chain.get())
	self.shoInfo(self.resnum.get(), self.chain.get())
