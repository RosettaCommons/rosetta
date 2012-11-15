#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/full_control/FullControl.py
## @brief  Full Control window
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
import glob
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
import modules.tools.analysis as analysis_tools
from modules.definitions.restype_definitions import *
from modules.tools import general_tools
from modules.tools import protocols as protocol_tools
from modules.protocols import for_design
from modules.protocols import loop_minimization
from modules.ScoreBase import *


class FullControlWindow():
    """
    Window for full control of protein backbone and sidechains
    """
    
    def __init__(self, ScoreObject, pose):
        
        #Main classes
        self.pose = pose
        self.design_protocols = for_design.design_protocols(ScoreObject, self.pose)
        self.loop_protocols = loop_minimization.Loop_Min(ScoreObject, self.pose)
        self.residue_definitions = definitions(); #Defines all residue type information.
        self.score_object = ScoreObject
        
        
        
        #GUI variables
        self.resnum = StringVar();
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
        
        
        
    def __exit__(self):
        exit()
        
    def shoInfo(self, res, chain):
        try:
            res = self.pose.pdb_info().pdb2pose(chain, int(res))
            self.restype_var.set(self.pose.residue(res).name()); #Exception comes from this! Segfualts on .psi etc.  This needs to be first. 
        except PyRosettaException:
            print "Residue does not exist in PDB.."
            return
        self.PhiVar.set("%.3f"%self.pose.phi(res)); self.PsiVar.set("%.3f"%self.pose.psi(res)); self.OmeVar.set("%.3f"%self.pose.omega(res))
        self.rotamer_energy.set("%.3f REU"%analysis_tools.return_energy(self.pose, res))
        self.rotamer_probabilty.set("%.3f"%analysis_tools.return_probability(self.pose, res))
        #Calculate Total Energy.
        self.residue_energy_total.set("%.3f REU"%self.score_base.ret_total_weighted_residue_energy(res))
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
    
    def makeWindow(self, main, r=0, c=0):
        
        try :
	    print self.pose.pdb_info().name()
	except AttributeError:
	    tkMessageBox.showwarning(message = 'Please Load a Pose...')
	    return
            
        try:
            self.score_base = ScoreBase(self.pose, self.score_object.score)
            
        except PyRosettaException:
            print "Please Load a pose"
            return
        
        print "Everything is dieing."
            
        self.main = main
        self.main.title("Full Control")
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
        self.button_Delta = Button(self.main, text = "Delta", command=lambda: self.pose.assign(protocol_tools.bbMover(self.pose, int(self.resnum.get()), self.chain.get(), float(self.PhiVar.get()), float(self.PsiVar.get()), float(self.OmeVar.get()))))
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
        self.listbox_mutate_restypes = Listbox(self.main)
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
        self.button_Pack.grid(row=r+15, column=c+1)
        self.listbox_mutate_restypes.grid(row=r+17, column=c+1, rowspan=6)
        self.entry_Etot.grid(row=r+23, column=c+0); self.entry_ResE.grid(row=r+23, column=c+2)
        self.label_Etot.grid(row=r+24, column=c+0); self.label_Terms.grid(row=r+24, column=c+1); self.label_Etype.grid(row=r+24, column=c+2)
        self.listbox_energy_terms.grid(row=r+25, column = 1, rowspan=6)
        self.listbox_mutate_restypes.bind("<Double-Button-1>", lambda event: self.mutRes())
        self.listbox_energy_terms.bind("<Double-Button-1>", lambda event: self.get_residue_energy(self.listbox_energy_terms.get(self.listbox_energy_terms.curselection())))
    
    #### Quick Min ####
        self.button_relax_residue.grid(row=r+17, column=c+2, sticky= W+E)
        self.button_backrub_residue.grid(row=r+17, column=c, sticky= W+E)
        self.button_relax_residue_and_neighbors.grid(row=r+18, column=c+2, sticky= W+E)
        #self.button_relax_residue_and_neighbors_bb.grid(row=r+18, column=c+2, sticky= W+E)
        self.button_backrub_residue_and_neighbors.grid(row=r+18, column=c, sticky= W+E)
        #Need to use a ScoreFxnControl Object.  Would like to use the currently selected Scoretype though - So - Pass the Object!!
        #Populating the Listboxes:
        (zeroscores, nonzeroscores) = self.score_object.scoreOption("Breakdown ScoreFxn")
        for type in nonzeroscores:
            self.listbox_energy_terms.insert(END, type)
            #For now, I will show this.  Later, I will use an object to hold all Residue and Design/CDR info.  Useful in the future.

        for type in self.residue_definitions.restype_info["All"]:
            self.listbox_mutate_restypes.insert(END, type)
    
    def get_residue_energy(self, type_weight_string):
        self.last_type_string_selected = type_weight_string
        type_string = type_weight_string.split(";")[0]
        weight = float(type_weight_string.split(";")[1])
        res = int(self.resnum.get())
        chain = self.chain.get()
        posenum = self.pose.pdb_info().pdb2pose(chain, res)
        emap = self.score_base.ret_residue_energy(posenum)
        e = weight*emap[self.score_base.fa_terms[type_string]]
        self.eterm_type.set(type_string)
        self.residue_energy_of_Eterm.set("%.3f REU"%e)
        
    #### 'CALLBACKS' ####
    def relax_residue(self):
        """
        Relaxes the given residue
        """
        
        self.loop_protocols.RelaxLoop(1, [self.resnum.get()+":"+self.resnum.get()+":"+self.chain.get()], 1)
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
