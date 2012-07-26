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


class FullControl():
    '''
    Window for full control of protein backbone.  Basic for now, maybe will expand in future. Needs a bit of fixing from previous iterations
    '''
    
    def __init__(self, ScoreObject, pose):
        
        #Main classes
        self.pose = pose
        self.design_protocols = for_design.design_protocols(ScoreObject, self.pose)
        self.loop_protocols = loop_minimization.Loop_Min(ScoreObject, self.pose)
        self.residue_definitions = definitions(); #Defines all residue type information.
        self.score_object = ScoreObject
        
        
        
        #GUI variables
        self.ResNumentry_Var = StringVar()
        self.Chainentry_Var = StringVar()
        self.PhiVar = StringVar()
        self.PsiVar = StringVar()
        self.OmeVar = StringVar()
        self.ResVar = StringVar()
        self.CutVar = StringVar()
        self.ResCur = StringVar(); #????
        self.RotProbVar = StringVar()
        self.RotEVar = StringVar()
        self.VacVar = StringVar()
        self.VacVar.set("A")
        self.ETermVar = StringVar()
        self.ETermVar.set("Residue Energy: ")
        self.TotResEVar = StringVar()
        self.ResVarE = StringVar()
        
        
        
    def __exit__(self):
        exit()
        
    def shoInfo(self, res, chain):
        try:
            res = self.pose.pdb_info().pdb2pose(chain, int(res))
            self.ResVar.set(self.pose.residue(res).name()); #Exception comes from this! Segfualts on .psi etc.  This needs to be first. 
        except PyRosettaException:
            print "Residue does not exist in PDB.."
            return
        self.PhiVar.set("%.3f"%self.pose.phi(res)); self.PsiVar.set("%.3f"%self.pose.psi(res)); self.OmeVar.set("%.3f"%self.pose.omega(res))
        self.RotEVar.set("%.3f REU"%analysis_tools.rotamers().retEn(self.pose, res))
        self.RotProbVar.set("%.3f"%analysis_tools.rotamers().retProb(self.pose, res))
        #Calculate Total Energy.
        self.TotResEVar.set("%.3f REU"%self.score_base.ret_total_weighted_residue_energy(res))
        if not self.ETermVar.get()=="Residue Energy: ":
            #get_residue_energy function should split so it accepts the energy type directly.
            self.get_residue_energy(self.last_type_string_selected)
        
    def nextRes(self, res, chain):
        res = int(res)+1
        self.ResNumentry_Var.set(res)
        self.shoInfo(res, chain)
    def prevRes(self, res, chain):
        res = int(res)-1
        self.ResNumentry_Var.set(res)
        self.shoInfo(res, chain)
    '''
    def addRes(self):
        num = self.ResNumentry_Var.get()
        chain = self.Chainentry_Var.get()
        add = num+":"+num+":"+chain
        LisLoop.append(add)
        Lmode.LisLoops.insert(END, add)
    '''
    
    def makeWindow(self, main, r=0, c=0):
        try:
            self.score_base = ScoreBase(self.pose, self.score_object.score)
        except PyRosettaException:
            print "Please Load a pose"
            exit()
            #return
            
        self.main = main
        self.main.title("Full Control")
        self.label_Fu = Label(self.main, text = "Full Control", font="Arial")
        self.button_Tree = Button(self.main, text = "Set Loop Fold Tree", command=lambda: self.pose.assign(tools.loops.initLoops().setLoopBreak(self.pose, VarStartLoop.get(), VarEndLoop.get(), VarLoopChain.get(), self.CutVar.get())))
        self.button_CCD = Button(self.main, text = "Close Loop by CCD", command=lambda: self.pose.assign(Protocols.LoopProtocols().etc().closeLoop(self.pose, 1, VarPym.get(), \
        tools.loops.initLoops().InitializeLoop(self.pose.pdb_info().pdb2pose(VarLoopChain.get(),int(VarStartLoop.get())), self.pose.pdb_info().pdb2pose(VarLoopChain.get(), int(VarEndLoop.get())), self.pose, VarLoopChain.get()), \
        tools.loops.initLoops().InitializeLoopMovemap(self.pose.pdb_info().pdb2pose(VarLoopChain.get(),int(VarStartLoop.get())), self.pose.pdb_info().pdb2pose(VarLoopChain.get(), int(VarEndLoop.get()))))))
        self.entry_Res = Entry(self.main, textvariable=self.ResNumentry_Var, justify=CENTER)
        self.entry_Cha = Entry(self.main, textvariable=self.Chainentry_Var, justify=CENTER)
        self.entry_Nam = Entry(self.main, textvariable=self.ResVar, justify=CENTER)
        self.label_Res = Label(self.main, text = "Residue #")
        self.label_Cha = Label(self.main, text = "Chain ID")
        self.label_Nam = Label(self.main, text = "Residue Name")
        self.button_Get = Button(self.main, text="Get Info", command=lambda: self.shoInfo(self.ResNumentry_Var.get(), self.Chainentry_Var.get()))
        self.button_PRes = Button(self.main, text="Previous Residue", command=lambda: self.prevRes(self.ResNumentry_Var.get(), self.Chainentry_Var.get()))
        self.button_NRes = Button(self.main, text = "Next Residue", command=lambda: self.nextRes(self.ResNumentry_Var.get(), self.Chainentry_Var.get()))
        #self.label_Ba = Label(self.main, text="BackBone", font = "Arial 12 bold")
        self.entry_Cut = Entry(self.main, textvariable=self.CutVar, justify=CENTER)
        self.CutVar.set("CutPoint")
        self.entry_Phi = Entry(self.main, textvariable=self.PhiVar, justify=CENTER)
        self.entry_Psi = Entry(self.main, textvariable=self.PsiVar, justify=CENTER)
        self.entry_Ome = Entry(self.main, textvariable=self.OmeVar, justify=CENTER)
        #self.AddLis = Button(self.main, text = "Add Res Loop List", command = lambda:self.addRes())
        self.label_Phi = Label(self.main, text="Phi")
        self.label_Psi = Label(self.main, text="Psi")
        self.label_Ome = Label(self.main, text="Omega")
        self.button_Delta = Button(self.main, text = "Delta", command=lambda: self.pose.assign(protocol_tools.bbMover(self.pose, int(self.ResNumentry_Var.get()), self.Chainentry_Var.get(), float(self.PhiVar.get()), float(self.PsiVar.get()), float(self.OmeVar.get()))))
        #self.button_previous_pose = Button(self.main, text = "Previous Pose", command = lambda: loaddownFromPList(len(pList)-pListTrack))
        #self.button_next_pose = Button(self.main, text = "Next Pose", command = lambda: loadupFromPList(len(pList)+pListTrack))
        #self.label_Rot = Label(self.main, text = "Rotamers", font = "Arial 12 bold")
        self.label_RotPr = Label(self.main, text = "Aprox Rotamer Probability")
        self.label_RotEn = Label(self.main, text = "Rotamer Energy")
        self.entry_rotamer_probability = Label(self.main, textvariable=self.RotProbVar, justify=CENTER, relief=SUNKEN)
        self.entry_rotamer_energy = Label(self.main, textvariable=self.RotEVar, justify=CENTER, relief=SUNKEN)
        self.SliPhi = Scale(self.main, from_=-180, to=180, orient=HORIZONTAL,variable = self.PhiVar)
        #command=lambda:p.assign(general_tools.protocols().bbMover(p, int(VarPym.get()), int(self.ResNumentry_Var.get()), self.Chainentry_Var.get(), float(self.PhiVar.get()), float(self.PsiVar.get()), float(self.OmeVar.get())))
        self.SliPsi = Scale(self.main, from_=-180, to=180, orient=HORIZONTAL, variable = self.PsiVar)
        #command=lambda:p.assign(general_tools.protocols().bbMover(p, int(VarPym.get()), int(self.ResNumentry_Var.get()), self.Chainentry_Var.get(), float(self.PhiVar.get()), float(self.PsiVar.get()), float(self.OmeVar.get())))
        '''
        self.label_change_rotamer = Label(self.main, text = "Change Rotamer")
        self.button_up_rotamer = Button(self.main, text = "Up Rotamer")
        self.button_down_rotamer = Button(self.main, text = "Down Rotamer")
        self.button_change_library = Button(self.main, text = "Change Library")
        self.button_update_library = Button(self.main, text = "Update")
        '''
        
        #self.label_Resi = Label(self.main, text = "Residues", font="Arial 12 bold")
        #self.label_Prot = Label(self.main, text = "Protocols", font="Arial 12 bold")
        #This is all of the mutagenesis and repacking stuff....
    #### Quick Protocols ####
        self.button_relax_residue = Button(self.main, text = "Relax Residue", command=lambda: self.relRes())
        self.button_backrub_residue= Button(self.main, text = "Backrub Residue", state=DISABLED)
        self.button_relax_residue_and_neighbors=Button(self.main, text= "Relax X-Res-X", command = lambda: self.relRes_Neighbors(False))
        self.button_relax_residue_and_neighbors_bb=Button(self.main, text= "Relax X-Res-X (BB only)", command = lambda: self.relRes_Neighbors(True))
        self.button_backrub_residue_and_neighbors=Button(self.main, text = "Backrub X-Res-X", state=DISABLED)
        self.button_Pack = Button(self.main, text = "Pack Rotamer", command=lambda: self.packRes())
        #self.button_Undu = Button(self.main, text= "Undo Mutagenesis", command = lambda: loadFromPList(len(pList)-1))
        self.LisMut = Listbox(self.main)
        #self.check_button_Pym = Checkbutton(self.main, text = "Show in Pymol", variable=VarPym)
        self.check_button_Vac = Checkbutton(self.main, text = "Pack in Vaccinity")
        self.entry_Vac = Label(self.main, textvariable = self.VacVar, justify=CENTER, relief=SUNKEN)
        self.label_Etot = Label(self.main, text = "Total Residue Energy")
        self.entry_Etot = Label(self.main, textvariable = self.TotResEVar, relief=SUNKEN)
        self.label_Terms = Label(self.main, text = "energy terms")
        self.entry_ResE= Label(self.main, textvariable = self.ResVarE, relief=SUNKEN)
        self.label_Eres = Label(self.main, textvariable = self.ETermVar)
        self.LisTerms = Listbox(self.main)
        
        
    #### GRID ####
        self.entry_Res.grid(row=r+5, column=c+0); self.entry_Cha.grid(row=r+5, column=c+1); self.entry_Nam.grid(row=r+5, column=c+2)
        self.label_Res.grid(row=r+6, column=c+0); self.label_Cha.grid(row=r+6, column=c+1); self.label_Nam.grid(row=r+6,column=c+2)
        self.button_PRes.grid(row=r+7, column=c+0); self.button_Get.grid(row=r+7, column=c+1, sticky=W+E); self.button_NRes.grid(row=r+7, column=c+2)
        self.entry_Phi.grid(row=r+8, column=c+0); self.entry_Psi.grid(row=r+8, column=c+1); self.entry_Ome.grid(row=r+8, column=c+2)
        self.SliPhi.grid(row=r+9, column=c+0); self.SliPsi.grid(row=r+9, column=c+1)
        self.label_Phi.grid(row=r+10, column=c+0); self.label_Psi.grid(row=r+10, column=c+1); self.label_Ome.grid(row=r+10, column=c+2)
        self.button_Delta.grid(row=r+11, column=c+1); 
        self.entry_rotamer_probability.grid(row=r+12, column=c+0); self.entry_rotamer_energy.grid(row=r+12, column=c+2)
        self.label_RotPr.grid(row=r+13, column=c+0); self.label_RotEn.grid(row=r+13, column=c+2)
        
    #### Mutagenesis, Repacking, Residue Energies ####
        self.button_Pack.grid(row=r+15, column=c+1)
        #self.check_button_Vac.grid(row=r+16, column=c+0); self.entry_Vac.grid(row=r+16, column=c+1); #self.check_button_Pym.grid(row=r+16, column=c+2)
        self.LisMut.grid(row=r+17, column=c+1, rowspan=6)
        self.entry_Etot.grid(row=r+23, column=c+0); self.entry_ResE.grid(row=r+23, column=c+2)
        self.label_Etot.grid(row=r+24, column=c+0); self.label_Terms.grid(row=r+24, column=c+1); self.label_Eres.grid(row=r+24, column=c+2)
        self.LisTerms.grid(row=r+25, column = 1, rowspan=6)
        self.LisMut.bind("<Double-Button-1>", lambda event: self.mutRes())
        self.LisTerms.bind("<Double-Button-1>", lambda event: self.get_residue_energy(self.LisTerms.get(self.LisTerms.curselection())))
    
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
            self.LisTerms.insert(END, type)
            #For now, I will show this.  Later, I will use an object to hold all Residue and Design/CDR info.  Useful in the future.

        for type in self.residue_definitions.restype_info["All"]:
            self.LisMut.insert(END, type)
    
    def get_residue_energy(self, type_weight_string):
        self.last_type_string_selected = type_weight_string
        type_string = type_weight_string.split(";")[0]
        weight = float(type_weight_string.split(";")[1])
        res = int(self.ResNumentry_Var.get())
        chain = self.Chainentry_Var.get()
        posenum = self.pose.pdb_info().pdb2pose(chain, res)
        emap = self.score_base.ret_residue_energy(posenum)
        e = weight*emap[self.score_base.fa_terms[type_string]]
        self.ETermVar.set(type_string)
        self.ResVarE.set("%.3f REU"%e)
        
    #### 'CALLBACKS' ####
    def relRes(self):
        self.loop_protocols.RelaxLoop(1, [self.ResNumentry_Var.get()+":"+self.ResNumentry_Var.get()+":"+self.Chainentry_Var.get()], 1)
        self.shoInfo(self.ResNumentry_Var.get(), self.Chainentry_Var.get())
    def relRes_Neighbors(self, bbonly=False):
        if bbonly:
            self.loop_protocols.relax_residue_and_neighbors(1, (self.ResNumentry_Var.get()), self.Chainentry_Var.get(), 1, True)
        else:
            self.loop_protocols.relax_residue_and_neighbors(1, (self.ResNumentry_Var.get()), self.Chainentry_Var.get())
        self.shoInfo(self.ResNumentry_Var.get(), self.Chainentry_Var.get())
    def packRes(self):
        '''
        Lovely Tkinter..
        '''
        self.design_protocols.pack_residue(self.ResNumentry_Var.get(), self.Chainentry_Var.get())
        self.shoInfo(self.ResNumentry_Var.get(), self.Chainentry_Var.get())
    def mutRes(self):
        self.design_protocols.mutateRes(self.ResNumentry_Var.get(), self.Chainentry_Var.get(), self.LisMut.get(self.LisMut.curselection()))
        self.shoInfo(self.ResNumentry_Var.get(), self.Chainentry_Var.get())
