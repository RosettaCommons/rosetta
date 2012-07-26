#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/analysis/AnalyzePyRosetta.py
## @brief  Was the main analysis GUI.  Needs a ton of work. 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
import glob
import tkFileDialog
import tkMessageBox
import tkSimpleDialog


#Needs to be completely rewritten!
class AnalyzePyRosetta():
    def __init__(self, main):
        self.main=main
        self.main.title("Decoy Analysis Window")
        self.outAlign = StringVar()
        self.outAlign.set(0)
        self.SeqOut=StringVar()
        self.VarSeqOut = StringVar()
        self.seqOut = StringVar()
        self.seqOut.set(0)
        self.VarAnalFile = StringVar()
        self.VarAnalTop = StringVar()
        self.AnalysisSettings = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
        self.AnalysisSettings[0].set(0); AnalysisSettings[1].set(0); AnalysisSettings[2].set(1); AnalysisSettings[3].set(0)
        self.AnalysisSettings[5].set("total_score")
        self.Etypes = ("total_score", "rmsd", "fa_atr", "fa_rep", "fa_sol", "fa_intra_rep", "pro_close", "fa_pair", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc", "dslf_ss_dst", "dslf_cs_ang", "dslf_ca_dih", "rama", "omega", "fa_dun", "p_aa_pp")
        self.EtypesRes = ("fa_atr", "fa_rep", "fa_sol", "fa_intra_rep", "pro_close", "fa_pair", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc", "dslf_ss_dst", "dslf_cs_ang", "dslf_ca_dih", "fa_dun", "p_aa_pp")

    #def kickAnalyze(self, list, var, pos):
        #list[pos]=var
    
    def analyzeSeq(self):
        if filename1.get()== "0":
            tkMessageBox.showerror(message="Please choose PDB file to compare results to....")

            return
        else:
            pass
        
        
        if not LisLoop:
            tkMessageBox.showerror(title="Loop", message="Please Specify the regions you would like to compare...")
            ans = tkMessageBox.askquestion(message = "Regions not Specified - Do you want to compare the whole protein?", default = tkMessageBox.NO)
            if ans == "no":
                return

        if self.outAlign.get()=="1":
            print "Saving File"
            if not self.VarSeqOut.get():
                tkMessageBox.showerror(message="Please Specify a filename")
                return
            elif not dirname1.get():
                tkMessageBox.showerror(message="Please speficy a directory with PDB's")
                return
            else:
                if not dirnameout.get()==0:
                    print "Saving to Decoy's Folder..."
                    tools.analysis.Decoys().AnalyzeSeq(LisLoop, AnalysisSettings, VarAnalFile.get(), dirname1.get(), dirname1.get()+"/"+self.VarSeqOut.get(), filename1.get(), self.seqOut.get())
                else:
                    print "Saving to Given Directory"
                    tools.analysis.Decoys().AnalyzeSeq(LisLoop, AnalysisSettings, VarAnalFile.get(), dirname1.get(), dirnameout.get()+"/"+self.VarSeqOut.get(), filename1.get(), self.seqOut.get())
        elif self.outAlign.get()=="0":
            fileout = 0
            tools.analysis.Decoys().AnalyzeSeq(LisLoop, AnalysisSettings, VarAnalFile.get(), dirname1.get(), fileout, filename1.get(), self.seqOut.get())
        
    def setType(self):
        type = self.ListAnal.get(self.ListAnal.curselection())
        AnalysisSettings[5].set(type)
        self.SeqOut.set(type)
    def setAndUpdate(self):
        self.setType()
        #tools.analysis.Decoys().findLowDecoys(dirname1.get(), VarAnalFile.get(), AnalysisSettings, LisLoop, dirnameout.get())
        self.LisAnalCombo.insert(END, self.ListAnal.get(self.ListAnal.curselection()))
    def remLis(self):
        self.LisAnalCombo.delete(self.LisAnalCombo.index(self.LisAnalCombo.curselection()))
    def setTk(self):

        
        self.label_Analy=Label(self.main, text="Analysis Menu")
        self.label_DIR = Label(self.main, text = "Choose Directory")
        self.button_DIR = Button(self.main, text= "Open", command=lambda: dirname1.set(tools.input.tk_get_directory(dirname1.get())))
        self.label_File = Label(self.main, text = "File Name Contains:")
        self.entry_File = Entry(self.main, textvariable=VarAnalFile, justify=CENTER)
        self.label_Top = Label(self.main, text ="Return Top __%:")
        self.label_Top2 =Label(self.main, text = "Sort by Score Type:")
        self.entry_Top = Entry(self.main, textvariable = AnalysisSettings[4], justify=CENTER)
        self.entry_Top2 = Entry(self.main, textvariable = AnalysisSettings[5], justify=CENTER)
        self.check_button_RMSD = Checkbutton(self.main, text= "Find Low RMSD", variable=AnalysisSettings[0])
        self.check_button_ENER = Checkbutton(self.main, text = "Find Low Energies", variable=AnalysisSettings[2])
        self.check_button_Lrmsd = Checkbutton(self.main, text = "Find Low Loop RMSD", variable=AnalysisSettings[1])
        self.label_Move = Label(self.main, text = "Move Top Scoring Poses to output DIR?")
        self.check_button_Move = Checkbutton(self.main, text = "Yes", variable = AnalysisSettings[3])
        self.button_Analy = Button(self.main, text = "Analyze", command=lambda: tools.analysis.Decoys().findLowDecoys(dirname1.get(), VarAnalFile.get(), AnalysisSettings, LisLoop, dirname1.get()))
        self.check_button_Seq = Checkbutton(self.main, text = "Output Alignment?", variable = self.outAlign)
        self.check_button_File = Checkbutton(self.main, text= "Output Full Seq Data?", variable = self.seqOut)
        self.entry_Seq = Entry(self.main, justify = CENTER, textvariable = self.VarSeqOut)
        
        self.button_Seq = Button(self.main, text = "Align and Compare Sequences", command = lambda: self.analyzeSeq())
        self.ListAnal = Listbox(self.main)
        self.LisAnalCombo = Listbox(self.main)
        self.label_Combo = Label(self.main, text = "Decoy Filter Selections")
        
    def shoTk(self):
        #self.label_Analy.grid(row = 0, column=1)
        self.label_DIR.grid(row=1, column=0); self.button_DIR.grid(row=1, column=1)
        self.label_File.grid(row=2, column=0); self.entry_File.grid(row=2, column=1)
        self.label_Top.grid(row=3, column=0); self.entry_Top.grid(row=3, column=1)
        self.label_Top2.grid(row=4, column=0); self.entry_Top2.grid(row=4, column=1)
        self.check_button_RMSD.grid(row=5, column=0); self.check_button_ENER.grid(row=5, column=1)
        self.check_button_Lrmsd.grid(row=6, column=0);
        self.label_Move.grid(row=7, column=0); self.check_button_Move.grid(row=7, column=1)
        self.button_Analy.grid(row=8, column=1)
        self.ListAnal.grid(row=2, column=2, rowspan=8)
        self.label_Combo.grid(row=1, column=2, columnspan=2)
        self.LisAnalCombo.grid(row=2, column=3, rowspan=8)
        self.button_Seq.grid(row =10, column = 0); self.check_button_Seq.grid(row=10, column=1); self.check_button_File.grid(row=10, column=2); self.entry_Seq.grid(row=10, column=3)
        self.ListAnal.bind("<ButtonRelease-1>", lambda event: self.setType())
        self.ListAnal.bind("<Double-Button-1>", lambda event: self.setAndUpdate())
        self.LisAnalCombo.bind("<Double-Button-1>", lambda event: self.remLis())
    def setLis(self):
        for type in Etypes:
            self.ListAnal.insert(END, type)