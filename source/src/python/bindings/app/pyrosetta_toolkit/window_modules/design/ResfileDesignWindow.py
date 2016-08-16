#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/design_window/resfile_design.py
## @brief  Main resfile design window.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os

#Tkinter Imports
from Tkinter import *
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.modules.definitions.restype_definitions import *
from app.pyrosetta_toolkit.window_main import global_variables

class ResfileDesignWindow:
    def __init__(self, main, DesignDic, pose):
        self.main=main
        self.main.title("Design Toolbox")

        self.pwd = self.location()[0]
        self.DesignDic = DesignDic; #Main dictionary for saving the resfile info. Can be empty upon construction of window.
        #current and design are labeled opposite . need to fix.
        self.current_accessible_sa = StringVar(); self.current_relative_mutability = StringVar(); self.current_surface_probability = StringVar()
        self.design_accessible_sa = StringVar(); self.design_relative_mutability = StringVar(); self.design_surface_probability = StringVar()
        self.current_residue_name = StringVar(); #Split three letter aa code.
        self.current_residue_name_full = StringVar(); #Full rosetta code including variants.

        self.current_residue = StringVar();
        self.current_chain = StringVar();
        self.pose = pose
        self.TypeCurSelection = StringVar(); #Current selection of restype
        self.residue_definitions = definitions(); #Defines all residue type information.


    def setTk(self):
        #self.entry_Res = Entry(self.main, textvariable=addRes, justify=CENTER)
        #self.entry_Chain = Entry(self.main, textvariable=DesignChain, justify=CENTER)

        self.label_Chain = Label(self.main, text = "Chain")
        self.shoRes = Label(self.main, textvariable=self.current_residue_name_full)
        self.label_Res = Label(self.main, text = "Residue or Start:End")
        self.listbox_restypes = Listbox(self.main)
        self.listbox_residues = Listbox(self.main)
        self.entry_residue = Entry(self.main, textvariable=self.current_residue, justify=CENTER)
        self.entry_chain = Entry(self.main, textvariable=self.current_chain, justify=CENTER)
        self.button_check_residue = Button(self.main, text="Check Residue", command=lambda: self.check_residue_callback())
        self.button_next_residue = Button(self.main, text="Next Res", command=lambda: self.nextRes())
        self.button_previous_residue= Button(self.main, text="Previous Res", command=lambda: self.prevRes())
        self.listbox_current_designs = Listbox(self.main)
        self.button_save_resfile = Button(self.main, text="Write Resfile", command=lambda: self.saveDesign())
        self.button_clear_resfile = Button(self.main, text="Clear Data", command = lambda: self.clearResData())
        self.button_ClearRes = Button(self.main, text = "Clear Residue", command=lambda: self.clearCurrentResidue())
        self.button_load_resfile = Button(self.main, text="Load Resfile")

        #### Mutability Data ####
        self.label_accessible_SA = Label(self.main, text="Accesible SA");
        self.label_relative_mutability = Label(self.main, text="Rel. Mutability");
        self.label_surface_probability = Label(self.main, text = "Surface Probability");

        self.active_label_current_accessible_sa = Label(self.main, textvariable=self.current_accessible_sa)
        self.active_label_current_relative_mutability = Label(self.main, textvariable=self.current_relative_mutability)
        self.active_label_current_surface_probability=Label(self.main, textvariable=self.current_surface_probability)

        self.active_label_design_accessible_sa = Label(self.main, textvariable=self.design_accessible_sa)
        self.active_label_design_relative_mutability = Label(self.main, textvariable=self.design_relative_mutability)
        self.active_label_design_surface_probability = Label(self.main, textvariable=self.design_surface_probability)

        #### Photo ####
        InterfacePhoto =PhotoImage(file = (self.pwd+"/Media/InterfaceDesignSmall.gif"))
        self.Photo = Label(self.main, image=InterfacePhoto)
        self.Photo.image = InterfacePhoto

    def shoTk(self):
        #self.entry_Res.grid(row=1, column=1)
        self.label_Chain.grid(row=2, column=0)
        self.listbox_restypes.grid(row=0, column=4, rowspan=10)
        self.listbox_residues.grid(row=0, column=5, rowspan=10)
        #self.entry_Chain.grid(row=3, column=1)
        self.label_Res.grid(row=1, column=0)
        self.shoRes.grid(row=0, column=1)
        self.entry_residue.grid(row=1, column=1)
        self.button_check_residue.grid(row=3, column=1)
        self.button_previous_residue.grid(row=3, column=0)
        self.button_next_residue.grid(row=3, column=2)
        self.listbox_current_designs.grid(row=12, column=4, rowspan=10)
        self.entry_chain.grid(row=2, column=1)
        self.button_save_resfile.grid(row=12, column=5, sticky=W+E)
        #self.button_load_resfile.grid(row=13, column=5, sticky=W+E); #Load resfile does not work yet.
        self.button_clear_resfile.grid(row=14, column=5, sticky=W+E)
        #self.button_ClearRes.grid(row=15, column=5, sticky=W+E)
        self.label_accessible_SA.grid(row=4, column=0); self.active_label_current_accessible_sa.grid(row=6, column=0); self.active_label_design_accessible_sa.grid(row=5, column=0)
        self.label_relative_mutability.grid(row=4, column=1); self.active_label_current_relative_mutability.grid(row=6, column=1); self.active_label_design_relative_mutability.grid(row=5, column=1)
        self.label_surface_probability.grid(row=4, column=2); self.active_label_current_surface_probability.grid(row=6, column=2); self.active_label_design_surface_probability.grid(row=5, column=2)
        #ListBox Behavior
        self.listbox_restypes.bind("<ButtonRelease-1>", lambda event: self.click_restype_callback(self.listbox_restypes, self.listbox_residues))
        self.listbox_residues.bind("<Double-Button-1>", lambda event: self.add_to_current_designs_callback(self.listbox_residues, self.listbox_current_designs))
        self.listbox_residues.bind("<ButtonRelease-1>", lambda event: self.updateInfo(self.listbox_residues))
        self.listbox_current_designs.bind("<Double-Button-1>", lambda event: self.remove_from_current_designs_callback(self.listbox_current_designs))
        #Photo
        self.Photo.grid(row =12, column = 6, rowspan=5)
        self.main.grid_columnconfigure(ALL, weight=1)
#### GUI FUNCTIONS ####

    #### 'CALLBACKS' ####

    def check_residue_callback(self):
        """
        What happens when you click check residue.  Should be a callback, but need both res and chain specified before update.
        """

        if not self.check_for_residue_existance():
            return

        self.listbox_current_designs.delete(0, END)
        self.current_residueFull = self.current_residue.get().split(":")

        if len(self.current_residueFull) > 1:
            self.listbox_restypes.delete(9)
            self.listbox_restypes.insert(9, "Conserved")
        else:
            res = self.current_residueFull[0] +":"+self.current_chain.get()
            #print res
            if not self.pose.total_residue():
                return

            resType = self.pose.pdb_info().pdb2pose(self.current_chain.get(), int(self.current_residue.get()))
            try:
                resName = self.pose.residue(resType).name()
            except PyRosettaException:
                print "Residue does not exist in PDB..."
                return


            self.current_residue_name_full.set(resName)
            three_letter_name = resName.split("_")[0]; #Fix for Rosetta designated chain endings, His_d
            self.current_residue_name.set(resName)

            #Fix for not having data for certain residue types:
            try:
                self.design_accessible_sa.set(self.residue_definitions.resinfo[three_letter_name][0])
                self.design_relative_mutability.set(self.residue_definitions.resinfo[three_letter_name][1])
                self.design_surface_probability.set(self.residue_definitions.resinfo[three_letter_name][2])
            except KeyError:
                print "No Data for restype"
                self.design_accessible_sa.set("")
                self.design_relative_mutability.set("")
                self.design_surface_probability.set("")

            con = "Conserved:" +three_letter_name
            self.listbox_restypes.delete(9)
            self.listbox_restypes.insert(9, con)
            if self.DesignDic.has_key(res):
                for residues in self.DesignDic[res]:
                    self.listbox_current_designs.insert(END, residues)

    def click_restype_callback(self, ListTypes, ListTypesFull):
        """
        What happens when you click a design type (conserved, polar, etc.)
        """

        #Ropulate if more then 1 residue
        #if self.current_residue.get().split(":")>1:
            #self.check_residue_callback()

        ListTypesFull.delete(0, END)
        type = ListTypes.get(ListTypes.curselection())
        self.TypeCurSelection.set(type)
        ListTypesFull.insert(END, "All")
        ListTypesFull.insert(END, "All+Self")
        try:
            for res in self.residue_definitions.restype_info[type]:
                ListTypesFull.insert(END, res)
        except KeyError:
            return

    def add_to_current_designs_callback(self, ListTypesFull, check_button_ckList):
        """
        What happens when you add a residue type to the design.  Updates other listboxes, etc.
        """

        if not self.check_for_residue_existance():
            return

        self.current_chainFull = self.current_residue.get().split(":")
        if len(self.current_chainFull) > 1:
            ResStart = int(self.current_chainFull[0]); ResEnd = int(self.current_chainFull[1])
            if (ListTypesFull.get(ListTypesFull.curselection()))=="All":
                for ResType in self.residue_definitions.restype_info[self.TypeCurSelection.get()]:
                    print ResType
                    for i in range(ResStart, ResEnd+1):
                        res = repr(i) + ":"+self.current_chain.get()
                        if self.DesignDic.has_key(res):
                            self.DesignDic[res].append(ResType)
                        else:
                            self.DesignDic[res] = [ResType, ]
            if (ListTypesFull.get(ListTypesFull.curselection()))=="All+Self":
                start = self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResStart)
                end =   self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResEnd)
                for ResType in self.residue_definitions.restype_info[self.TypeCurSelection.get()]:
                    for i in range(ResStart, ResEnd+1):
                        res = repr(i) + ":"+self.current_chain.get()
                        if self.DesignDic.has_key(res):
                            self.DesignDic[res].append(ResType)
                        else:
                            self.DesignDic[res] = [ResType, ]
                        #Gets number according to rosetta
                        resType = self.pose.pdb_info().pdb2pose(self.current_chain.get(), i)
                        #Gets name according to rosetta
                        resType = self.pose.residue(resType).name()
                        resType = resType.split("_")[0];
                        for y in self.residue_definitions.restype_info["All"]:
                            z = y.split(":")
                            if resType == z[1]:
                                self.DesignDic[res].append(y)
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="NATRO":
                for i in range(ResStart, ResEnd+1):
                        res = repr(i) + ":"+self.current_chain.get()
                        self.DesignDic[res] = ["NATRO", ]
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="NATAA":
                for i in range(ResStart, ResEnd+1):
                        res = repr(i) + ":"+self.current_chain.get()
                        self.DesignDic[res] = ["NATAA", ]

            elif (ListTypesFull.get(ListTypesFull.curselection()))=="ALLAA":
                for i in range(ResStart, ResEnd+1):
                        res = repr(i) + ":"+self.current_chain.get()
                        self.DesignDic[res] = ["ALLAA", ]
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="All Conserved Mutations":
                start = self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResStart)
                end =   self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResEnd)
                print repr(start); print repr(end)
                for i in range(start, end+1):
                    print i
                    resType = i
                    resType = self.pose.residue(resType).name()
                    resType = resType.split("_")[0];
                    resType = "Conserved:"+resType
                    res = self.pose.pdb_info().pose2pdb(i).split()[0]+":"+self.current_chain.get()
                    for types in self.residue_definitions.restype_info[resType]:
                        if self.DesignDic.has_key(res):
                            self.DesignDic[res].append(types)
                        else:
                            self.DesignDic[res] = [types, ]
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="All Conserved Mutations+Self":
                print "4"
                start = self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResStart)
                end =   self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResEnd)
                for i in range(start, end+1):
                    resType = i
                    resType = self.pose.residue(resType).name()
                    resType = resType.split("_")[0];
                    res = self.pose.pdb_info().pose2pdb(i).split()[0]+":"+self.current_chain.get()
                    resType = "Conserved:"+resType
                    for types in self.residue_definitions.restype_info[resType]:
                        if self.DesignDic.has_key(res):
                            self.DesignDic[res].append(types)
                        else:
                            self.DesignDic[res] = [types, ]

                    resTypei = self.pose.residue(i).name()
                    for y in self.residue_definitions.restype_info["All"]:
                        z = y.split(":")
                        if resTypei == z[1]:
                            self.DesignDic[res].append(y)
            else:
                ResType= ListTypesFull.get(ListTypesFull.curselection())
                for i in range(ResStart, ResEnd+1):
                    res = repr(i) + ":"+ self.current_chain.get()
                    if self.DesignDic.has_key(res):
                        self.DesignDic[res].append(ResType)
                    else:
                        self.DesignDic[res] = [ResType, ]
                    #self.updateInfo(ListTypesFull)
            check_button_ckList.delete(0, END)
            res = repr(ResStart) + ":" + self.current_chain.get()
            for residues in self.DesignDic[res]:
                check_button_ckList.insert(END, residues)
        else:
            resType = self.pose.pdb_info().pdb2pose(self.current_chain.get(), int(self.current_residue.get()))
            res = self.pose.residue(resType).name()

            self.current_residue_name_full.set(res)
            res= res.split("_")[0];
            self.current_residue_name.set(res)
            try:
                self.design_accessible_sa.set(self.residue_definitions.resinfo[self.pose.residue(resType).name()][0])
                self.design_relative_mutability.set(self.residue_definitions.resinfo[self.pose.residue(resType).name()][1])
                self.design_surface_probability.set(self.residue_definitions.resinfo[self.pose.residue(resType).name()][2])
            except KeyError:
                pass
            res = self.current_residue.get() +":" + self.current_chain.get()
            if (ListTypesFull.get(ListTypesFull.curselection()))=="All":
                for ResType in self.residue_definitions.restype_info[self.TypeCurSelection.get()]:
                    if self.DesignDic.has_key(res):
                        self.DesignDic[res].append(ResType)
                    else:
                        self.DesignDic[res] = [ResType, ]
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="All+Self":
                for ResType in self.residue_definitions.restype_info[self.TypeCurSelection.get()]:
                    if self.DesignDic.has_key(res):
                        self.DesignDic[res].append(ResType)
                    else:
                        self.DesignDic[res] = [ResType, ]
                for y in self.residue_definitions.restype_info["All"]:
                    z = y.split(":")
                    if self.current_residue_name.get() == z[1]:
                        self.DesignDic[res].append(y)

            elif (ListTypesFull.get(ListTypesFull.curselection()))=="NATRO":
                self.DesignDic[res] = ["NATRO", ]
            elif (ListTypesFull.get(ListTypesFull.curselection()))=="NATAA":
                self.DesignDic[res] = ["NATAA", ]

            elif (ListTypesFull.get(ListTypesFull.curselection()))=="ALLAA":
                self.DesignDic[res] = ["ALLAA", ]
            else:
                ResType= ListTypesFull.get(ListTypesFull.curselection())
                if self.DesignDic.has_key(res):
                    self.DesignDic[res].append(ResType)
                else:
                    self.DesignDic[res] = [ResType, ]
                self.updateInfo(ListTypesFull)
            check_button_ckList.delete(0, END)
            for residues in self.DesignDic[res]:
                check_button_ckList.insert(END, residues)

    def remove_from_current_designs_callback(self, check_button_ckList):
        """
        What happens when you double click the bottom listbox to remove a residue.
        """

        self.current_residueFull = self.current_residue.get().split(":")
        if len(self.current_residueFull) > 1:
            ResStart = int(self.current_residueFull[0]); ResEnd = int(self.current_residueFull[1])
            for i in range(ResStart, ResEnd +1):
                res = repr(i) + ":" +self.current_chain.get()
                self.DesignDic[res].remove(check_button_ckList.get(check_button_ckList.curselection()))
            check_button_ckList.delete(0, END)
            res = repr(ResStart) + ":" +self.current_chain.get()
            for residues in self.DesignDic[res]:
                check_button_ckList.insert(END, residues)
        else:
            res = self.current_residue.get() +":" + self.current_chain.get()
            self.DesignDic[res].remove(check_button_ckList.get(check_button_ckList.curselection()))
            check_button_ckList.delete(0, END)
            for residues in self.DesignDic[res]:
                check_button_ckList.insert(END, residues)

    def check_for_residue_existance(self):
        """
        Checks to make sure residue/residues exist before adding to the residue checkbox.
        Returns False if something is wrong.
        """
        if not self.pose.total_residue():
            print "No pose Loaded."
            return False

        if not self.current_chain.get() or not self.current_residue.get():
            print "Chain or residue not set"
            return False

        current_region = self.current_residue.get().split(":")
        if len(current_region)>1:
            ResStart = int(current_region[0]); ResEnd = int(self.current_region[1])


            if self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResStart)==0 or self.pose.pdb_info().pdb2pose(self.current_chain.get(), ResEnd)==0:
                print "Region not found in pose"
                return False
        else:
            if self.pose.pdb_info().pdb2pose(self.current_chain.get(), int(self.current_residue.get())) ==0:

                print "Residue not found in pose"
                return False

        #If everythig is good then:
        return True

    def updateInfo(self, ListTypesFull):
        Restype = ListTypesFull.get(ListTypesFull.curselection())
        typefull = Restype.split(":")
        self.current_residueFull = self.current_residue.get().split(":")

        if (len(self.current_residueFull) < 2) and (len(typefull)==3):
            self.current_accessible_sa.set(self.residue_definitions.resinfo[typefull[1]][0])
            self.current_relative_mutability.set(self.residue_definitions.resinfo[typefull[1]][1])
            self.current_surface_probability.set(self.residue_definitions.resinfo[typefull[1]][2])

    #### Other GUI Functions ####

    def nextRes(self):
        """
        Update functions when residue is +1
        """

        self.current_residueFull = self.current_residue.get().split(":")
        if len(self.current_residueFull) < 2:
            res = int(self.current_residue.get())+1
            self.current_residue.set(res)
            self.check_residue_callback()
        else:
            print "Please Change selection to get info..."

    def prevRes(self):
        """
        Update functions when residue is -1
        """

        self.current_residueFull = self.current_residue.get().split(":")
        if len(self.current_residueFull) < 2:
            res = int(self.current_residue.get())-1
            self.current_residue.set(res)
            self.check_residue_callback()
        else:
            print "Please Change selection to get info..."

    def setTypes(self):
        for type in self.residue_definitions.restypes:
            self.listbox_restypes.insert(END, type)

    def clearResData(self):
        self.DesignDic.clear()
        self.check_residue_callback()
        print "Design Data Cleared..."
        print self.DesignDic

    def clearCurrentResidue(self):
        pass


#### FUNCTIONS ####
    def saveDesign(self):
        out = tkFileDialog.asksaveasfilename(initialdir = global_variables.current_directory, title ="Save As...")
        if not out: return
        output_tools.save_resfile_w_designdic(self.pose, self.DesignDic, out)
        outSP = out.split("/")
        length = len(outSP)
        folder = outSP[length-2]
        if folder == "FragSets_Designs":
            prot.LisProt1Frag.insert(END, outSP[length-1])

    def location(self):
      """
      Allows the script to be self-aware of it's path.
      So that it can be imported/ran from anywhere.
      """

      p = os.path.abspath(__file__)
      pathSP = os.path.split(p)
      return pathSP
