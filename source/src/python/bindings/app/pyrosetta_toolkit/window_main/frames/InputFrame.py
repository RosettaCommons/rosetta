#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_main/input.py
## @brief  Handles the Regional input of the main GUI window
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog

#Toolkit Imports
from app.pyrosetta_toolkit.modules.Region import Region
from app.pyrosetta_toolkit.modules.tools import sequence as sequence_tools
from app.pyrosetta_toolkit.modules.tools import input as input_tools
from app.pyrosetta_toolkit.modules.tools import general_tools as gen_tools

#from app.pyrosetta_toolkit. import main_window
from app.pyrosetta_toolkit.window_main.IO.GUIInput import GUIInput

class InputFrame(Frame):

    def __init__(self, main, toolkit, input_class, **options):
        Frame.__init__(self, main, **options)
        self.toolkit = toolkit
        self.pose = self.toolkit.pose

        self.input_class = input_class
        self.pwd= self.location()[0]

        self.create_GUI_objects()
        self.grid_GUI_objects()

        #Ignore this.  It is for Komodo Autocomplete.
        if 0:
            self.main = Tk()
            #self.toolkit = main_window()
            self.input_class = GUIInput()

    def create_GUI_objects(self):

        ############## LOOPS ##############
        #self.label_Loop = Label(self, text="Region Selection", font=("Arial"))
        self.AddLoopButton = Button(self, text = "Add Region", command=lambda: self.add_region())
        self.RmLoopButton = Button(self, text = "Remove Region", command = lambda:self.remove_region())
        self.StartLoopLabel=Label(self, text="Start of Region:")
        self.EndLoopLabel=Label(self, text="End of Region:")
        self.ChainIDLabel=Label(self, text="Chain ID:")
        self.loops_listbox = Listbox(self)
        self.loops_scroll = Scrollbar(self)
        self.loops_listbox.config(yscrollcommand = self.loops_scroll.set); self.loops_scroll.config(command = self.loops_listbox.yview)
        self.StartLoopEntry=Entry(self, textvariable=self.input_class.region_start)
        self.EndLoopEntry=Entry(self, textvariable=self.input_class.region_end)
        self.ChainIDEntry=Entry(self, textvariable=self.input_class.region_chain)
        ##################################


        ############ SEQUENCE ############

        self.ShoSeqButton = Button(self, text="Show Sequence", command=lambda: self.show_sequence())

        ##################################


           ### Photo ###

        DesignPhoto =PhotoImage(file = (self.pwd+ "/media/RosettaLogo.gif"))
        self.Photo = Label(master=self, image=DesignPhoto)
        self.Photo.image = DesignPhoto



    def grid_GUI_objects(self):

        ############# LOOPS #############
        #self.label_Loop.grid(row=11, column=0, columnspan=2, pady=15)

        self.loops_listbox.bind("<Double-Button-1>", lambda event: self.color_region_from_listbox())
        self.loops_listbox.bind('<ButtonRelease-1>', lambda event: self.insert_region_into_entry_from_listbox())
        self.StartLoopLabel.grid(row=16, column=0); self.StartLoopEntry.grid(row=16, column=2)
        self.EndLoopLabel.grid(row=17, column=0); self.EndLoopEntry.grid(row=17, column=2)
        self.ChainIDLabel.grid(row=18, column=0); self.ChainIDEntry.grid(row=18, column=2)
        self.AddLoopButton.grid(row=20, column=2, sticky = W+E)
        self.RmLoopButton.grid(row=21, column=2, sticky=W+E)
        self.loops_listbox.grid(row=22, column=0, rowspan=6, columnspan = 1, sticky = W+E, padx=3);
        self.loops_scroll.grid(row=22, column=1, rowspan=6, sticky=E+N+S)
        ################################


        ########### SEQUENCE ###########
        self.ShoSeqButton.grid(row=19, column=2, sticky=E+W)
        ################################

        self.Photo.grid(row=22, column=2, rowspan=6, columnspan=1, sticky=W+E, padx=3)



#### FUNCTIONS ####
    ##ToDO:  Next refactor will have a map be the main region container.  [loop_string]:region

    def load_loop(self):
        loops_as_strings = self.input_class.load_loop()
        if not loops_as_string:return
        for loop_string in loops_as_strings:
            loop_stringSP = loop_string.split(":")
            self.input_class.region_start.set(loop_stringSP[0])
            self.input_class.region_end.set(loop_stringSP[1])
            self.input_class.region_chain.set(loop_stringSP[2])
            self.add_region()



    def add_region(self):
        """
        Adds region to loop_string and regions.
        Sets GUIInput to have a region variable to indicate last added/current region.
        """

        if not self.check_region_can_be_added(False): return

        looFull = self.input_class.region_start.get()+ ":"+ self.input_class.region_end.get()+":"+self.input_class.region_chain.get().upper()



        region = self.return_region_from_entry()
        if not region.region_exists(self.pose):return
        self.input_class.region_sequence.set(region.get_sequence(self.input_class.pose))
        self.input_class.region = region

        #Add Region to Regions
        self.input_class.regions.add_region(region)
        self.input_class.loops_as_strings.append(looFull)
        self.loops_listbox.insert(END, looFull)

        if self.toolkit.pymol_class.auto_send_region_colors.get():
            self.color_region_from_entry()

    def remove_region(self):
        try:
            #Current
            self.input_class.loops_as_strings.remove(self.loops_listbox.get(self.loops_listbox.curselection()))

            #Replacement of loops_as_string
            region_string = self.loops_listbox.get(self.loops_listbox.curselection())

            self.input_class.regions.remove_region(region_string)

            self.loops_listbox.delete(self.loops_listbox.curselection())
            self.input_class.region_sequence.set("")

        except TclError:
            print "Please select a region to remove."
            return



    def check_region_can_be_added(self, silence=True):
        """
        Makes sure region is set in entry boxes.
        """

        if not self.pose.total_residue():
            print "Please load a pose"
            return False

        elif not self.input_class.region_chain.get():
            if not silence:
                print "Missing chain component of the region"
            return False

        else:
            try:
                int(self.input_class.region_start.get())
            except ValueError:

                self.input_class.region_start.set("")
                if self.input_class.region_end.get() and not silence:
                    print "Assuming N-Terminal Region."


            try:
                int(self.input_class.region_end.get())
            except ValueError:

                self.input_class.region_start.set("")
                if self.input_class.region_start.get() and not silence:
                    print "Assuming C-Terminal Region."

        if not self.input_class.region_start.get() and not self.input_class.region_end.get() and not silence:
            print "Assuming Full Chain"


        return True

    def return_region_from_entry(self):
        """
        Returns a region of the currently set region parameters in the region entry boxes.
        """
        if not self.input_class.region_chain.get():return
        looFull = self.input_class.region_start.get()+ ":"+ self.input_class.region_end.get()+":"+self.input_class.region_chain.get().upper()


        return gen_tools.loop_string_to_region(looFull)

    def insert_region_into_entry_from_listbox(self):
        try:
            region_string = self.loops_listbox.get(self.loops_listbox.curselection())
            regionSP = region_string.split(":")
            self.input_class.region_start.set(regionSP[0])
            self.input_class.region_end.set(regionSP[1])
            self.input_class.region_chain.set(regionSP[2])
        except TclError:
            pass

    def show_sequence(self):
        if not self.pose.total_residue():return
        if not self.input_class.region_chain.get():
            self.input_class.region_sequence.set(self.input_class.pose.sequence())

        else:
            if not self.check_region_can_be_added(): return
            region = self.return_region_from_entry()
            if not region.region_exists(self.pose):return
            self.input_class.region_sequence.set(region.get_sequence(self.input_class.pose))

    def color_region_from_entry(self):
        """
        Colors region in entry box.
        """

        if not self.check_region_can_be_added(): return
        region = self.return_region_from_entry()
        if not region.region_exists(self.pose):return
        self.toolkit.pymol_class.color_region(region.get_rosetta_start(self.pose), region.get_rosetta_end(self.pose))

    def color_region_from_listbox(self):
        """
        Colors region from selection in listbox.
        """
        self.insert_region_into_entry_from_listbox()
        self.color_region_from_entry()

    def reinit_regions_into_listbox(self):
        """
        Uses loops_as_strings to reintroduce regions into the listbox upon new pose.
        Main method here called by GUIInput.
        """
        self.loops_listbox.delete(0, END)

        for loop_string in self.input_class.loops_as_strings:
            self.loops_listbox.insert(END, loop_string)

    def color_residue(self):
        resnum = self.input_class.residue_rosetta_resnum.get()
        self.toolkit.pymol_class.color_residue(int(resnum))

    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported/ran from anywhere.
        """

        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP
