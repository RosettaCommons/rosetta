#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/options_system/options_creator.py
## @brief  Window to load and set rosetta global options
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *

#Python Imports
import os

#Tkinter Imports
from Tkinter import *
import tkFileDialog

#Toolkit Imports
from app.pyrosetta_toolkit.window_main import global_variables

class OptionSystemManager:
    def __init__(self):
        """
        Interacts with the rosetta options system by re-initing every time.
        To get and set options That do not have to be always present use:
        (These will be lost when protocols are run as protocols need to reinitialize rosetta if run on multiple processors)
        rosetta.basic.options.get_datatype_option('string_without-')
        rosetta.basic.options.set_datatype_option('string_without-', value)
        Else: use OptionSystemManager.add_option('string')
        """
        self.opts = []
        self.find_rosetta_database()
        #These are expanded options from what PyRosetta __init__ has.
        self.basic_options = ["-ex1", "-ex2aro", "-add_orbitals", "-run:seed_offset 1000"]
        self.extra_options = []
        self.pwd = self.location()[0]
        
        self.current_directory = global_variables.current_directory
        
        self.reset_options()
        if os.path.exists(self.pwd+"/settings.txt"):
            print "Loading personal option settings..."
            self.load_settings(self.pwd+"/settings.txt")
            
        self.common_options = [
            'Enter Custom:',
            '-restore_pre_talaris_2013_behavior',
            '-ex2',
            '-use_input_sc',
            '-linmem_ig 20',
            '-double_lazy_ig',
            '-ignore_unrecognized_res',
            '-out:file:silent',
            '-relax:constrain_relax_to_start_coords',
            '-relax:coord_constrain_sidechain',
            '-relax:dualspace',
            '-nonideal'
            #'-fixbb:minimize_sidechains',
            #'-fixbb:min_pack'
        ]

        self.common_option = StringVar(); self.common_option.set(self.common_options[0])
        self.custom_option = StringVar()
        
    def setTk(self, main):
        main.title("Options System Manager")
        self.main = main

    
    #### Set ####
        self.common_option_menu = OptionMenu(self.main, self.common_option, *self.common_options)
        self.custom_option_entry= Entry(self.main, textvariable=self.custom_option, justify=CENTER)
        self.add_option_button = Button(self.main, text = "Add option", command = lambda: self.extend())
        
    #### Save/Load ####
        self.save_defaults_button = Button(self.main, text = "Save Defaults",command = lambda: self.save_settings(self.pwd+"/settings.txt"))
        self.load_defaults_button = Button(self.main, text="Load Defaults",command = lambda: self.load_settings(self.pwd+"/settings.txt"))
        self.save_button =Button(self.main, text = "Save as...", command = lambda: self.save_settings(""))
        self.load_button=Button(self.main, text = "Load ...", command = lambda: self.load_settings(""))
        self.clear_button=Button(self.main, text = "Clear Defaults", command = lambda: self.clear_defaults())
        self.clear_current_button = Button(self.main, text = "Clear Current", command = lambda: self.reset_options())
        
    ### Callbacks ###
        #self.common_option.trace_variable('w', self.common_option_disable)
        
    def shoTk(self, r=0, c=0):
        self.common_option_menu.grid(row = r, column=c, sticky=W+E)
        self.custom_option_entry.grid(row=r+1, column=c, sticky=W+E)
        self.add_option_button.grid(row=r, column=c+1, rowspan=2, sticky=W+E)
        
        self.save_defaults_button.grid(row=r, column=c+2, sticky=W+E)
        self.load_defaults_button.grid(row=r+1, column=c+2, sticky=W+E)
        self.save_button.grid(row = r, column=c+3, sticky=W+E)
        self.load_button.grid(row= r+1, column=c+3, sticky=W+E)
        self.clear_button.grid(row=r+2, column=c+2, columnspan=2, sticky=W+E)
        self.clear_current_button.grid(row = r+3, column=c+2, columnspan = 2, sticky=W+E)
    
        
    def print_current_options(self):
        print "Options:  "+" ".join(self.opts)
        
    def add_option(self, option_string):
        """
        Used to add an option outside of the GUI.
        """
        self.common_option.set(option_string)
        self.extend()
            
    def common_option_disable(self, name, index, mode):
        """
        Callback for common_option.  If not set to enter custom, disables option entry.
        """
        varValue = self.common_option.get()
        if varValue==self.common_options[0]:
            self.custom_option_entry.config(state=NORMAL)
        else:
            self.custom_option_entry.config(state=DISABLED)
            
    def find_rosetta_database(self):
        '''
        Directly from __init__ file in /rosetta
        Author: Sergey Lyskov
        '''

        self.database = rosetta.rosetta_database_from_env()

            
    def location(self):
        '''
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        '''
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP
    
    def load_settings(self, path):
        if path != self.pwd+"/settings.txt":
            path = tkFileDialog.askopenfilename(initialdir=self.current_directory, title='Option Settings')
            if not path:return
        if not os.path.exists(path):
            print "No settings to load..."
            return
        self.current_directory = os.path.dirname(path)
        settings_opts =[]
        FILE = open(path, 'r')
        for line in FILE:
            line = line.strip()
            settings_opts.append(line)
        self.opts = self.opts+settings_opts
        self.extra_options = self.extra_options+settings_opts

        rosetta.init(" ".join(self.opts))
        FILE.close()
        print "Settings loaded..."
        
    def clear_defaults(self):
        os.system('rm '+self.pwd+"/settings.txt")
        self.reset_options()
        print "Defaults Cleared..."
    
    def save_settings(self, path):
        '''
        Saves settings file which looks like the @settings file.  asks for filename and directory if it is not default
        '''
        
        if path != self.pwd+"/settings.txt":
            path = tkFileDialog.asksaveasfilename(initialdir=self.current_directory, title='Option Settings')
            
        if not path:
            return
        self.current_directory = os.path.dirname(path)
        FILE = open(path, 'w')
        for option in self.extra_options:
            FILE.write(option+"\n")
        FILE.close()
        print "Settings saved...."
    
    def re_init(self):
        """
        Used by multiprocessing protocols to reinit rosetta with arguments held here.
        ANY arguments set by rosetta.basic.options.set_xxx_option() MAY be annihilated.
        """
        rosetta.init(" ".join(self.opts))
        
    def reset_options(self):
        '''
        Resets self.opts and self.args, extends with basic options that are nessessary, initializes rosetta.
        '''
        
        self.opts = []
        print self.basic_options
        self.extra_options = []
        self.opts = self.opts+self.basic_options
        rosetta.init(" ".join(self.opts))
        
    def extend(self):
        '''
        calls extend options.  Used to determine if common option or custom option.
        '''
        if self.common_option.get()==self.common_options[0]:
            self.extend_options(self.custom_option.get())
        else:
            self.extend_options(self.common_option.get())
            
    def extend_options(self, string):
        if not string[0]=='-':string = "-"+string
        
        self.opts.append(string)
        self.extra_options.append(string)
        rosetta.init(" ".join(self.opts))
        print "Options extended.."
        
if __name__ == '__main__':
    test = OptionSystemManager()
    main = Tk()
    test.setTk(main)
    test.shoTk()
    main.mainloop()
