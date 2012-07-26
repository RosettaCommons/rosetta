#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/options_system/options_creator.py
## @brief  Window to load and set rosetta global options
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
import tkFileDialog
import os


class Option_System_Manager:
    def __init__(self, current_directory=False):
        
        self.opts = []
        self.args = rosetta.utility.vector1_string()
        self.find_rosetta_database()
        self.basic_options = ["app", "-database", self.database, "-ex1", "-ex2aro"]
        self.extra_options = []
        self.pwd = self.location()[0]
        
        ## Used due to main Gui
        if not current_directory:
            self.current_dirctory = self.pwd
        else:
            self.current_directory = current_directory
        ## Used due to main Gui
        
        self.reset_options()
        if os.path.exists(self.pwd+"/settings.txt"):
            print "Loading personal option settings..."
            self.load_settings(self.pwd+"/settings.txt")
            
        self.common_options = [
            'Custom:',
            '-dun10',
            '-dun08',
            '-ex2',
            '-ignore_unrecognized_res',
            '-out:file:silent'
        ]

        
    def setTk(self, main):
        self.main = main
        self.common_option = StringVar(); self.common_option.set(self.common_options[0])
        self.custom_option = StringVar()
    
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
    def find_rosetta_database(self):
        '''
        Directly from __init__ file in /rosetta
        Writer: Sergey Lyskov
        '''
        if os.path.isdir('rosetta_database'):
            self.database = os.path.abspath('rosetta_database')
            print 'Found rosetta_database at %s, using it...' % self.database
    
        elif 'PYROSETTA_DATABASE' in os.environ:
            self.database = os.path.abspath( os.environ['PYROSETTA_DATABASE'] )
            print 'PYROSETTA_DATABASE environment variable was set to: %s... using it...' % self.database
    
        elif os.path.isdir(os.environ['HOME'] + '/rosetta_database'):
            self.database = os.path.abspath(os.environ['HOME'] + '/rosetta_database')
            print 'Found rosetta_database at home folder, ie: %s, using it...' % self.database
    
        elif sys.platform == "cygwin" and os.path.isdir('/rosetta_database'):
            self.database = os.path.abspath('/rosetta_database')
            print 'Found rosetta_database at root folder, ie: %s, using it...' % self.database
    
        else:
            print 'Could not find rosetta_database! Check your paths or set PyRosetta environment vars. Exiting...'
            sys.exit(1)
            
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

        self.args.extend(list(settings_opts))
        rosetta.core.init(self.args)
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
        
    def reset_options(self):
        '''
        Resets self.opts and self.args, extends with basic options that are nessessary, initializes rosetta.
        '''
        
        self.opts = []
        print self.basic_options
        self.extra_options = []
        self.opts = self.opts+self.basic_options
        self.args = rosetta.utility.vector1_string()
        self.args.extend(self.opts)
        rosetta.core.init(self.args)
    def extend(self):
        '''
        calls extend options.  Used to determine if common option or custom option.
        '''
        if self.common_option.get()==self.common_options[0]:
            self.extend_options(self.custom_option.get())
        else:
            self.extend_options(self.common_option.get())
            
    def extend_options(self, string):
        self.opts.append(string)
        opts = [string]
        self.extra_options.append(string)
        self.args.extend(opts)
        rosetta.core.init(self.args)
        print "Options extended.."
        
if __name__ == '__main__':
    test = Option_System_Manager()
    main = Tk()
    test.setTk(main)
    test.shoTk()
    main.mainloop()