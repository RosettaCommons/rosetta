#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/RosettaSetup.py
## @brief  Dialog window that loads when rosettaprotocols runs.  Used to set default paths.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkSimpleDialog
import tkFileDialog
import tkMessageBox
import os
import os.path





class SetupRosettaPaths(tkSimpleDialog.Dialog):
    """
    Window to setup rosetta paths for the system.
    Saves a settings file, which is used by subsequent rosetta tools.
    """
            
    def body(self, main):
        #INIT
        self.pwd = (self.location())[0]
        self.database = StringVar()
        self.applications = StringVar()
        self.source = StringVar()
        #self.fragmentpicker = StringVar()
    
        #check_button_ck if setting file exist.  If so, load into variables.
        if os.path.exists(self.pwd+"/PATHSETTINGS.txt"):
            FILE = open(self.pwd+"/PATHSETTINGS.txt", 'r')
            self.option_array = {"DATABASE":self.database, "APPLICATIONS":self.applications, "ROSETTA_SOURCE":self.source}
            for line in FILE:
                lineSP = line.split()
                if self.option_array.has_key(lineSP[0]):
                    variable = self.option_array[lineSP[0]]
                    variable.set(lineSP[1])
                else:
                    print lineSP[0]+" not found in file. Please fix manually."
                    FILE.close()
                    exit()
            FILE.close()
            
        #Check each:
        self.checkPaths(False); #False sets self.result to False.
        
        
        self.main = main
        #self.main.title("Setup Rosetta Paths")
        self.column = 0
        self.row = 0
        
        #LABELS
        self.databaseLabel = Label(self.main, text = "Path to rosetta_database:")
        self.sourceLabel = Label(self.main, text = "Path to rosetta_source:")
        self.applicationLabel = Label(self.main, text = "Path to Applications:")
        #self.fragmentLabel = Label(self.main, text = "Path to Fragment Picker:")
        #ENTRY
        self.databaseEntry = Entry(self.main, textvariable = self.database)
        self.sourceEntry = Entry(self.main, textvariable = self.source)
        self.applicationEntry = Entry(self.main, textvariable = self.applications)
        #self.fragmentEntry = Entry(self.main, textvariable = self.fragmentpicker)
        #BUTTON
        self.databaseButton = Button(self.main, text = "Choose Path", command = lambda:self.database.set(tkFileDialog.askdirectory(initialdir=self.pwd, title="Choose rosetta_database Directory")))
        self.applicationButton = Button(self.main, text = "Choose Path", command = lambda:self.applications.set(tkFileDialog.askdirectory(initialdir = self.pwd, title = "Choose Application Directory")))
        self.sourceButton = Button(self.main, text = "Choose Path", command = lambda: self.source.set(tkFileDialog.askdirectory(initialdir = self.pwd, title = "Choose rosetta_source Directory")))
        #self.fragmentButton = Button(self.main, text = "Choose Path", command = lambda:self.fragmentpicker.set(tkFileDialog.askdirectory(initialdir=self.pwd, title="Choose Fragment Picker Directory")))
        #SAVE
        self.saveButton = Button(self.main, text = "Load Settings", command = lambda: self.loadSettings())
        
        
        #Sets up menu.
        self.sourceLabel.grid(column=self.column, row=self.row, sticky=W+E)
        self.sourceEntry.grid(column=self.column+1, row=self.row, sticky=W+E)
        self.sourceButton.grid(column=self.column+2, row=self.row, stick=W+E)
        
        self.databaseLabel.grid(column = self.column, row = self.row+1, sticky = W+E)
        self.databaseEntry.grid(column = self.column+1, row = self.row+1, sticky = W+E)
        self.databaseButton.grid(column =self.column+2, row = self.row+1, sticky = W+E)
        

        self.applicationLabel.grid(column = self.column, row=self.row+2, sticky = W+E)
        self.applicationEntry.grid(column = self.column+1, row = self.row+2, sticky = W+E)
        self.applicationButton.grid(column =self.column+2, row = self.row+2, sticky = W+E)

        #self.fragmentLabel.grid(column = self.column, row=self.row+2, sticky = W+E)
        #self.fragmentEntry.grid(column = self.column+1, row = self.row+2, sticky = W+E)
        #self.fragmentButton.grid(column = self.column+2, row=self.row+2, sticky = W+E)
        
        self.saveButton.grid(column = self.column, row = self.row+3, columnspan = 2, sticky = W+E)
        
        EngPhoto =PhotoImage(file = (self.pwd+"/../media/RosettaLogo.gif"))
        self.Photo = Label(self.main, image=EngPhoto)
        self.Photo.image = EngPhoto
        self.Photo.grid(row =self.row, column = self.column+3, rowspan = 5)
        
        
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        """
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP
    
    def loadSettings(self):
        """
        Loads settings from an alternative location, saves rosettasettings.txt.
        """
        f = tkFileDialog.askopenfilename(initialdir=self.pwd)
        if not f:return
        FILE = open(f, 'r')
        for line in FILE:
            FILE = open(f, 'r')
            for line in FILE:
                lineSP = line.split()
                if self.option_array.has_key(lineSP[0]):
                    variable = self.option_array[lineSP[0]]
                    variable.set(lineSP[1])
                else:
                    print lineSP[0]+" not found in file. Please fix manually."
                    FILE.close()
                    exit()
        FILE.close()
        
    def apply(self):
        self.checkPaths()
        FILE = open(self.pwd+"/PATHSETTINGS.txt", mode = 'w')
        FILE.write("ROSETTA_SOURCE\t"+self.source.get()+'\n')
        FILE.write("DATABASE\t"+self.database.get()+'\n')
        FILE.write("APPLICATIONS\t"+self.applications.get()+'\n')

        #FILE.write("FRAGMENTPICKER\t"+fragmentpicker+'\n')
        print "Settings Saved..."
        FILE.close()
        
    def checkPaths(self, result=True):
        """
        This is now hardcoded for where RosettaPathSetup Resides.  Update this if the file is moved.
        First check Settings file.  If those paths not found, will try to use relative paths.
        """
        self.result = result
        
        if not os.path.exists(self.database.get()):
            print "Using relative path for database"
            self.database.set(self.get_relative_path(self.pwd, 4)+'/rosetta_database')
            if not os.path.exists(self.database.get()):
                self.database.set("NA")
                self.result = False
        if not os.path.exists(self.applications.get()):
            print "Using relative path for applications"
            self.applications.set(self.get_relative_path(self.pwd, 3)+'/bin')
            if not os.path.exists(self.applications.get()):
                self.applications.set("NA")
                self.result = False
        if not os.path.exists(self.source.get()):
            print "Using relative path for source"
            self.source.set(self.get_relative_path(self.pwd, 3))
            if not os.path.exists(self.source.get()):
                self.source.set("NA")
                self.result = False
        
        return
    
    def get_relative_path(self, path, num_back):
        pathSP = path.split("/")
        directories = len(pathSP)
        return "/".join(pathSP[:(directories-num_back)])
        
if __name__ == '__main__':
    MainWindow = Tk()
    SetupWindow = SetupRosettaPaths()
    SetupWindow.shoWindow(MainWindow, 0, 0)
    MainWindow.mainloop()
