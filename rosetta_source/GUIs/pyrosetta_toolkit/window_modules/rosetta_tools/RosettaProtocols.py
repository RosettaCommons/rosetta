#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/RosettaProtocols.py
## @brief  Main window for settup up Rosetta config files.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkMessageBox
import RosettaSetup
import RosettaClusterSetup
#import tools.pdbs
import tkFileDialog
import tkSimpleDialog
import pickle
import os
import tkFont
import re
import tools
from ..clean_pdb import fix_pdb_window



class ProtocolSetup():
    '''
    This allows you to setup a rosetta protocol.  Save it's options file, or run it.
    Very useful to see all the applications, see possible options, and see descriptions for each option.
    However, manually currated, and the whole program needs to be rewritten to have it otherwise...
    
    '''
    
    
    def __init__ (self, main):
        self.main = main
        self.pwd = (self.location())[0]
        self.defaultdir = tools.loaddefaultdir(self.pwd)

        self.basic_OPTIONS = [
            '--Packing--',
            '--Score--',
            '--Refine--',
            '--Run--',
            '--Input--',
            '--Output--',
            '--Docking--',
            '--Loops--'
            
        ]
        self.doc_types = [
            'Description',
            'References',
            'MetaData',
            'Examples',
            'Algorithm',
            'Inputs',
            'Outputs',
            'Tips',
            'Limitations',
            'Analysis'
        ]
        
        self.specifOPTIONS = dict() #?
        self.appDOC = dict() #?
        self.chosen_doc_type = StringVar(); #Which doc type user has chosen to look at
        self.chosen_doc_type.set(self.doc_types[0])
    
    ####Initialize Objects
        self._init_main_variables()
        self._init_objects_and_settings()
        self.loadRosettaSettings()
        
    def _init_main_variables(self):

        self.appRoot = StringVar()
        self.info_type = StringVar(); #Type of the info that is shown.  Part of the radiobutton on the menu.
        self.info_type.set("currated"); #Default type.
        #AppList and Settings
        self.array_of_applications = []
        self.last_app_clicked = StringVar(); #This keeps track of the app which was last clicked...
        self.option = StringVar(); #Set option in OptionMenu
        self.option.set("Options")
        
        #check_button_ck if setting file exist.  If so, load into variables.
        
        
    def _init_objects_and_settings(self):
        self.setup_pathbuilder_objects()
        
    def setTk(self):
        self.listbox_applications = Listbox(self.main)
        self.option_menu_options = OptionMenu(self.main, self.option, "Options")
        self.button_add_option = Button(self.main, text = "Add Option", command = lambda: self.textbox_cmd_options.insert(1.0,self.option.get()+"\n"))
        
        self.button_save_config = Button(self.main, text = "Save Configuration", command = lambda: self.saveConfiguration())
        self.button_load_config = Button(self.main, text = "Load Configuration", command = lambda: self.loadConfiguration())
        self.button_run_config  = Button(self.main, text = "Run  Configuration", command = lambda: self.runConfiguration())
        
        self.button_show_all_options = Button(self.main, text="Show all options for app", command = lambda:os.system(self.application_directory.get()+"/"+self.last_app_clicked.get()+'.'+self.appRoot.get()+ " -help"))
    #### PHOTO ####
        EngPhoto =PhotoImage(file = (self.pwd+"/RosettaLogo.gif"))
        self.Photo = Label(self.main, image=EngPhoto)
        self.Photo.image = EngPhoto
        
    #### Documentation ####
        self.option_menu_doc_types = OptionMenu(self.main, self.chosen_doc_type, *self.doc_types)
        self.author = Button(self.main, text = "Author", command = lambda:self.textbox_cmd_optionsHelp.insert(1.0, self.appDOC[self.last_app_clicked.get()]['AUTHOR']+"\n\n"))
        
        self.check_button_path_builder=Checkbutton(self.main, variable=self.path_builder_check, text=" Show Path Builder?")
        self.documentation_frame = Frame(self.main, bd=3, relief=GROOVE);helpfont = tkFont.Font(size=12)
        self.documentation_textbox = Text(self.documentation_frame,wrap="word", height = 20, width=55, background = 'white', font = helpfont)

    #### Config ####
        scroll_cmd_options = Scrollbar(self.main)
        self.frame_cmd_options = Frame(self.main, bd=4, relief=GROOVE)
        self.textbox_cmd_options = Text(self.frame_cmd_options, height=8, width=55)
        self.textbox_cmd_options.configure(yscrollcommand=scroll_cmd_options.set)
        scroll_cmd_options.config(command = self.textbox_cmd_options.yview)
        
    def shoTk(self, r, c):
        '''
        Grid the Tk objects to the window.
        '''
        self._r_ = r; self._c_ = c

        self.Photo.grid(row =self._r_, column = self._c_+4)
        self.listbox_applications.grid(row = self._r_, column = self._c_, sticky = W+E)
        #self.optionEntry.grid(column = self._c_+1, row = self._r_, columnspan = 3)
        self.option_menu_options.grid(column = self._c_+1, row = self._r_, sticky = W+E)
        self.button_show_all_options.grid(column=self._c_+1, row=self._r_+1, sticky=W+E, columnspan=2)
        self.button_add_option.grid(column=self._c_+2, row = self._r_, sticky = W+E)
        self.button_save_config.grid(column = self._c_, row = self._r_+1, sticky = W+E)
        self.button_load_config.grid(column = self._c_, row = self._r_+2, sticky = W+E)
        self.button_run_config.grid(column = self._c_, row = self._r_+3, sticky = W+E)
        #self.input.grid(column=self._c_, row=self._r_+4, sticky = W+E)
        #self.output.grid(column=self._c_, row=self._r_+5, sticky=W+E)
        
        self.author.grid(column = self._c_+8, row = self._r_, sticky=W+E)
        self.option_menu_doc_types.grid(column = self._c_+5, columnspan = 2, row = self._r_+5)
        self.textbox_cmd_options.grid(row=self._r_+1, column=self._c_+1, columnspan=4, rowspan=4, sticky=W+E)
        self.frame_cmd_options.grid(row = self._r_+1, column=self._c_+1, columnspan =4, rowspan=4, sticky=W+E)
        self.documentation_textbox.grid(row=self._r_, column=self._c_+5, columnspan = 2, rowspan = 5, sticky = W+E)
        self.documentation_frame.grid(row = self._r_, column=self._c_+5, columnspan=2, rowspan = 5, sticky = W+E, padx=6, pady=7)
        
        self.check_button_path_builder.grid(row = self._r_+6, column = self._c_)
        self.listbox_applications.bind("<ButtonRelease-1>", lambda event: self.populateOptionMenu(self.listbox_applications.get(self.listbox_applications.curselection())))
        
        #Initialize everything first:
        self.sho_manually_currated_info()
        self.set_tracers()
    def setMenu(self, main):
        self.main = main
        self.MenBar = Menu(self.main)
        self.Help = Menu(self.main, tearoff=0)
        #self.Help.add_command(label = "General")
        #self.Help.add_command(label = "Running")
        #self.Help.add_command(label = "Parallel Runs")
        self.Help.add_command(label = "List all Possible Options for app", command = lambda: os.system(self.application_directory.get()+"/"+self.last_app_clicked.get()+'.'+self.appRoot.get()+ " -help"))
        self.MenBar.add_cascade(label = "Help", menu = self.Help)
        
        self.MenRepopulate = Menu(self.main, tearoff=0)
        self.MenRepopulate.add_radiobutton(label = "APP", variable = self.info_type, value="app")
        self.MenRepopulate.add_radiobutton(label = "ALL", variable = self.info_type, value="all")
        self.MenRepopulate.add_radiobutton(label = "DOXYGEN", variable = self.info_type, value="doxygen")
        self.MenRepopulate.add_radiobutton(label = "Currated", variable = self.info_type, value="currated")
        self.MenBar.add_cascade(label = "RePopulate", menu=self.MenRepopulate)
        self.Setup = Menu(self.main, tearoff=0)
        #Does this work?
        self.Setup.add_command(label = "Rosetta Paths", command = lambda: self.loadRosettaSettings())
        self.Setup.add_command(label = "Set default directory", command = lambda: tools.setdefaultdir(self.defaultdir, self.pwd))
        self.Setup.add_separator()
        self.Setup.add_command(label = "Setup PDB(s) for Rosetta", command = lambda: fix_pdb_window.FixPDB().runfixPDBWindow(self.main, 0, 0))
        self.Setup.add_command(label = "Make list of PDBs from a directory", command = lambda: self.makeList())
        
        self.MenBar.add_cascade(label = "Setup", menu=self.Setup)
        self.Cluster = Menu(self.main, tearoff=0)
        self.Cluster.add_command(label = "Setup Cluster Run", foreground='red', command = lambda: self.shoClusterSetup())
        
        self.Cluster.add_command(label = "Start multi-core Run")
        self.MenBar.add_cascade(label = "Parallel", menu=self.Cluster)
        self.main.config(menu=self.MenBar)
        
    def set_tracers(self):
        '''
        Sets any tracers that the program needs to be aware of.
        '''
        self.option.trace_variable('w', self.helpCallback)
        self.chosen_doc_type.trace_variable('w', self.appCallback)
        self.path_builder_check.trace_variable('w', self.pathBuildCallback)
        self.info_type.trace_variable('w', self.info_typeCallback)
#### CALLBACK ####
    def helpCallback(self, name, index, mode):
        '''
        More like option call back.  Inserts option info into documentation_textbox.
        '''
        
        varValue = self.option.get()
        try:
            self.documentation_textbox.insert(1.0, self.specifOPTIONS[self.last_app_clicked.get()][varValue]+"\n\n")
        except KeyError:
            pass
        
    def appCallback(self, name, index, mode):
        varValue = self.chosen_doc_type.get(); varValue = varValue.upper()
        try:
            #Quick fix as placeholder.
            if self.info_type.get()=="app":
                self.sho_app_help_options()
            else:
                self.documentation_textbox.insert(1.0, self.appDOC[self.last_app_clicked.get()][varValue]+"\n\n")
        except KeyError:
            self.documentation_textbox.insert(1.0, "Documentation not found...\n\n")
        except AttributeError:
            pass
        
    def pathBuildCallback(self, name, index, mode):
        varValue = self.path_builder_check.get()
        if self.path_builder_check.get()==0:
            self.delPathBuilder()
        else:
            self.shoPathBuilder()
    
    def info_typeCallback(self, name, index, mode):
        '''
        Callback for the info type - currated, doxygen info, options_rosetta, etc.
        '''
        varValue = self.info_type.get()
        if varValue == "currated":
            self.sho_manually_currated_info()
        elif varValue == "all":
            self.sho_all_rosetta_options()
        elif varValue == "doxygen":
            self.sho_doxygen_documentation()
        elif varValue == "app":
            self.sho_app_help_options()
        else:
            print "Callback failed..."
        return
    
    def populateOptionMenu(self, app):
        '''
        What happens after you double click an application.
        '''
        
        self.last_app_clicked.set(app)
        self.option_menu_options["menu"].delete(0, END)
        apOPTIONS = []
        if self.specifOPTIONS.has_key(app):
            for keys in self.specifOPTIONS[app]:
                apOPTIONS.append(keys)
            apOPTIONS.sort()
            for i in apOPTIONS:
                self.option_menu_options["menu"].add_command(label=i, command=lambda temp = i: self.option_menu_options.setvar(self.option_menu_options.cget("textvariable"), value = temp))
        else:
            print "No Options found.  You are on your own...."
            noneList = ["Not Found","Refer to RosettaCommons"]
            for i in noneList:
                self.option_menu_options["menu"].add_command(label=i, command=lambda temp = i: self.option_menu_options.setvar(self.option_menu_options.cget("textvariable"), value = temp))
        #This is where we use put the description of the protocol into the menu.
        try:
            self.documentation_textbox.insert(1.0, self.appDOC[app]['DESCRIPTION']+"\n\n")
        except KeyError:
            self.documentation_textbox.insert(1.0, "No Documentation Found\n\n")
            
#### SHOW OPTION TYPES #### 
    def sho_manually_currated_info(self):
        '''
        Loads and shows manually currated info if possible.
        '''
        self.specifOPTIONS = pickle.load(open(self.pwd+"/Rosetta3-3.p")); #APP:Descriptions
        self.appDOC = pickle.load(open(self.pwd+"/Rosetta3-3Apps.p")); #No Idea
        self.read_applications_from_directory(self.application_directory.get()); #Populate array_of_applications
        self.populate_applications(self.array_of_applications)
        
    def sho_all_rosetta_options(self):
        '''
        Loads data from options_rosetta.py
        '''
        pass
        #Need to figure out how to parse the options class.
        option_path = os.path.join(self.source_directory.get(), "src", "basic", "options","options_rosetta.py")
        ROSETTA_OPTIONS = open(option_path, 'r')
        
        ROSETTA_OPTIONS.close()
        pass
    
    def sho_doxygen_documentation(self):
        '''
        Shows doxygen_documentation
        '''
        pass
    
    def sho_app_help_options(self):
        '''
        uses app -help to parse and print all info.
        Repopulates the specifOPTIONS dictionary.
        '''
        #Note, I have heard popopen fails in 3.0 and less. So -
        
        print "Reading all available options for each application.  This may take a few minutes..."
        #Read applications, parse each app into dictionary.
        self.read_applications_from_directory(self.application_directory.get())
        try:
            self.app_help_options
            self.populateOptionMenu(self.last_app_clicked.get())
            return
        except AttributeError:
            pass
        for app in self.array_of_applications:
            app_path = self.application_directory.get()+"/"+app+'.'+self.appRoot.get()
            if os.path.exists(app_path):
                os.system(app_path+" -help > temp_options.txt")
                OPTIONS = open("temp_options.txt", 'r')
                self.specifOPTIONS[app] = dict(); #Reset the option dictionary.
                option_type = ""
                for line in OPTIONS:
                    line = line.strip()
                    lineSP = line.split("|")
                    if len(lineSP)<3 or re.search("option group", lineSP[3]) or re.search("Option Setting", lineSP[3]):
                        continue
                    
                    #Here we make the optiongroup with the option.
                    if re.search(":", lineSP[0]):
                        option_type = lineSP[0].strip()

                    opt = option_type+lineSP[0].strip()+" "+lineSP[1].strip()
                    desc =lineSP[2].strip()+" "+lineSP[3].strip()
                    self.specifOPTIONS[app][opt]=desc
                OPTIONS.close()
                os.system('rm temp_options.txt')
        self.populateOptionMenu(self.last_app_clicked.get())
        self.app_help_options = self.specifOPTIONS ; #This is so that we do not have to reload.
    #### FUNCTIONS ####

    def read_applications_from_directory(self, Apppath):
        '''
        Reads all applications from the directory.
        Seperates .linuxgcc.
        -appends the appLIST, sets appRoot
        -Adds the basic options in globalOPTIONS to appLIST
        -sorts the appLIST.
        **DOES NOT UPDATE LISTBOX, etc.**
        '''
        
	#This exception handling does not work!  Fix!
        if os.path.exists(Apppath):
            list = os.listdir(Apppath)
        else:
            os.system('rm '+self.pwd+"/RosettaSettings.txt")
            self.loadRosettaSettings
            list = os.listdir(Apppath)
        for apps in list:
            appsSP = apps.split(".")
            if appsSP[1]=="default":
                self.array_of_applications.append(appsSP[0])
                self.appRoot.set(appsSP[2])
        for app in self.basic_OPTIONS:
            self.array_of_applications.append(app)
        self.array_of_applications.sort()
    

    def populate_applications(self, array):
        '''
        Populates the application listbox with an array of strings.
        '''
        
        for string in array:
            self.listbox_applications.insert(END, string)
            
        
    def loadRosettaSettings(self):
        '''
        Load Rosetta settings text file.
        '''
        
        setup_window = RosettaSetup.SetupRosettaPaths(self.main)
        if not setup_window.result:
            print "Please check paths and try again..."
            second_setup_window = RosettaSetup.SetupRosettaPaths(self.main)
            if not second_setup_window.result:
                print "Directories are still wrong...exiting..."
                exit()

        self.database_directory = StringVar()
        self.application_directory = StringVar()
        self.source_directory = StringVar()
        
        FILE = open(self.pwd+"/RosettaSettings.txt", 'r')
        for line in FILE:
            lineSP = line.split()
            if lineSP[0]=="DATABASE":
                self.database_directory.set(lineSP[1])
            elif lineSP[0]=="APPLICATIONS":
                self.application_directory.set(lineSP[1])
            elif lineSP[0]=="ROSETTA_SOURCE":
                self.source_directory.set(lineSP[1])
        FILE.close()
            
            
           
    
    def saveConfiguration(self, save_as_temp = False):
        '''
        Save cmd configuration to file.
        if save_as_tmp, does not do the file dialog.
        '''
        
        # Add database and app path to config before saving.
        #check_button_cks to make sure the correct app is selected, so we don't screw up when running on the cluster.
        if tkMessageBox.askquestion(message="Is this the currect application:  "+self.last_app_clicked.get()+"?", default=tkMessageBox.YES)=="no":
            print "Save file aborted, returning..."
            return
        if not save_as_temp:
            FILE = tkFileDialog.asksaveasfile(initialdir = self.defaultdir)
        else:

            FILE = open(self.pwd+"/temp_settings_"+self.last_app_clicked.get()+".txt", 'w')
        config = self.textbox_cmd_options.get(1.0, END)
        
        config = '#'+self.application_directory.get()+"/"+self.last_app_clicked.get()+'.'+self.appRoot.get()+"\n"+'-in:path:database '+self.database_directory.get()+"\n"+config
        #config = self.toolKitInterfaceConfigChange(config)
        FILE.write(config)
        FILE.close()
        return
    
    def loadConfiguration(self):
        '''
        Load a rosetta cmd config file.
        '''
        
        FILE = tkFileDialog.askopenfile(initialdir = self.defaultdir)
        config = FILE.read()
        config = config.strip()
        #Parse config, take out database and app type. Set curselection to app type.
        configSP = config.split("\n")
        apppath = configSP[0]; configSP.pop(0)
        app = os.path.split(apppath)[1].split('.')[0]
        app = app.replace("#", "")
        print app
        for stuff in configSP:
            if re.search("database", stuff):
                ind = configSP.index(stuff)
                configSP.pop(ind)
                break
        config = "\n".join(configSP)
        
        
        #set curselection to app type:
        self.listbox_applications.selection_set(self.array_of_applications.index(app))
        #self.listbox_applications.selection_set()
        self.textbox_cmd_options.insert(1.0, config)
        return
    
    def runConfiguration(self):
        '''
        Used for quick run.  Like conversions, scoring, etc.
        popopen doesn't work well, so just freezes system? Needs work.
        '''
        self.saveConfiguration(True)
        os.system(self.application_directory.get()+"/"+self.last_app_clicked.get()+'.'+self.appRoot.get()+" @"+self.pwd+"/temp_settings_"+self.last_app_clicked.get()+".txt")
#### Other Windows ####

    def setup_pathbuilder_objects(self):
        '''
        Initialized pathbuilder Tk objects.
        '''
        self.root = StringVar()
        self.filename = StringVar()
        self.path_builder_check = IntVar(); #variable for if path builder is shown or not.
        self.path_builder_check.set(0)
        self.rootLabel = Label(self.main, text = "Root")
        self.rootEntry = Entry(self.main, textvariable = self.root)
        self.openDIR = Button(self.main, text="Find Root", command = lambda: self.root.set(self.getRoot()))
        self.fileLabel = Label(self.main, text = "FileName")
        self.fileEntry = Entry(self.main, textvariable = self.filename)
        self.insertBoth = Button(self.main, text="Insert Both", command = lambda: self.insertDIR_File())
        self.insertDIR  = Button(self.main, text="Insert DIR", command = lambda: self.insertDIR_())
        self.insertOpen = Button(self.main, text="Search + Insert File", command = lambda: self.insertSearch())
        self.pathWidgets = [
            self.rootLabel,
            self.rootEntry,
            self.openDIR,
            self.fileLabel,
            self.fileEntry,
            self.insertBoth,
            self.insertDIR,
            self.insertOpen,
        ]
        
    def shoPathBuilder(self):
        self.rootLabel.grid(row = self._r_+5, column=self._c_+2, sticky=W+E)
        self.rootEntry.grid(row = self._r_+5, column=self._c_+1)
        self.openDIR.grid(row = self._r_+5,column=self._c_+3)
        self.fileLabel.grid(row = self._r_+6, column=self._c_+2)
        self.fileEntry.grid(row = self._r_+6, column=self._c_+1)
        self.insertBoth.grid(row=self._r_+7, column=self._c_+2, sticky=W+E)
        self.insertDIR.grid(row = self._r_+7, column=self._c_+1, sticky=W+E)
        self.insertOpen.grid(row = self._r_+8, column=self._c_+1, columnspan=2, sticky=W+E)
    
    def delPathBuilder(self):
        for widget in self.pathWidgets:
            widget.grid_forget()
    def getRoot(self):
        root = tkFileDialog.askdirectory(initialdir = self.defaultdir)
        return root
    
    def insertDIR_File(self):
        self.textbox_cmd_options.insert(1.0, " "+self.root.get()+"/"+self.filename.get())
        
    def insertDIR_(self):
        self.textbox_cmd_options.insert(1.0, " "+self.root.get())
        
    def insertSearch(self):
        filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
        self.textbox_cmd_options.insert(1.0, " "+filepath)
            
    def shoClusterSetup(self):
        WinClusSetup = Toplevel(self.main)
        WinClusSetup.title("Cluster Setup")
        Clus = RosettaClusterSetup.ClusterSetup()
        Clus.shoWindow(WinClusSetup, 0, 0)
        
        
#### OTHER FUNCTIONS ####

    def makeList(self):
        '''
        Makes a list of PDB's from a directory.  Does not walk directory.  This can be an option later.
        '''
        
        directory = tkFileDialog.askdirectory(initialdir = self.defaultdir)
        if directory ==None:
            return
        FILES = os.listdir(directory)
        NEWFILE = open(directory+"/PDBLIST.txt", 'w')
        for names in FILES:
            if re.search(".pdb", names) and (re.search("\._", names)== None) and (re.search("~", names)== None):#This shows how stupid python/linux can be...:
                print names
                p = os.path.join(directory, names)
                NEWFILE.write(p+"\n")
        print "Written..."
        print "File saved as 'PDBLIST.txt' in directory specified."
        NEWFILE.close()
        return
    
    def makeRecursiveList(self, directory):
        '''
        Makes a list of PDB's from a directory.  Walks Directory
        DOES NOT WORK>
        '''
        directory = tkFileDialog.askdirectory(initialdir = self.defaultdir)
        filenum = 1
        if contains=="all":
            contains = ""
        #print outFolder
        AllFiles = []
        for root, dirs, files in os.walk(inFolder, topdown=True):
            #print "Root" + root
            for f in files:
                if re.search(contains, f):
                    print "File_"+repr(filenum)+"_"+f
                    p = os.path.join(root, f)
                    AllFiles.append(p)
        
    def getfile(self):
        '''
        simply gets a filename...wish I could do this within the command syntax of tkinter...
        '''
        filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
        return filepath
    
    def location(self):
        '''
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        '''
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP

if __name__ == '__main__':
    MainWindow = Tk()
    MainWindow.title("Window to Rosetta")
    SetupWindow = ProtocolSetup(MainWindow)
    SetupWindow.setTk()
    SetupWindow.shoTk(0, 0)
    SetupWindow.setMenu(MainWindow)
    MainWindow.mainloop()
