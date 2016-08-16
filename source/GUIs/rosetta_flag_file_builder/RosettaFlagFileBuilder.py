#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/RosettaFlagFileBuilder.py
## @brief  Main window for settup up Rosetta config files.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import pickle
import os
import tkFont
import re
import tools
import webbrowser
import sys

#Tkinter Imports
from Tkinter import *
import tkMessageBox
import tkFileDialog
import tkSimpleDialog

#Project Imports
from settings import RosettaPathSetup
from modules.DoxygenParser import DoxygenParser
import QsubClusterSetup

class RosettaFlagFileBuilder():
    """
    This allows you to setup a rosetta run.  Save it's options file, or run it.
    Useful to see all the applications, see possible options, and see descriptions for each option.
    Work is under way to enable better parsing of app documentation.
    
    """
    
    
    def __init__ (self, main):
        main.title("Rosetta Flag File Builder")
        self.main = main
        self.main.grid_columnconfigure(0, weight=1)
        
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
            #"Metadata",
            "Code and Demo",
            "References",
            "Purpose",
            "Algorithm",
            "Limitations",
            "Modes",
            "Input Files",
            "Options",
            "Tips",
            "Expected Outputs",
            "Post Processing",
            "new_stuff",
        ]
        
        self.appOPTIONS = dict(); #App specific options [string app][string option]:[string description]
        self.appDOC = dict(); #App specific documentation
        
        self.chosen_doc_type = StringVar(); #Which doc type user has chosen to look at
        self.chosen_doc_type.set("Purpose")
        self.app_binaries = []; #List of binaries in the /bin directory.
    
    ####Initialize Objects####
        self._init_main_variables()
        self._init_objects_and_settings()
        self.__load_rosetta_settings__()
        self.parser = DoxygenParser(self.source_directory.get())
        
    def _init_main_variables(self):

        self.appRoot = StringVar()
        self.info_type = StringVar(); #Type of the info that is shown.  Part of the radiobutton on the menu.
        #AppList and Settings
        self.array_of_applications = []
        self.last_app_clicked = StringVar(); #This keeps track of the app which was last clicked...
        self.option = StringVar(); #Set option in OptionMenu
        self.option.set("Options")
        
        #check_button_ck if setting file exist.  If so, load into variables.
        
        
    def _init_objects_and_settings(self):
        self.setup_pathbuilder_objects()
        
    def setTk(self):
        self.listbox_applications = Listbox(self.main, bd=4, relief=GROOVE)
        self.option_menu_options = OptionMenu(self.main, self.option, "Options")
        self.button_add_option = Button(self.main, text = "Add Option", command = lambda: self.textbox_cmd_options.insert(1.0,self.option.get()+"\n"))
        
        self.button_show_all_options = Button(self.main, text="Show all options for app", command = lambda:os.system(self.application_directory.get()+"/"+self.appDOC[self.last_app_clicked.get()]["AppName"]+'.'+self.appRoot.get()+ " -help"))
    
    #### PHOTO ####
        EngPhoto =PhotoImage(file = (self.pwd+"/media/RosettaLogo.gif"))
        self.Photo = Label(self.main, image=EngPhoto)
        self.Photo.image = EngPhoto
        
    #### Documentation ####
        self.option_menu_doc_types = OptionMenu(self.main, self.chosen_doc_type, *self.doc_types)
        self.author = Button(self.main, text = "Author", command = lambda:self.documentation_textbox.insert(1.0, self.appDOC[self.last_app_clicked.get()]['Metadata']+"\n\n"))
        
        self.check_button_path_builder=Checkbutton(self.main, variable=self.path_builder_check, text=" Show Path Builder?")
        self.documentation_frame = Frame(self.main, bd=3, relief=GROOVE);helpfont = tkFont.Font(size=11)
        self.documentation_textbox = Text(self.documentation_frame,wrap="word", height = 20, width=58, background = 'white', font = helpfont)

    #### Config ####
        scroll_cmd_options = Scrollbar(self.main)
        self.frame_cmd_options = Frame(self.main, bd=4, relief=GROOVE)
        self.textbox_cmd_options = Text(self.frame_cmd_options, height=8, width=55)
        self.textbox_cmd_options.configure(yscrollcommand=scroll_cmd_options.set)
        scroll_cmd_options.config(command = self.textbox_cmd_options.yview)
        
    def shoTk(self, r, c):
        """
        Grid the Tk objects to the window.
        """
        self._r_ = r; self._c_ = c

        self.Photo.grid(row =self._r_, column = self._c_+4)
        self.listbox_applications.grid(row = self._r_, column = self._c_, sticky = W+E)
        #self.optionEntry.grid(column = self._c_+1, row = self._r_, columnspan = 3)
        self.option_menu_options.grid(column = self._c_+1, row = self._r_, sticky = W+E)
        #self.button_show_all_options.grid(column=self._c_, row=self._r_+1+1+1+1, sticky=W+E)
        self.button_add_option.grid(column=self._c_+2, row = self._r_, sticky = W+E)
        
        self.author.grid(column = self._c_+8, row = self._r_, sticky=W+E)
        self.option_menu_doc_types.grid(column = self._c_+5, columnspan = 2, row = self._r_+5)
        self.textbox_cmd_options.grid(row=self._r_+1, column=self._c_+1, columnspan=4, rowspan=4, sticky=W+E)
        self.frame_cmd_options.grid(row = self._r_+1, column=self._c_+1, columnspan =4, rowspan=4, sticky=W+E)
        self.documentation_textbox.grid(row=self._r_, column=self._c_+5, columnspan = 2, rowspan = 5, sticky = W+E)
        self.documentation_frame.grid(row = self._r_, column=self._c_+5, columnspan=2, rowspan = 5, sticky = W+E, padx=6, pady=7)
        
        self.check_button_path_builder.grid(row = self._r_+6, column = self._c_)
        self.listbox_applications.bind("<ButtonRelease-1>", lambda event: self.__populate_option_menu__(self.listbox_applications.get(self.listbox_applications.curselection())))
        
        #Initialize everything first.  Order matters:
        self.info_type.set("doxygen"); #Default type.
        self.__show_doxygen__()
        self.__set_tracers__()
        self.read_applications_from_directory(self.application_directory.get()); #To set self.binaries and approot
        
        
    def setMenu(self, main):
        self.main = main
        self.MenBar = Menu(self.main)
        
        #self.Help.add_command(label = "General")
        #self.Help.add_command(label = "Running")
        #self.Help.add_command(label = "Parallel Runs")

        self.file_menu = Menu(self.main, tearoff=0)
        self.file_menu.add_command(label = 'Save Configuration', command = lambda: self.saveConfiguration())
        self.file_menu.add_command(label = 'Load Configuration', command = lambda: self.loadConfiguration())
        self.file_menu.add_separator()
        self.file_menu.add_command(label = 'Run Configuration', command = lambda: self.runConfiguration())
        self.MenBar.add_cascade(label = "File", menu = self.file_menu)
        self.MenRepopulate = Menu(self.main, tearoff=0)
        self.MenRepopulate.add_radiobutton(label = "All Available Options", variable = self.info_type, value="full")
        self.MenRepopulate.add_radiobutton(label = "Doxygen Documentation", variable = self.info_type, value="doxygen")
        self.MenRepopulate.add_radiobutton(label = "Manually Currated", variable = self.info_type, value="currated")
        self.MenBar.add_cascade(label = "RePopulate", menu=self.MenRepopulate)
        self.Setup = Menu(self.main, tearoff=0)

        self.Setup.add_command(label = "Rosetta Paths", command = lambda: self.__load_rosetta_settings__())
        self.Setup.add_command(label = "Set default directory", command = lambda: tools.setdefaultdir(self.defaultdir, self.pwd))
        self.Setup.add_separator()
        #self.Setup.add_command(label = "Setup PDB(s) for Rosetta", command = lambda: fix_pdb_window.FixPDB().runfixPDBWindow(self.main, 0, 0))
        self.Setup.add_command(label = "Make list of PDBs from a directory", command = lambda: self.makeList())
        
        self.MenBar.add_cascade(label = "Setup", menu=self.Setup)
        self.Cluster = Menu(self.main, tearoff=0)
        self.Cluster.add_command(label = "Setup  Qsub Cluster Run (JD1)", command = lambda: self.shoClusterSetup())
        #self.Cluster.add_command(label = "Setup   MPI Cluster Run (JD2)", command = lambda: self.shoMPIClusterSetup())
        #self.Cluster.add_command(label = "Start multi-core Run", foreground='red')
        self.MenBar.add_cascade(label = "Parallel", menu=self.Cluster)
        self.help = Menu(self.main, tearoff=0)
        
        self.help.add_command(label="Rosetta Manual", command = lambda: webbrowser.open("http://www.rosettacommons.org/manual_guide"))
        self.help.add_command(label="Rosetta Forums", command = lambda: webbrowser.open("http://www.rosettacommons.org/forum"))
        self.help.add_command(label="Rosetta BugTracker", command = lambda: webbrowser.open("http://bugs.rosettacommons.org"))
        self.help.add_separator()
        self.help.add_command(label = "List all Possible Options for app", command = lambda: os.system(self.application_directory.get()+"/"+self.last_app_clicked.get()+'.'+self.appRoot.get()+ " -help"))
        self.MenBar.add_cascade(label = "Help", menu = self.help)
        
        
        self.main.config(menu=self.MenBar)
        
        self.documentation_textbox.insert(1.0, "Efforts are under way to standardize and better parse the formatting, but please refer to rosettacommons for full production runs.\n")
        self.documentation_textbox.insert(1.0, "This application attempts to help in formatting a Rosetta flag file and exploring documentation.  However, numerous applications have special formatting for Doxygen documention, including embeded HTML.\n")
        
    
    def saveConfiguration(self, save_as_temp = False):
        """
        Save cmd configuration to file.
        if save_as_tmp, does not do the file dialog.
        """
        
        # Add database and app path to config before saving.
        #Checks to make sure the correct app is selected and NAMED right due to doxygen bs - so we don't screw up when running.
        app = self.appDOC[self.last_app_clicked.get()]["AppName"]

        if not os.path.exists(self.application_directory.get()+"/"+app+'.'+self.appRoot.get()):
            app= tkSimpleDialog.askstring(title="Continue?", prompt="Application not found.  Please double check name: ", initialvalue=app)
    
        if not app:return
        if not save_as_temp:
            filename = tkFileDialog.asksaveasfilename(initialdir = self.defaultdir)
            if not filename:return
            self.defaultdir = os.path.dirname(filename)
            FILE = open(filename, 'w')
        else:

            FILE = open(self.pwd+"/temp_settings_"+app+".txt", 'w')
        config = self.textbox_cmd_options.get(1.0, END)
        
        config = '#'+self.application_directory.get()+"/"+app+'.'+self.appRoot.get()+"\n"+'-in:path:database '+self.database_directory.get()+"\n"+config
        #config = self.toolKitInterfaceConfigChange(config)
        FILE.write(config)
        FILE.close()
        return app
    
    def loadConfiguration(self):
        """
        Load a rosetta cmd config file.
        If the first line in the file is #[my_application_path], then we set everything in the window to this application.
        """
        
        filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
        if not filepath:return
        self.defaultdir = os.path.dirname(filepath)
        
        FILE = open(filepath, 'r')
        config = FILE.read()
        config = config.strip()
        #Parse config, take out database and app type. Set curselection to app type.
        configSP = config.split("\n")
        apppath = configSP[0]; configSP.pop(0)
        if re.search("#", apppath):
            
            app = os.path.split(apppath)[1].split('.')[0]
        
            app = app.replace("#", "")
            
            app_found = False
            for p in self.appDOC:
                if app==self.appDOC[p]["AppName"]:
                    app = p
                    app_found=True
            
            print app
            
            if app_found:
                self.listbox_applications.selection_set(self.array_of_applications.index(app))
                self.__populate_option_menu__(app)
            
        for option in configSP:
            if re.search("database", option):
                ind = configSP.index(option)
                configSP.pop(ind)
                break
        config = "\n".join(configSP)
        
        self.textbox_cmd_options.delete(1.0, END)
        self.textbox_cmd_options.insert(1.0, config)
        return
    
    def runConfiguration(self):
        """
        Used for quick run.  Like conversions, scoring, etc.
        """
        app = self.saveConfiguration(True)
        #processors = tkSimpleDialog.askinteger(title='processors', prompt="Processors to use", initialvalue=1)
        #if not processors: return
        
        os.system(self.application_directory.get()+"/"+app+'.'+self.appRoot.get()+" @"+self.pwd+"/temp_settings_"+app+".txt")
        os.remove(self.pwd+"/temp_settings_"+app+".txt")
        
    def __set_tracers__(self):
        """
        Sets any tracers that the program needs to be aware of.
        """
        self.option.trace_variable('w', self.__option_doc_callback__)
        self.chosen_doc_type.trace_variable('w', self.__app_doc_callback__)
        self.path_builder_check.trace_variable('w', self.__path_builder_callback__)
        self.info_type.trace_variable('w', self.__doc_type_callback__)

#### CALLBACK ####
    def __option_doc_callback__(self, name, index, mode):
        """
        Inserts option documentation info into documentation_textbox.
        """
        
        varValue = self.option.get()
        
        try:
            self.documentation_textbox.insert(1.0, self.appOPTIONS[self.last_app_clicked.get()][varValue]+"\n\n")
        except KeyError:
            pass
        
    def __app_doc_callback__(self, name, index, mode):
        """
        Repopulate function.
        """
        
        varValue = self.chosen_doc_type.get();
        if self.info_type.get()=="currated":
            varValue = varValue.upper()
        try:
            #Quick fix.
            if self.info_type.get()=="app":
                self.__show_app_help_options__()
            else:
                self.documentation_textbox.insert(1.0, self.appDOC[self.last_app_clicked.get()][varValue]+"\n\n")
        except KeyError:
            self.documentation_textbox.insert(1.0, "Documentation not found...\n\n")
        except AttributeError:
            pass
        
    def __path_builder_callback__(self, name, index, mode):
        """
        Callback for path_builder checkbox
        """
        varValue = self.path_builder_check.get()
        if self.path_builder_check.get()==0:
            self.delPathBuilder()
        else:
            self.shoPathBuilder()
    
    def __doc_type_callback__(self, name, index, mode):
        """
        Callback for the info type variable, set through the menu - currated, doxygen info, options_rosetta, etc.
        """
        varValue = self.info_type.get()
        if varValue == "currated":
            self.__show_manually_currated__()
        elif varValue == "doxygen":
            self.__show_doxygen__()
        elif varValue == "full":
            self.__show_app_help_options__()
        else:
            print "Callback failed..."
        return
    
    def __populate_applications__(self, array):
        """
        Populates the application listbox with an array of strings.
        """
        self.listbox_applications.delete(0, END); #Just added. Double check.
        for string in array:
            self.listbox_applications.insert(END, string)
        
            
        
    def __populate_option_menu__(self, app):
        """
        What happens after you click an application.
        """
        
        self.last_app_clicked.set(app)
        self.option_menu_options["menu"].delete(0, END)
        apOPTIONS = []
        if self.appOPTIONS.has_key(app):
            for keys in self.appOPTIONS[app]:
                apOPTIONS.append(keys)
            apOPTIONS.sort()
            for i in apOPTIONS:
                self.option_menu_options["menu"].add_command(label=i, command=lambda temp = i: self.option_menu_options.setvar(self.option_menu_options.cget("textvariable"), value = temp))
        else:
            print "No Options found.  Refer to RosettaCommons"
            noneList = ["Not Found","Refer to RosettaCommons"]
            for i in noneList:
                self.option_menu_options["menu"].add_command(label=i, command=lambda temp = i: self.option_menu_options.setvar(self.option_menu_options.cget("textvariable"), value = temp))
        #This is where we use put the description of the protocol into the menu.
        try:
            if self.info_type.get()=="currated":
                self.documentation_textbox.insert(1.0, self.appDOC[app]['DESCRIPTION']+"\n\n")
            else:
                self.documentation_textbox.insert(1.0, self.appDOC[app]['Purpose']+"\n\n")
        except KeyError:
            self.documentation_textbox.insert(1.0, "No Documentation Found\n\n")
            
#### SHOW OPTION/DOCUMENTATION TYPES ####
    def __read_general_options__(self):
        
        option_path=self.pwd+"/AppOptions/Rosetta3-3"; #This will change once I everything else works.
        directorylist =os.listdir(option_path)
        for f in directorylist:
            if re.search(".txt", f) and (re.search("--", f)) and (not re.search("\._", f)) and (not re.search("~", f)) :
                FILE = open(option_path+"/"+f, 'r')
                fileSP = f.split("."); app = fileSP[0]
                self.appOPTIONS[app]=dict()
                for line in FILE:
                    if line == "\n":
                        break
                    lineSP = line.split("==")
                    self.appOPTIONS[app][lineSP[0]]=lineSP[1]
                FILE.close()
    
    def __show_manually_currated__(self):
        """
        Loads and shows manually currated info if possible.
        """
        self.appOPTIONS = pickle.load(open(self.pwd+"/option_binaries/Rosetta3-3.p")); #APP:Descriptions
        self.appDOC = pickle.load(open(self.pwd+"/option_binaries/Rosetta3-3Apps.p")); #APP:Documentation
        for app in self.appDOC:
            self.appDOC[app]["AppName"]=app
        self.array_of_applications= self.read_applications_from_directory(self.application_directory.get()); #Populate array_of_applications
        self.__populate_applications__(self.array_of_applications)
    
    def __show_doxygen__(self):
        """
        Shows doxygen_documentation
        """
        #Note - Most AppNames now match, so running/saving should work for the most part.
        app_map = self.parser.get_name_map()
        self.appOPTIONS = dict()
        for app_name in app_map:
            self.appDOC[app_name]=dict()
            data = self.parser.get_app_data(app_name)
            
            #Get all sections for documentation
            for section in data.sections:
                self.appDOC[app_name][section]=self.parser.get_app_section(app_name, section)
            
            #Get and split option data
            options = self.parser.get_app_section(app_name, "Options")
            if not options:continue
            option_array = options.split("\n")
            
            #NOTE: This may need tweeking
            self.appOPTIONS[app_name]=dict()
            for option in option_array:
                optionSP = option.split()
                #len 2 has no description, while length 3 does
                stripped = option.strip()
                if not stripped:continue
                if optionSP[0]=="@li":
                    try:
                        op = optionSP[1]
                        #Only want real options.  Not words on a new line.
                        if not op[0]=="-":
                            continue
                    except IndexError:
                        continue
                    
                    try:
                        desc = " ".join(optionSP[2:])
                        if self.appOPTIONS[app_name].has_key(op):
                            desc = desc+self.appOPTIONS[app_name][op]
                        self.appOPTIONS[app_name][op]=desc
                    except IndexError:
                        self.appOPTIONS[app_name][op]=""
                else:
                    #We don't get defaults in the option menu, which blows, but this looks like the only way - as some don't have defaults!
                    #Here, we should parse options_rosetta to get a description and default. That will be for later. 
                    #Some options list arn't ACTUAL c++ rosetta options.  So, it could get pretty nasty.
                    try:
                        op = optionSP[0]
                        #Only want real options.  Not words on a new line.
                        if not op[0]=="-":
                            continue
                        desc = " ".join(optionSP[1:])
                        if self.appOPTIONS[app_name].has_key(op):
                            continue
                            #desc = desc+self.appOPTIONS[app_name][op]
                        self.appOPTIONS[app_name][op]=desc
                    except IndexError:
                        pass
                        
        #self.array_of_applications = sorted(self.appOPTIONS); #ONLY App Documentation with options
        self.array_of_applications = sorted(self.appDOC); #All App Documentation
        self.__populate_applications__(self.array_of_applications)
        self.__read_general_options__()
        self.app_help_options = self.appOPTIONS ; #This is so that we do not have to reload.
        
    def __show_app_help_options__(self):
        """
        uses app -help to parse and print all info.
        Repopulates the specifOPTIONS dictionary.
        """
        #Note, I have heard popopen fails in 3.0 and less. So -
        
        print "Reading all available options for each application.  This may take a few minutes..."
        #Read applications, parse each app into dictionary.
        self.array_of_applications=self.read_applications_from_directory(self.application_directory.get())
        #try:
           # self.app_help_options
           # self.__populate_option_menu__(self.appDOC[self.last_app_clicked.get()]["AppName"])
            #return
        #except AttributeError:
            #print "Could not identify app"
        for a in self.appDOC:
            app = self.appDOC[a]["AppName"]
            app_path = self.application_directory.get()+"/"+app+'.'+self.appRoot.get()
            if os.path.exists(app_path):
                os.system(app_path+" -help > temp_options.txt")
                OPTIONS = open("temp_options.txt", 'r')
                self.appOPTIONS[a] = dict(); #Reset the option dictionary.
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
                    self.appOPTIONS[a][opt]=desc
                OPTIONS.close()
                os.system('rm temp_options.txt')
        self.__populate_option_menu__(self.last_app_clicked.get())
        self.app_help_options = self.appOPTIONS ; #This is so that we do not have to reload.
    
    
    #### FUNCTIONS ####

    def read_applications_from_directory(self, Apppath):
        """
        Reads all applications from the directory.
        Seperates .linuxgcc.
        -appends the appLIST, sets appRoot
        -Adds the basic options in globalOPTIONS to appLIST
        -sorts the appLIST.
        **DOES NOT UPDATE LISTBOX, etc.**
        """
        self.app_binaries=[]
        array_of_applications = []; #Reset the array.
	#If there is nothing in the directory, we seqfault without error...
        if os.path.exists(Apppath):
            list = os.listdir(Apppath)
        else:
            os.system('rm '+self.pwd+"/settings/PATHSETTINGS.txt")
            self.__load_rosetta_settings__
            list = os.listdir(Apppath)
        for apps in list:
            appsSP = apps.split(".")
            if appsSP[1]=="default":
                array_of_applications.append(appsSP[0])
                self.appRoot.set(appsSP[2])
                self.app_binaries.append(appsSP[0])
        for app in self.basic_OPTIONS:
            array_of_applications.append(app)
        array_of_applications.sort()
        
        return array_of_applications
    


    def __load_rosetta_settings__(self):
        """
        Load Rosetta settings text file.
        """
        
        setup_window = RosettaPathSetup.SetupRosettaPaths(self.main) 
        if not setup_window.result:
            #print "Please check paths and try again..."
            self.result = False
            self.main.destroy()
            
        else:
            self.result = True
        self.database_directory = StringVar()
        self.application_directory = StringVar()
        self.source_directory = StringVar()
        
        FILE = open(self.pwd+"/settings/PATHSETTINGS.txt", 'r')
        for line in FILE:
            lineSP = line.split()
            if lineSP[0]=="DATABASE":
                self.database_directory.set(lineSP[1])
            elif lineSP[0]=="APPLICATIONS":
                self.application_directory.set(lineSP[1])
            elif lineSP[0]=="ROSETTA_SOURCE":
                self.source_directory.set(lineSP[1])
        FILE.close()
            
            
           
    

        
#### Other Windows ####

    def setup_pathbuilder_objects(self):
        """
        Initialized pathbuilder Tk objects.
        """
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
        if not root:return ""
        self.defaultdir=root
        return root
    
    #These do not correspond to where the position is, but they should.
    def insertDIR_File(self):
        self.textbox_cmd_options.insert("1.end", " "+self.root.get()+"/"+self.filename.get()+" ")
        
    def insertDIR_(self):
        self.textbox_cmd_options.insert("1.end", " "+self.root.get())
        
    def insertSearch(self):
        filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
        if not filepath:return
        self.defaultdir = os.path.dirname(filepath)
        
        self.textbox_cmd_options.insert("1.end", " "+filepath)
            
    def shoClusterSetup(self):
        WinClusSetup = Toplevel(self.main)
        WinClusSetup.title("Qsub Cluster Setup")
        Clus = QsubClusterSetup.QsubClusterSetup()
        Clus.shoWindow(WinClusSetup, 0, 0)
        
        
#### OTHER FUNCTIONS ####

    def makeList(self):
        """
        Makes a list of PDB's from a directory. Use PyRosetta Toolkit for that.
        """
        
        directory = tkFileDialog.askdirectory(initialdir = self.defaultdir)
        if not directory:return
        self.defaultdir=directory
        
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
    
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        """
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP

if __name__ == '__main__':
    MainWindow = Tk()
    MainWindow.title("Window to Rosetta")
    SetupWindow = RosettaFlagFileBuilder(MainWindow)
    if not SetupWindow.result:
        sys.exit()
    SetupWindow.setTk("Paths not found or canceled")
    SetupWindow.shoTk(0, 0)
    SetupWindow.setMenu(MainWindow)
    MainWindow.mainloop()
