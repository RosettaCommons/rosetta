#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/parse_doxygen.py
## @brief  Functions to parse doxygen app documention for the GUI
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import glob
import re

#Tkinter Imports
from Tkinter import *

class DoxygenParser:
    """
    This class is specifically for parsing Rosetta app doxygen documentation.
    """
    def __init__(self, rosetta_source_directory):
        self.doxygen_directory = rosetta_source_directory+"/doc/apps/public"
        self.app_names = dict() #[string app name]:[path]
        self.all_app_data = []; #[Array of Data classes]
        self.name_map = dict(); #[string app_name]:[Data class]
        self.__get_info__(); #Main setup method
    
    

    def get_app_data(self, app_name):
        """
        Returns a data object.  Which just has a bunch of StringVars that can be accessed through their variable name or by dictionary matching a string.
        """
        
        return self.name_map[app_name]
    
    def get_app_section(self, app_name, section):
        """
        Returns section string for a particular app
        """
        return self.name_map[app_name].get_section(section)
    
    
    
    def get_name_map(self):
        """
        Returns name map dictionary.
        """
        return self.name_map
    
    def get_all_app_data(self):
        """
        Returns array of all Data objects found.
        """
        return self.all_app_data
    
    def __get_info__(self):
        """
        Parses all app doxygen documentation.
        """
        files = glob.glob(self.doxygen_directory+"/*.dox")
        for path in sorted(files):
            
            #Skip GUI doxygens
            if re.search("pyrosetta_toolkit", path) or re.search("RosettaFlagFileBuilder", path):
                continue
            d = Data(path)
            self.all_app_data.append(d)
            self.name_map[d.app_name]=d
            
        #Grab for carbohydrates.  Maybe all of them will be organized soon...
        self.app_types = os.walk(self.doxygen_directory).next()[1]
        for type in self.app_types:
            dir_path = self.doxygen_directory+"/"+type
            files = glob.glob(dir_path+"/*.dox")
            for path in sorted(files):
                if re.search("pyrosetta_toolkit", path) or re.search("RosettaFlagFileBuilder", path):
                    continue
                d = Data(path)
                self.all_app_data.append(d)
                self.name_map[d.app_name]=d
       

class Data:
    """
    Create this for each doxygen file you want info on.  Container and parser.  Holds all variables as stringvars from doxygen file.
    Also has info in a dictionary, which you can access through it's get method
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.app_name = os.path.basename(filepath).split(".")[0]
        #Variables.  Will all be long strings.  All we need to do is put them in the GUI.
        #First we initialize them
        self.metadata=StringVar(); self.code_demos = StringVar(); self.refs = StringVar()
        self.purpose =StringVar(); self.algorithm =  StringVar(); self.inputs=StringVar()
        self.options=StringVar(); self.tips = StringVar(); self.outputs = StringVar()
        self.postprocess = StringVar(); self.new_stuff = StringVar();
        self.AppName = StringVar(); #The real app name from the file.
        self.limits = StringVar(); self.modes = StringVar()
        #Now we have a dictionary for matching them.  This would be simple in perl.
        #Since for the lower case, people don't seem to be following the convention (fp_algorithm?), we use the cap form and match the line.
        self.sections = {
            "Metadata":self.metadata,
            "Code and Demo":self.code_demos,
            "References":self.refs,
            "Purpose":self.purpose,
            "Algorithm":self.algorithm,
            "Input Files":self.inputs,
            "Options":self.options,
            "Tips":self.tips,
            "Expected Outputs":self.outputs,
            "Post Processing":self.postprocess,
            "new_stuff":self.new_stuff,
            "AppName":self.AppName,
            "Limitations":self.limits,
            "Modes":self.modes
        }
       
        self.__parse_sections__()
    
    
    def get_section(self, section):
        """
        Access to it's dictionary.  Return specific string for the section you are interested in.
        *This is how you should be interacting with this class.
        """
        return self.sections[section].get()
        
    def __parse_sections__(self):
        """
        Runs parse section for each section
        """
        for section in self.sections:
            self.__parse_section__(section)
    
    def __parse_section__(self, section):
        """
        Parse an individual section, and set variable in self.sections.  Yes this takes longer to read every file once for every section.  But, a few seconds here doesn't matter.
        Keeps original "\n" charactors within section.  Formatting may be aweful.  Deals with options differently.
        Code is not easily readable.  Niether is the doc!
        """
        #print "Finding "+section
        FILE = open(self.filepath, 'r')
        while True:
            #Don't want empty lines being split.  So we continue if one is found.  Will have empty lines in the section once @section identifier is found
            line = FILE.readline()
            if line=="\n":
                continue
            if not line:break
            lineSP = line.split()
            if len(lineSP)<=1:
                continue
            if lineSP[0]=='@section' and re.search(section, line):
                #Grab info on next lines until next section is found.  Then break.
                section_string = ""
                new_section = False

                while not new_section:
                    line_string = FILE.readline()
                    if not line_string:break
                    line_stringSP = line_string.split()
                    if len(line_stringSP)<=1:
                        if not re.search('@verbatim', line_string) and not re.search('@endverbatim', line_string) and section=='Options':
                            continue
                        
                    #if line_stringSP[0]=="@section": new_section=True
                    if re.search('@section', line_string): new_section=True
                    else:
                        #We only want verbatim options or LI.  We then parse it in protocol builder.  Which is a bitch because everyone follows something different!
                        if section=='Options':
                            verbatim=False
                            #print line_stringSP[0]
                            if line_stringSP[0]=='@li':
                                section_string = section_string+line_string
                                continue
                            elif re.search('@verbatim', line_string):
                                verbatim=True
                            while verbatim:
                                line_string = FILE.readline()
                                if re.search('@endverbatim', line_string):
                                    verbatim=False
                                    continue
                                if re.search('@section', line_string):
                                    new_section=True
                                    verbatim=False
                                    break
                                section_string = section_string+line_string
                        else:
                            #So that verbatims are not in the result string
                            if not re.search("verbatim", line_string):
                                line_string = line_string.replace('@li', "")
                                section_string = section_string+line_string

                    #else:new_section=False
                #Next section found, returning.
                FILE.close()
                #print "Section: "+section
                #print section_string
                self.sections[section].set(section_string)
                return
            
            elif lineSP[0]=='@page' and section=="AppName":
                self.sections["AppName"].set(lineSP[1])
                return
                

if __name__ == '__main__':
    """
    Testing
    """
    Tk()
    source_dir = ""
    parser = DoxygenParser(source_dir)
    op = parser.get_app_section("FloppyTail", "Options")
    #Test not matching:
    print op
    """
    data = parser.get_all_app_data()
    list = os.listdir(source_dir+"/bin")
    apps=[]
    found = []
    not_found = []
    
    for app in list:
        appSP = app.split(".")
        if appSP[1]=="default":
            apps.append(appSP[0])
    for d in data:
        found_bool = False
        for app in apps:
            if re.search(d.get_section("AppName"), app):
                found_bool=True
                break
        if found_bool:
            found.append(d.app_name)
        else:
            not_found.append(d.app_name)
    
    for app in sorted(found):
        print "Found: "+app
    for app in sorted(not_found):
        print "Not Found: "+app
    """    
        
#@section metadata Metadata

#@section code_demos Code and Demo

#@section refs References

#@section purpose Purpose

#@section algorithm Algorithm

#@limits Limitations

#@modes Modes

#@section inputs Input Files

#@section ft_options Options

#@section tips Tips

#@section outputs Expected Outputs

#@section postprocess Post Processing

#@section new_stuff New things since last release

        
        
        
        
