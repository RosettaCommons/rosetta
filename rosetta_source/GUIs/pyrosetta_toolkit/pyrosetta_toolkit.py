#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/pyrosetta_toolkit.py
## @brief  Main window for the toolkit.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import sys
from rosetta import *
from Tkinter import *
from os import listdir
from os import getcwd
from shutil import copy
from os import remove
from os import environ
from os import path
from os import system
import glob
import signal
import tkFileDialog
import tkMessageBox
import tkSimpleDialog

def location():
   '''
   Allows the script to be self-aware of it's path.
   So that it can be imported/ran from anywhere.
   '''
        
   p = os.path.abspath(__file__)
   pathSP = os.path.split(p)
   return pathSP


sys.path.append(location()[0]); #Allows all Window_Modules to use Modules.

from modules import tools

#from Modules import *
import functools
#import tools.pdbs

#import Window_Modules
#from Window_Modules.ScoreFunction.ScoreFxnControl import *
from window_modules import *
from window_main.input import *
from window_main.menu import *
from window_main.output import *
from window_main.quick_min import *
from window_main.simple_analysis import *






'''
print CDR_Analysis.__file__
#Um...
opts = ['app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE']), '-ex1', '-ex2aro', '-dun08', '-ramaneighbors true']
args = rosetta.utility.vector1_string()
args.extend(opts)
rosetta.core.init(args)
'''

class main_window:
   def __init__(self):
      '''
      Initializes the main window.
      Sets common global variables.
      '''
      self._tk_ = Tk() #Used to be c!
      self.pose = Pose()
      self.native_pose = Pose()
      self.pwd = self.location()[0]
      self.DesignDic = dict()
      
      
   ### Init ###
      self._init_objects()
      self._init_input_variables()
      self._init_output_variables()
      self._init_protocol_variables()
      self._init_menu_variables()
      
   ### Photo ###
      
      DesignPhoto =PhotoImage(file = (self.pwd+ "/media/DesignSmall.gif"))
      self.Photo = Label(master=self._tk_, image=DesignPhoto, relief = SUNKEN)
      self.Photo.image = DesignPhoto
      
   
   def run(self):
      self._tk_.title("PyRosetta Toolkit")
      self.show_gui()
      self._tk_.mainloop()
      
   def _init_objects(self):
      self.ScoreObject = ScoreFxn(); #Main Score Function Object. Holds Score.  Controls switching scorefunctions, etc.
      self.PyMOLObject = AdvancedPyMOL(self.pose); #PyMOL Object for advanced visualization.
      self.FullControlObject = FullControl(self.ScoreObject, self.pose); #Handles full control of protein.  This way, everything is saved...which is sorta cool.
   
   def _init_input_variables(self):
      self.input_class = Input(self._tk_, self)
      self.input_class.current_directory.set(self.pwd)
      self.current_directory= self.input_class.current_directory
      self.filename = self.input_class.filename
      self.PDBLIST = self.input_class.PDBLIST
      self.loops_as_strings = self.input_class.loops_as_strings
      
      
   def _init_output_variables(self):
      self.output_class = Outputs(self._tk_, self)
      self.outdir = self.output_class.outdir
      self.outname = self.output_class.outname
      
      
   def _init_protocol_variables(self):
      self.protocol_class = Minimization(self._tk_, self)
      self.simple_analysis_class = simple_analysis(self._tk_, self)
      
   def _init_menu_variables(self):
      self.menu_class = Menus(self._tk_, self)
      
   def location(self):
      '''
      Allows the script to be self-aware of it's path.
      So that it can be imported/ran from anywhere.
      '''
        
      p = os.path.abspath(__file__)
      pathSP = os.path.split(p)
      return pathSP
   
   def show_gui(self):
      '''
      Shows each piece of the main GUI.
      Does not do anything with the Window Modules, just each individual main component of the main window.
      NOTE: The main window was seperated as it was 1500 lines and hard to edit.
      
      Window Modules should be initialized through the Menus class in /main_window/menu.py
      '''
      
      self.Photo.grid(row = 12, column = 2, rowspan=12)
   ## Input ##
      self.input_class.setTk();   self.input_class.shoTk()
   ## Output ##
      self.output_class.setTk();  self.output_class.shoTk()
   ## Menu ##
      self.menu_class.setTk();    self.menu_class.shoTk()
   ## quick_min ##
      self.protocol_class.setTk();self.protocol_class.shoTk(20,3)
   ## Simple Analysis ##
      self.simple_analysis_class.setTk(); self.simple_analysis_class.shoTk(11,2)




if __name__ == '__main__':
   rosetta.init()
   main_window_class = main_window()
   main_window_class.run()
