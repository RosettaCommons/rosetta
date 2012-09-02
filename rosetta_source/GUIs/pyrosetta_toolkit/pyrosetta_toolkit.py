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

import tkFont
from rosetta import *

from Tkinter import *
from Tkinter import Frame as FrameTk
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
   """
   Allows the script to be self-aware of it's path.
   So that it can be imported/ran from anywhere.
   """
        
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
from window_main.InputFrame import *
from window_main.menu import *
from window_main.OutputFrame import *
from window_main.QuickProtocolsFrame import *
from window_main.SimpleAnalysisFrame import *





"""
print CDR_Analysis.__file__
#Um...
opts = ['app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE']), '-ex1', '-ex2aro', '-dun08', '-ramaneighbors true']
args = rosetta.utility.vector1_string()
args.extend(opts)
rosetta.core.init(args)
"""

class main_window:
   def __init__(self):
      """
      Initializes the main window.
      Sets common global variables.
      """
      self._tk_ = Tk()
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

   ### TextBox ###
      
      self.output_frame = FrameTk(self._tk_, bd=3, relief=GROOVE)
      outfont = tkFont.Font(size=12)
      self.output_textbox= Text(self.output_frame,wrap="word", height=7,width=103,font = outfont)
      
      sys.stdout = self; #Set stdout to be redirected to textbox using the write function override.
      
      
   def write(self, text):
      self.output_textbox.insert(1.0, text)
      
   def run(self):
      self._tk_.title("PyRosetta Toolkit")
      self.show_gui()
      self._tk_.mainloop()
      
   def _init_objects(self):
      self.ScoreObject = ScoreFxn(); #Main Score Function Object. Holds Score.  Controls switching scorefunctions, etc.
      self.PyMOLObject = AdvancedPyMOL(self.pose); #PyMOL Object for advanced visualization.
      self.FullControlObject = FullControl(self.ScoreObject, self.pose); #Handles full control of protein.  This way, everything is saved...which is sorta cool.
   
   def _init_input_variables(self):
      self.input_class = InputFrame(self._tk_, self)
      self.input_class.current_directory.set(self.pwd)
      self.current_directory= self.input_class.current_directory
      self.filename = self.input_class.filename
      self.PDBLIST = self.input_class.PDBLIST
      self.loops_as_strings = self.input_class.loops_as_strings
      
      
   def _init_output_variables(self):
      self.output_class = OutputFrame(self._tk_, self)
      self.outdir = self.output_class.outdir
      self.outname = self.output_class.outname
      
      
   def _init_protocol_variables(self):
      self.protocol_class = QuickProtocolsFrame(self._tk_, self, bd=2, relief=RAISED)
      self.simple_analysis_class = SimpleAnalysisFrame(self._tk_, self, bd=2, relief=RAISED)
      
   def _init_menu_variables(self):
      self.menu_class = Menus(self._tk_, self)
      
   def location(self):
      """
      Allows the script to be self-aware of it's path.
      So that it can be imported/ran from anywhere.
      """
        
      p = os.path.abspath(__file__)
      pathSP = os.path.split(p)
      return pathSP
   
   def show_gui(self):
      """
      Shows each piece of the main GUI.
      Does not do anything with the Window Modules, just each individual main component of the main window.
      These Inhereit from the Frame class.  See one of these for an example.
      Window Modules should be initialized through the Menus class in /main_window/menu.py
      """
      #6x4 Grid
      
      #Grid:
      self.menu_class.setTk();    self.menu_class.shoTk()
      self.input_class.grid(row=0, column=0, rowspan=6, pady=3);
      self.output_class.grid(row=0, column=1, rowspan=2, pady=3, sticky=W);
      self.simple_analysis_class.grid(row=2, column=1, sticky=W)
      self.protocol_class.grid(row=5, column=1, sticky=W, padx=5)
      
         ### Text Output ###
      self.output_textbox.grid(column=0, row = 7, rowspan=3, columnspan=3,sticky=W+E)
      self.output_frame.grid(column=0, row=7, rowspan=3,  columnspan=3, sticky=W+E, pady=3, padx=6)
      #self.Photo.grid(row = 0, column = 2, rowspan=4)
      
      """
      #Pack:
      self.menu_class.setTk();    self.menu_class.shoTk()
      self.input_class.pack(side=LEFT, padx=3, pady=3)
      self.output_class.pack(padx=3, pady=3)
      
      self.simple_analysis_class.pack(padx=3, pady=3)
      self.protocol_class.pack(padx=3, pady=3)
      #self.output_textbox.pack(side=BOTTOM, padx=3, pady=3)
      #self.output_frame.pack(side=BOTTOM, padx=3, pady=3)
      """
class MainTracer(rosetta.basic.PyTracer):
   def __init__(self, textbox):
        rosetta.basic.PyTracer.__init__(self)
        self.textbox=textbox
   def output_callback(self, s):
      s = "  "+s
      self.textbox.insert(1.0, s)
      #print s
   def setTextBox(textbox):
      self.textbox = textbox
   
      



if __name__ == '__main__':
   rosetta.init()

   rosetta.init()
   main_window_class = main_window()
   TR = MainTracer(main_window_class.output_textbox)
   rosetta.basic.Tracer.set_ios_hook(TR, rosetta.basic.Tracer.get_AllChannels_string())
   rosetta.init(extra_options="-mute all")
   
   main_window_class.run()
