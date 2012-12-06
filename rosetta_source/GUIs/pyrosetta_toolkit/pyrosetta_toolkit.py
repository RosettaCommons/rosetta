#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/pyrosetta_toolkit.py
## @brief  Main window for the toolkit.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#Rosetta Imports
from rosetta import *
from modules import tools
import functools

#Python Imports
from os import listdir
from os import getcwd
from shutil import copy
from os import remove
from os import environ
from os import path
from os import system
import glob
import signal
import sys

#Tk Imports
from Tkinter import *
from Tkinter import Frame as FrameTk
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
import tkFont

#Toolkit Imports
from window_main import global_variables
from window_main.menu import *
from window_main.frames.InputFrame import InputFrame
from window_main.frames.OutputFrame import OutputFrame
from window_main.frames.QuickProtocolsFrame import QuickProtocolsFrame
from window_main.frames.SimpleAnalysisFrame import SimpleAnalysisFrame
from window_main.IO.GUIInput import GUIInput
from window_main.IO.GUIOutput import GUIOutput

from window_modules.pymol_integration.PyMOL import AdvancedPyMOL
from window_modules.scorefunction.ScoreFxnControl import ScoreFxn
from window_modules.full_control.FullControlWindow import FullControlWindow


def location():
   """
   Allows the script to be self-aware of it's path.
   So that it can be imported/ran from anywhere.
   """
        
   p = os.path.abspath(__file__)
   pathSP = os.path.split(p)
   return pathSP
sys.path.append(location()[0]); #Allows all Window_Modules to use Modules.


class main_window:
   def __init__(self):
      """
      Initializes the main window.
      Sets common global variables.
      """
      self._tk_ = Tk()
      self.pose = Pose()
      self.native_pose = Pose()
      self.current_directory = global_variables.current_directory = self.location()[0]
      self.DesignDic = dict()
      
      
   ### Init ###
      self._initialize_GUI()
      self._initialize_Frames()

   ### TextBox ###
      
      self.textbox_frame = FrameTk(self._tk_, bd=3, relief=GROOVE)
      outfont = tkFont.Font(family="Helvetica", size=11)
      self.output_textbox= Text(self.textbox_frame,wrap="word", height=7,width=103,font = outfont)
      
      self.old_stdout = sys.stdout
      
      self.output_class.terminal_output.trace_variable('w', self.output_tracer)
      self.output_class.terminal_output.set(0)
      print "For additional protocol options, please use the Option System Manager.\n"
      print "Please see RosettaCommons for full documentation and references for all protocols and tools utilized in the GUI\n"
      
   def _initialize_GUI(self):
      """
      Creates object for the GUI
      """
      self.options_class = OptionSystemManager(global_variables.current_directory)
      self.input_class = GUIInput(self)
      self.output_class = GUIOutput(self)
      
      self.score_class = ScoreFxn(); #Main Score Function Object. Holds Score.  Controls switching scorefunctions, etc.
      self.pymol_class = AdvancedPyMOL(self.pose); #PyMOL Object for advanced visualization.
      self.fullcontrol_class = FullControlWindow(self.score_class, self.pose, self.input_class, self.output_class); #Handles full control of protein.  This way, everything is saved...which is sorta cool.
      

      
   def _initialize_Frames(self):
      """
      Creates the Frame Objects that will go in the main window
      """
      self.input_frame = InputFrame(self._tk_, self, self.input_class, bd=1, relief=SUNKEN)
      self.output_frame = OutputFrame(self._tk_, self, self.output_class, bd=1, relief = SUNKEN)
      
      self.protocol_frame = QuickProtocolsFrame(self._tk_, self, self.output_class, bd=1, relief=SUNKEN)
      self.simple_analysis_frame = SimpleAnalysisFrame(self._tk_, self, bd=2, relief=SUNKEN)
      self.menu_class = Menus(self._tk_, self)
      

   
   def show_gui(self):
      """
      Shows each piece of the main GUI.
      Does not do anything with the Window Modules, just each individual main component of the main window.
      These Inhereit from the Frame class.  See one of these for an example.
      Window Modules should be initialized through the Menus class in /window_main/menu.py
      """
      #6x4 Grid Pain in the ass.  At some point, everything will move to Qt - Either in Python or C++
      
      #Grid:
      self.menu_class.setTk();    self.menu_class.shoTk()
      
      self.input_frame.grid(row=1, column=0, rowspan=7, padx=15, pady=15);
      self.output_frame.grid(row=0, column=1, rowspan=2, pady=3);
      
      self.protocol_frame.grid(row=3, column=1, rowspan=4, padx=5)
      self.simple_analysis_frame.grid(row=0, column=0, padx=5, pady=5)
      
      
      ### Text Output ###
      self.output_textbox.grid(column=0, row = 9, rowspan=2, columnspan=3,sticky=W+E)
      self.textbox_frame.grid(column=0, row=9, rowspan=2,  columnspan=3, sticky=W+E, pady=3, padx=6)
      #self.Photo.grid(row = 0, column = 2, rowspan=4)
      
      """
      #Pack:
      self.menu_class.setTk();    self.menu_class.shoTk()
      self.input_class.pack(side=LEFT, padx=3, pady=3)
      self.output_class.pack(padx=3, pady=3)
      
      self.simple_analysis_frame.pack(padx=3, pady=3)
      self.protocol_frame.pack(padx=3, pady=3)
      #self.output_textbox.pack(side=BOTTOM, padx=3, pady=3)
      #self.output_frame.pack(side=BOTTOM, padx=3, pady=3)
      """
      
   def run(self):
      self._tk_.title("PyRosetta Toolkit")
      self.show_gui()
      self._tk_.mainloop()
      
   def redirect_stdout_to_textbox(self):
      sys.stdout = self; #Set stdout to be redirected to textbox using the write function override.
      print "Redirect stdout to textbox"
      
   def redirect_stdout_to_default(self):
      sys.stdout = self.old_stdout
      print "Redirect stdout to default"
      
   def write(self, text):
      self.output_textbox.insert(1.0, text)

   def output_tracer(self, name, index, mode):
      """
      Controls where stdout goes.  Textbox or Terminal.
      Does not override Tracer for now.
      """
      
      varvalue = self.output_class.terminal_output.get()
      
      if (varvalue):
         self.redirect_stdout_to_default()
      else:
         self.redirect_stdout_to_textbox()
         
   def location(self):
      """
      Allows the script to be self-aware of it's path.
      So that it can be imported/ran from anywhere.
      """
        
      p = os.path.abspath(__file__)
      pathSP = os.path.split(p)
      return pathSP
   
class MainTracer(rosetta.basic.PyTracer):
   def __init__(self, textbox):
        rosetta.basic.PyTracer.__init__(self)
        self.textbox = textbox
   def output_callback(self, s):
      s = " "+s
      self.textbox.insert(1.0, s)
      #print s
   
      



if __name__ == '__main__':
   rosetta.init()

   rosetta.init()
   main_window_class = main_window()
   TR = MainTracer(main_window_class.output_textbox)
   rosetta.basic.Tracer.set_ios_hook(TR, rosetta.basic.Tracer.get_AllChannels_string(), False)
   #rosetta.init(extra_options="-mute all")
   
   main_window_class.run()
