#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/QsubClusterSetup.py
## @brief  Window for running rosetta executables on a cluster using QSUB.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import pickle
import os
import tkFont
import re
import math
from settings import QsubSettings

class QsubClusterSetup():
    """
    This class is used to make it easier to run Jobs on a cluster using qsub.  Use the functions in pyrosetta toolkit to combine resultant directories and rename.
    If you know your application is JD2, use MPIClusterSetup.
    Can run by itself. Use SSH -X to set X11 display forwarding in order to use the window on the cluster.
    """
    
    def __init__(self):
        self.pwd = (self.location())[0]
        self.jobs = StringVar()
        self.jobs.set("0")
        self.struc_per_job= StringVar()
        self.struc_per_job.set("0")  

        self.queue = StringVar()
        
        self.jobname = StringVar()
        self.jobtotal = StringVar()
        self.jobtotal.set("0")
        
        self.settings = QsubSettings.QsubSettings()
        self.QUEUEMETHOD = self.settings.QUEUE_LIST; #Various Queues for the cluster to use.
        if self.QUEUEMETHOD:
            self.queue.set(self.QUEUEMETHOD[0])
    def shoWindow(self, main, row, column):
        """
        Sets Tk variables and shows the main window.
        """
        
        self.main = main
        self.row = row; self.column = column
        
        self.job_label = Label(self.main, text="# of Jobs")
        self.struc_per_job_label = Label(self.main, text = "# of Structures/Job")
        self.total_struc_label= Label(self.main, textvariable = self.jobtotal)
        self.job_name_label = Label(self.main, text="Job name")
        self.queue_label = Label(self.main, text = "Queue")
        self.job_entry = Entry(self.main, textvariable = self.jobs, justify=CENTER); 
        self.struc_per_job_entry = Entry(self.main, textvariable = self.struc_per_job, justify=CENTER);
        self.job_name_entry = Entry(self.main, textvariable = self.jobname, justify=CENTER);
        self.queueOptions = OptionMenu(self.main, self.queue, *self.QUEUEMETHOD)
        self.specQ = Button(self.main, text="Check Q Avalability:", command = lambda: self.checkspecificQ()).grid(row=self.row+2, column=self.column+3, sticky=W+E)
        self.allQ = Button(self.main, text="Show all Q", command = lambda: self.checkallQ()).grid(row=self.row+3, column=self.column+3, sticky=W+E)
        self.runbutton_ = Button(self.main, text = "Run Config", command = lambda: self.runScript())
        self.calcbutton_ = Button(self.main, text = "calc", command = lambda: self.calcTotal()); #Couldn't get the call back to work...

        
        self.qstat = Button(text = "qstat", command = lambda: self.checkQ())
        
        #Show all the items.
        self.job_label.grid(row = self.row, column=self.column); self.job_entry.grid(column = self.column+1, row = self.row)
        self.struc_per_job_label.grid(row = self.row+1, column=self.column); self.struc_per_job_entry.grid(column=self.column+1, row=self.row+1)
        self.total_struc_label.grid(row = self.row+2, column = self.column+1)
        self.job_name_entry.grid(row = self.row, column = self.column+2)
        self.job_name_label.grid(row = self.row, column = self.column+3)
        self.queueOptions.grid(row= self.row+1, column=self.column+2, sticky = W+E)
        self.queue_label.grid(row=self.row+1, column=self.column+3)
        self.runbutton_.grid(row = self.row+5, column = self.column, columnspan = 3, sticky=W+E)
        self.calcbutton_.grid(row = self.row+2, column=self.column)
        
        self.frameHelp = Frame(self.main, bd=3, relief=GROOVE)
        self.textHelp = Text(self.frameHelp,wrap="word", height = 3, background = 'white')
        self.frameHelp.grid(row = self.row+8, column = self.column, columnspan = 4, sticky = W+E, pady = 3)
        self.textHelp.grid(row = self.row+8, column = self.column, columnspan = 4)
        self.qstat.grid(row = self.row+5, column=self.column+3)
        
        self.shoMenu()
    def shoMenu(self):
        """
        Sets up the GUI's menu, shows it.
        """
        
        self.MenBar = Menu(self.main)
        self.calc = Menu(self.main, tearoff=0)
        self.calc.add_command(label = "Calculate Runtime", command = lambda: self.calcRunTime())
        self.MenBar.add_cascade(label = "Caclulate", menu=self.calc)
        self.appMen = Menu(self.main, tearoff=0)
        #self.appMen.add_command(label = "App Specific Options")
        self.appMen.add_command(label = "Setup Cluster Settings", command = lambda: self.shoSetup())
        self.appMen.add_separator()
        self.appMen.add_command(label = "Repeat for List of Proteins", command = lambda: self.runScript(True))
        self.MenBar.add_cascade(label = "Options", menu = self.appMen)
        self.main.config(menu=self.MenBar)
        
    def shoSetup(self):
        """
        Shows the QsubSettings window to setup paths and queues for qsub.
        """
        
        win = Toplevel(self.main)
        self.settings.setTk(win)
        self.settings.shoTk(win, 0, 0)
    
    def checkspecificQ(self):
        """
        checks the specific Queue for availability
        """
        #self.textHelp.delete(1.0, END)
        #self.textHelp.insert(1.0, os.system('qstat -q '+self.queue.get()))
        os.system('qstat -q '+self.queue.get())
        
    def checkQ(self):
        os.system('qstat')
        return
    
    def checkallQ(self):
        """
        checks all Q available...
        """
        os.system(self.settings.maui_showq_path.get())
        
    def runScript(self, repeat=False):
        """
        Asks user for paths to config and optionally a PDBLIST
        Runs qsub_script.py by passing commandline arguments.
        """
        
        config = tkFileDialog.askopenfilename(title="config", initialdir = self.pwd)
        
        if repeat:
            filename = tkFileDialog.askopenfilename(title="PDBList", initialdir = self.pwd)
            if not filename:
                return
            
        if config:
            args = []
            args.append("--config "+config)
            args.append("--jobs "+self.jobs.get())
            args.append("--stru "+self.struc_per_job.get())
            args.append("--qsub "+self.settings.qsub_path.get())
            args.append("--queue "+self.queue.get())
            if self.jobname.get():
                args.append("--jobname "+self.jobname.get())
            #args.append("--debug")
            
            if self.settings.qsub_temp.get():
                args.append("--qsubtemp "+self.settings.qsub_temp.get())
            if self.settings.tempscripts.get():
                args.append("--tempscripts "+self.settings.tempscripts.get())
            if repeat:
                args.append("--pdb_list "+filename)
                
            arg_string = " ".join(args)
            command= 'python '+ self.pwd+"/qsub_script.py "+arg_string
            
            print "Running command: "
            print command
            os.system(command)
            self.textHelp.delete(1.0, END)
            self.textHelp.insert(1.0, "Cluster run started...")
        return

    def calcRunTime(self):
        """
        Calculates approximate runtime, graphs it if numpy and matplotlib can be imported.
        """
        
        time =  tkSimpleDialog.askfloat(title = "Time", prompt = "Approx Time(min) per Pose: ")
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            #totalNodes=tkSimpleDialog.askinteger(title = "Nodes", prompt = "Total Nodes on cluster")
            totalNodes = 100
            x_axes = []; y_axes = []
            m = ((time/60)*self.calcTotal())
            for x in range(1, totalNodes+1):
                y = (((time/60)*self.calcTotal())/x)/24
                x_axes.append(x),y_axes.append(y)
            plt.xlabel('Nodes')
            plt.ylabel('Days')
            plt.plot(x_axes, y_axes, 'o')
            plt.grid(True)
            plt.show()
            return
        
        except ImportError:
            self.textHelp.delete(1.0, END)
            nodes = tkSimpleDialog.askinteger(title = "Nodes", prompt="Approx Number of Nodes: ")
        
            time = time/60
            TotalCpuTime = self.calcTotal() * time
            TotalCpuTimeDays = TotalCpuTime/24
            self.textHelp.insert(1.0, "Total CPU time is: "+repr(TotalCpuTime)+" Hours or "+repr(TotalCpuTimeDays)+" Days")
            TotalTime = TotalCpuTime/nodes
            TotalTimeDays = TotalTime/24
            self.textHelp.insert(1.0, "Total Time is: "+repr(TotalTime)+ " Hours or "+ repr(TotalTimeDays)+ " Days"+"\n")
            
    def calcTotal(self):
        try:
            jobs = int(self.jobs.get())
            stru = int(self.struc_per_job.get())
            total = str(jobs*stru)
            self.jobtotal.set(total)
            return int(total)
        except TypeError:
            return 0
        except UnboundLocalError:
            return 0
            
        
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
    MainWindow.title("Qsub Cluster Setup")
    SetupWindow = QsubClusterSetup()
    SetupWindow.shoWindow(MainWindow, 0, 0)
    MainWindow.mainloop()
