#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/RosettaClusterSetup.py
## @brief  Window for running rosetta configs on cluster.  Needs work to be general.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import pickle
import os
import tkFont
import re
import math

class ClusterSetup():
    def __init__(self):
        self.pwd = (self.location())[0]
        self.jobs = StringVar()
        self.struc= StringVar()
        self.SEEDMETHOD = [
            'Job x Total Structures',
            'Job',
            'Custom'
        ]
        self.seed = StringVar()
        self.seed.set(self.SEEDMETHOD[0])
        self.QUEUEMETHOD = [
            'nmf',
            'apache',
            'batch',
            'glide',
            'dna',
            'q32'
        ]
        
        #Need to add a trace to get this to work correctly...
        self.queueproc = {
            'nmf':2,
            'apache':8,
            'batch':2,
            'glide':2,
            'dna':2,
            'q32':2
        }
        self.queue = StringVar()
        self.queue.set(self.QUEUEMETHOD[0])
        self.jobname = StringVar()
        self.jobtotal = StringVar()
        self.jobtotal.set("0")
        #self.proc = StringVar()
        #self.proc.set(self.queueproc['nmf'])
    def shoWindow(self, main, row, column):
        self.main = main
        self.row = row; self.column = column
        
        self.jlab = Label(self.main, text="# of Jobs")
        self.slab = Label(self.main, text = "# of Structures/Job")
        self.total= Label(self.main, textvariable = self.jobtotal)
        self.name = Label(self.main, text="Job name")
        #self.seedlabel = Label(self.main, text = "Name Method")
        self.queuelabel = Label(self.main, text = "Queue")
        #self.nodelabel = Label(self.main, text = 'proc/Node')
        #self.nodeentry = Entry(self.main, textvariable = self.proc, justify=CENTER)
        self.jent = Entry(self.main, textvariable = self.jobs, justify=CENTER); 
        self.sent = Entry(self.main, textvariable = self.struc, justify=CENTER);
        self.nent = Entry(self.main, textvariable = self.jobname, justify=CENTER);
        
        #self.seedOptions = OptionMenu(self.main, self.seed, *self.SEEDMETHOD)
        '''
        Seed uses a constant 1 million + 10 for each job being kicked off on different cpus.
        This should be an option in the future...
        -Also note that qsub directories should be an option in the future as well, as it is specific for the cluster.
        -Next, we want to add the ability to use JD2 options for those applications that require it.
        '''
        self.queueOptions = OptionMenu(self.main, self.queue, *self.QUEUEMETHOD)
        self.specQ = Button(self.main, text="Check Q Avalability:", command = lambda: self.checkspecificQ()).grid(row=self.row+4, column=self.column+1, sticky=W+E)
        self.allQ = Button(self.main, text="Show all Q", command = lambda: self.checkallQ()).grid(row=self.row+7, column=self.column, sticky=W+E)
        self.savebutton_ = Button(self.main, text = "Save Script", command = lambda: self.saveScript())
        self.runbutton_ = Button(self.main, text = "Run Config", command = lambda: self.runScript())
        #self.runCurbutton_ = Button(self.main, text = "Run Current Config", command = lambda: self.runCurrent())
        self.calcbutton_ = Button(self.main, text = "calc", command = lambda: self.calcTotal()); #Couldn't get the call back to work...
        #Small Menu.
        self.MenBar = Menu(self.main)
        self.calc = Menu(self.main, tearoff=0)
        self.calc.add_command(label = "Calculate Runtime", command = lambda: self.calcRunTime())
        self.MenBar.add_cascade(label = "Caclulate", menu=self.calc)
        self.appMen = Menu(self.main, tearoff=0)
        self.appMen.add_command(label = "App Specific Options")
        self.appMen.add_command(label = "Cluster Specific Options")
        self.appMen.add_separator()
        self.appMen.add_command(label = "Repeat for List of Proteins", command = lambda: self.listRepeat())
        self.MenBar.add_cascade(label = "Options", menu = self.appMen)
        self.main.config(menu=self.MenBar)
        
        self.qstat = Button(text = "check Q", command = lambda: self.checkQ())
        #Show all the items.
        self.jlab.grid(row = self.row, column=self.column); self.jent.grid(column = self.column+1, row = self.row)
        self.slab.grid(row = self.row+1, column=self.column); self.sent.grid(column=self.column+1, row=self.row+1)
        self.total.grid(row = self.row+2, column = self.column+1)
        self.nent.grid(row = self.row+2, column = self.column+2)
        self.name.grid(row = self.row+2, column = self.column+3)
        #self.seedOptions.grid(row = self.row+3, column=self.column+2, sticky=W+E)
        #self.seedlabel.grid(row = self.row+3, column=self.column+3)
        self.queueOptions.grid(row= self.row+4, column=self.column+2, sticky = W+E)
        self.queuelabel.grid(row=self.row+4, column=self.column+3)
        #self.nodelabel.grid(row=self.row+3, column=self.column+3)
        #self.nodeentry.grid(row=self.row+3, column=self.column+2)
        self.savebutton_.grid(row=self.row+7, column=self.column+1, columnspan = 2, sticky=W+E)
        self.runbutton_.grid(row = self.row+6, column = self.column+1, sticky=W+E)
        self.runCurbutton_.grid(row = self.row+6, column=self.column+2, sticky=W+E)
        self.calcbutton_.grid(row = self.row+2, column=self.column)
        
        self.frameHelp = Frame(self.main, bd=3, relief=GROOVE)
        self.textHelp = Text(self.frameHelp,wrap="word", height = 3, background = 'white')
        self.frameHelp.grid(row = self.row+8, column = self.column, columnspan = 4, sticky = W+E, pady = 3)
        self.textHelp.grid(row = self.row+8, column = self.column, columnspan = 4)
        self.qstat.grid(row = self.row+9, column=self.column+3)
    
    def listRepeat(self):
        '''
        Allows us to give a list, and repeat the config on a different processor for each protein in list.
        '''
        #First, we ask for the list.
        filename = tkFileDialog.askopenfilename(title="PDBList", initialdir = self.pwd)
        self.saveScript(F=1, repeat = filename)
        config = tkFileDialog.askopenfilename(title="Config", initialdir = self.pwd)
        command= 'python '+ self.pwd+"/tempScript.py "+config
        os.system(command)
        self.textHelp.delete(1.0, END)
        self.textHelp.insert(1.0, "Cluster run started...")
        return
    
    def checkspecificQ(self):
        '''
        checks the specific Queue for availability
        '''
        #self.textHelp.delete(1.0, END)
        #self.textHelp.insert(1.0, os.system('qstat -q '+self.queue.get()))
        os.system('qstat -q '+self.queue.get())
    def checkQ(self):
        os.system('qstat')
        return
    def checkallQ(self):
        '''
        checks all Q available...
        '''
        os.system('/usr/local/maui/bin/showq')
    def saveScript(self, F=0, repeat=0):
        '''
        Script needs to take in a rosetta protocol file.
        It needs to write a file for each job for the cluster, or at least kick it off.
        It needs take in the protocol file and concatonate it to what it was using.
        I could do this in perl, but I will try to do it in python.
        '''
        if repeat!=0:
            List = repeat
            repeat = True
        if F==0:
            F = tkFileDialog.asksaveasfilename(initialdir = self.pwd)+".py"
        else:
            F = self.pwd + "/tempScript.py"
        FILE = open(F, 'w')
        FILE.write("import os\n"+"import sys\n"+"import re\n")
        #Here, we give the script some variables.
        variables = dict()
        variables={#Surprisingly, this actually works quite nicely.
            'jobs':self.jobs.get(),
            'stru':self.struc.get(),
            'jobname':self.jobname.get(),
            'queue':self.queue.get(),
        }
        
        

        
        if repeat:
            variables['List']=List
        for key in variables:
            FILE.write(key+"="+"\""+variables[key]+"\""+"\n")
        #Fixes integers.
        FILE.write("\n\n"+"jobs=int(jobs)\n"+"stru=int(stru)\n\n")
        
        #Makes directories for the jobname (tempscripts and qsubtemp need options)
        FILE.write("if not os.path.exists('/common/madsci/Modeling/qsubtemp/tempscripts/'+jobname):\n")
        FILE.write("\tos.mkdir('/common/madsci/Modeling/qsubtemp/tempscripts/'+jobname)\n")
        FILE.write("if not os.path.exists('/common/madsci/Modeling/qsubtemp/'+jobname):\n")
        FILE.write("\tos.mkdir('/common/madsci/Modeling/qsubtemp/'+jobname)\n\n")
        #This sets up the list if repeating on a number of PDBs
        if repeat:
            FILE.write("FILE=open(List, mode='r')\n")
            FILE.write("pdbs=0\n")
            FILE.write("pdbList = []\n")
            FILE.write("for line in FILE:\n")
            FILE.write("\tline = line.rstrip('\\n')\n")
            FILE.write("\tpdbList.append(line)\n")
            FILE.write("\tpdbs+=1\n")
            FILE.write("FILE.close()\n")

        FILE.write("if len(sys.argv)<2:\n"+"\tprint \"No Config file specified. \"\n"+"\tos.abort()\n")
        FILE.write("config = sys.argv[1]\n")
        FILE.write("protocol = open(config, 'r')\n")
        FILE.write("prot = protocol.readline()\n")
        FILE.write("prot = prot.rstrip('\\n')\n")
        FILE.write("protocol.close()\n")
        FILE.write("offset = 0\n")
        
        #Here, we parse the script file, remove any outpaths (because we do that manually at the end of the job...), and set outpath for bash shell.
        FILE.write("configSP = prot.split()\n")
        FILE.write("outpath = None\n")
        FILE.write("for stuff in configSP:\n")
        FILE.write("\tif re.search('-out:path', stuff):\n")
        FILE.write("\t\tind = configSP.index(stuff)\n")
        FILE.write("\t\tconfigSP.pop(ind)\n")
        FILE.write("\t\toutpath = configSP[ind]\n")
        FILE.write("\t\tconfigSP.pop(ind)\n")
        
        FILE.write("startProt = ' '.join(configSP)\n")
        
        #Here, we loop to create each shell script for each job.
        if repeat:
            #Here, we repeat for each PDB:
            FILE.write("rootjob = jobname\n")
            FILE.write("x = 0\n")
            FILE.write("for pdbPATH in pdbList:\n")
            FILE.write("\tname = os.path.split(pdbPATH)[1].split('.')[0]\n")
            FILE.write("\tjobname = name+rootjob\n")
            
            
            #Here is the inner loop for jobs, so that you can do many jobs per pdb, or many structures per job per pdb.
            FILE.write("\tfor i in range(1, jobs+1):\n")
            FILE.write("\t\tprot = startProt\n")
            FILE.write("\t\tscript = '/common/madsci/Modeling/qsubtemp/tempscripts/'+rootjob+'/'+jobname+'_'+repr(i)+'.in'\n")
            FILE.write("\t\tJOBSCRIPT = open(script, 'w')\n")
            FILE.write("\t\tjran = str(1000000+offset)\n")
            #FILE.write("\t\tif method == 'Job x Total Structures':\n")
            #FILE.write("\t\t\tsuff = str(i)+str(stru)\n")
            #FILE.write("\t\telif method == 'job':\n")
            #FILE.write("\t\tsuff = str(i)\n")
            #FILE.write("\t\telse:\n")
            #FILE.write("\t\t\tprint 'Custom method not yet implemented.'\n")
            FILE.write("\t\tsuff = str(i)\n")
            FILE.write("\t\tprot = prot + ' -in:file:s '+pdbPATH+' -out:nstruct '+repr(stru)+ ' -constant_seed -jran '+jran+' -out:prefix '+jobname+' -out:suffix '+suff\n")
            FILE.write("\t\tJOBSCRIPT.write('#PBS -l nodes=1:ppn=1\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('mkdir /scratch/'+jobname+'_'+repr(x)+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('cd /scratch/'+jobname+'_'+repr(x)+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write(prot+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('mkdir '+outpath+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('mkdir '+outpath+'/raw'+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('cp -r /scratch/'+jobname+'_'+repr(x)+' '+outpath+'/raw'+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.write('rm -r /scratch/'+jobname+'_'+repr(x)+'\\n')\n")
            FILE.write("\t\tJOBSCRIPT.close()\n")
            FILE.write("\t\tj=repr(i)+'_'+rootjob\n")
            FILE.write("\t\tprint 'Kicking Job number '+repr(x)+'_'+jobname\n")
            FILE.write("\t\toffset = offset+10\n")
           #FILE.write("\t\tprint j\n")
            FILE.write("\t\trun = '/usr/local/bin/qsub -d /common/madsci/Modeling/qsubtemp/'+rootjob+' -N '+jobname+' -V -q '+queue+' '+script\n")
            FILE.write("\t\tos.system(run)\n")
            FILE.write("\t\tx+=1\n")
            self.textHelp.delete(1.0, END)
            self.textHelp.insert(1.0, "File Saved.")
        else:
            FILE.write("rootjob = jobname\n")
            FILE.write("for i in range(1, jobs+1):\n")
            FILE.write("\tprot = startProt\n")
            FILE.write("\tscript = '/common/madsci/Modeling/qsubtemp/tempscripts/'+rootjob+'/'+jobname+'_'+repr(i)+'.in'\n")
            FILE.write("\tJOBSCRIPT = open(script, 'w')\n")
            FILE.write("\tjran = str(1000000+offset)\n")
            #FILE.write("\tif method == 'Job x Total Structures':\n")
            #FILE.write("\t\tsuff = str(i)+str(stru)\n")
            #FILE.write("\telif method == 'job':\n")
            #FILE.write("\t\tsuff = str(i)\n")
            #FILE.write("\telse:\n")
            #FILE.write("\t\tprint 'Custom method not yet implemented.'\n")
            FILE.write("\tsuff = str(i)\n")
            FILE.write("\tprot = prot + ' -out:nstruct '+repr(stru)+ ' -constant_seed -jran '+jran+' -out:prefix '+jobname+' -out:suffix '+suff\n")
            FILE.write("\tJOBSCRIPT.write('#PBS -l nodes=1:ppn=1\\n')\n")
            FILE.write("\tJOBSCRIPT.write('mkdir /scratch/'+jobname+'_'+repr(i)+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write('cd /scratch/'+jobname+'_'+repr(i)+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write(prot+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write('mkdir '+outpath+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write('mkdir '+outpath+'/raw'+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write('cp -r /scratch/'+jobname+'_'+repr(i)+' '+outpath+'/raw'+'\\n')\n")
            FILE.write("\tJOBSCRIPT.write('rm -r /scratch/'+jobname+'_'+repr(i)+'\\n')\n")
            FILE.write("\tJOBSCRIPT.close()\n")
            FILE.write("\tprint 'Kicking Job number '+repr(i)\n")
            FILE.write("\toffset = offset+10\n")
            FILE.write("\trun = '/usr/local/bin/qsub -d /common/madsci/Modeling/qsubtemp/'+rootjob+' -N '+jobname+' -V -q '+queue+' '+script\n")
            FILE.write("\tos.system(run)\n")
            self.textHelp.delete(1.0, END)
            self.textHelp.insert(1.0, "File Saved.")
        FILE.close()
    def runScript(self):
        self.saveScript(F=1)
        config = tkFileDialog.askopenfilename(initialdir = self.pwd)
        command= 'python '+ self.pwd+"/tempScript.py "+config
        os.system(command)
        self.textHelp.delete(1.0, END)
        self.textHelp.insert(1.0, "Cluster run started...")
        return
    def runCurrent(self):
        #First, we have to check if there even is a current script to run.  How do we do this?
        pass
    def calcRunTime(self):
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
            stru = int(self.struc.get())
            total = str(jobs*stru)
            self.jobtotal.set(total)
            return int(total)
        except TypeError:
            pass
            
        
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
    MainWindow.title("Rosetta Cluster Setup")
    SetupWindow = ClusterSetup()
    SetupWindow.shoWindow(MainWindow, 0, 0)
    MainWindow.mainloop()
