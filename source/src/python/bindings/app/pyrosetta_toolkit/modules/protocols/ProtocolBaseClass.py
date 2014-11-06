#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/ProtocolBaseClass.py
## @brief  Simple base class for running protocols in the GUI. Includes multiprocessing. 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
import rosetta.basic.options

#Python Imports
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Queue
from multiprocessing import Process
import time
import os.path

#Tkinter Imports



#Toolkit Imports
#from window_main.IO.GUIInput import GUIInput
from app.pyrosetta_toolkit.window_main.IO.GUIOutput import GUIOutput
from app.pyrosetta_toolkit.window_modules.scorefunction.ScoreFxnControl import ScoreFxnControl

class ProtocolBaseClass:
    def __init__(self, pose, score_class, input_class, output_class):
        self.pose = pose
        
        self.score_class = score_class
        self.input_class = input_class
        self.output_class = output_class
        
        #Ignore these.  These are a trick to hint Komodo IDE about what types each of these variables are. And it fucking works!!!!!!!! Hallelujah!!!
        #They never actually run!
        if 0:
            #self.input_class = GUIInput()
            self.output_class = GUIOutput()
            self.score_class = ScoreFxnControl()
        
    def run_protocol(self, mover):
        """
        Runs the protocol using the multiprocessing module.  Assumes the protocol has a mover associated with it.  Runs apply.
        If your protocol has no mover, You can create one in python from the Rosetta Base class.  Or, simply overwrite this method.
        """
        
        #if create_master:
        #    master = multiprocessing.Process(name="master", target=self.run_protocol(mover, False))
        #    master.start()
        #    if not master.is_alive: self.output_class.terminal_output.set(0)
        #    return
        
        self.pdb_name = self.output_class.outdir.get()+"/"+self.output_class.outname.get()
        start_energy_score = self.score_class.score(self.pose)
        
        self.output_class.terminal_output.set(1); #Redirect to stdout. For multiprocessing and major Rosetta output.
        if self.output_class.processors.get()>1:
            
            #Multiprocessing is done manually due to mover being unpicklable - Hence no Queue or Pool objects.
            #First, we have an array of jobs:
            workers = []
            for i in range(1, self.output_class.decoys.get()+1):
                outname = self.pdb_name+"_"+repr(i)+".pdb"
                worker = Process(name = repr(i), target=self._run_mover, args=(mover, outname))
                workers.append(worker)
            total_allowed_jobs = self.output_class.processors.get()
            print "Total allowed jobs: "+repr(total_allowed_jobs)
            total_running_jobs = 0
            job_complete=False
            
            #Check if any PDB's already exist.  Delete them or pop the workers depending on the overwrite option:
            for worker in workers:
                if os.path.exists(self.pdb_name+"_"+worker.name+".pdb"):
                    if self.output_class.overwrite.get():
                        os.remove(self.pdb_name+"_"+worker.name+".pdb")
                        print "Overwriting "+self.pdb_name+"_"+worker.name+".pdb"
                    else:
                        workers.pop(workers.index(worker))
                        print self.pdb_name+"_"+worker.name+".pdb already exists.  Skipping.  "
                        
            #Run the protocol
            while not job_complete:
                
                time.sleep(5)
                for worker in workers:
                    if worker.is_alive():
                        pass
                        #print "Worker is alive"
                        #total_running_jobs+=1; #Increment total_running_jobs
                    elif os.path.exists(self.pdb_name+"_"+worker.name+".pdb"):
                        if worker.exitcode!=0:
                            print "%s.exitcode = %s" %(worker.name, worker.exitcode)
                        
                        workers.pop(workers.index(worker)); #If the job is done, pop it.
                        total_running_jobs-=1
                        print "Total running jobs: "+repr(total_running_jobs)
                        print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                        
                if len(workers)==0:
                    job_complete=True
                    break
                

                
                if total_running_jobs<total_allowed_jobs:
                    for worker in workers:
                        if ((not worker.is_alive())and (not os.path.exists(self.pdb_name+"_"+worker.name+".pdb"))):
                            print "Starting Worker"
                            try:
                                worker.start()
                            except AssertionError:
                                continue
                            total_running_jobs+=1
                            print "Total running jobs: "+repr(total_running_jobs)
                            print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                            if total_running_jobs>=total_allowed_jobs: break
            
                if total_running_jobs==0:
                    job_complete=True
                        
            while len(multiprocessing.active_children()) != 0: time.sleep(1)
            
        else:
            mc = MonteCarlo(self.pose, self.score_class.score, self.output_class.kT)
            for i in range(1, self.output_class.rounds.get()+1):
                
                mover.apply(self.pose)
                if self.output_class.use_boltzmann.get():
                    if mc.boltzmann(self.pose):
                        print "MC: Pose Accepted"                  
                    
                    
                if self.output_class.rounds.get()>1:
   
                    print "Round"+repr(i)+": "+ repr(self.score_class.score(self.pose))+" REU"
            
            if self.output_class.recover_low.get():
                mc.recover_low(self.pose)
            
            mc.show_scores()
            mc.show_counters()
            print "Start: "+ repr(start_energy_score)+" REU"
            print "End: "+repr(self.score_class.score(self.pose))+" REU"
            
            #Here we output the decoy if the user wishes, overwriting if nessessary.
            if self.output_class.decoys.get()==1:
                
                outpath = self.pdb_name+"_1.pdb"
                if os.path.exists(outpath):
                    if self.output_class.overwrite.get():
                        os.remove(outpath)
                        self.output_pose(self.pose, outpath)
                        print "PDB exists.  Overwriting."
                    else:
                        print "PDB exists.  Skipping output."
                else:
                    self.output_pose(self.pose, outpath)
                print "Setting decoy as current pose."    
                
        self.output_class.terminal_output.set(0); #Reset output to textbox
        if self.output_class.decoys.get()>1:
            print "Original decoy still loaded."
            
        print "Job Complete.  Any decoys written to output directory."
        return
    
        
        
    def _run_mover(self, mover, outputname):
        """
        Used for multiprocessing.  
        """
        #Reinitialize Rosetta to reinitialize options, and specifically, the SEED.
        self.input_class.options_manager.re_init()
        
        p = Pose(); #Copy it so that each process is working on a different pose object. Probably unnessessary.
        p.assign(self.pose)
        print outputname
        start = self.score_class.score(p)
        mc = MonteCarlo(p, self.score_class.score, self.output_class.kT)
        for x in range(1, self.output_class.rounds.get()+1):
            mover.apply(p)
            if self.output_class.use_boltzmann.get():
                if mc.boltzmann(p):
                    print "MC: Pose Accepted"
            
        
        if self.output_class.recover_low.get():
            mc.recover_low(p)
            
        self.output_pose(p, outputname)
        
        print "Start: " +repr(start)+" REU"
        print "End: " +repr(self.score_class.score(p))+" REU"    
            
            
    def output_pose(self, p, outputname):
        """
        Output pose function for ProtocolBaseClass.
        """
        p.dump_pdb(outputname)
        
        score_tag = ".fasc"
        if not p.is_fullatom():
            score_tag = ".sc"

        scorefile = self.pdb_name + score_tag
        output_scorefile(p, self.input_class.pdb_path.get(), outputname, scorefile, self.score_class.score, self.output_class.decoys.get(), self.pose)
