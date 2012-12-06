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

#Python Imports
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Queue
from multiprocessing import Process
import time
import os.path

#Tkinter Imports



#Toolkit Imports


class ProtocolBaseClass:
    def __init__(self, pose, score_class, input_class, output_class):
        self.pose = pose
        
        self.score_class = score_class
        self.input_class = input_class
        self.output_class = output_class
        
        self.py_jd = PyJobDistributor(output_class.outdir.get()+"/"+output_class.outname.get(), output_class.decoys, score_class.score)
        #self.workers Pool of workers.  Added to class at run_protocol.
        
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
        
        self.py_jd.nstruct = self.output_class.decoys.get()
        self.py_jd.scorefxn = self.score_class.score
        self.pdb_name = self.output_class.outdir.get()+"/"+self.output_class.outname.get()
        start_energy_score = self.score_class.score(self.pose)
        self.output_class.terminal_output.set(1); #Redirect to stdout. For multiprocessing and major Rosetta output.

        if self.output_class.auto_write.get():
            
            #First, we have an array of jobs:
            workers = []
            for i in range(1, self.output_class.decoys.get()+1):
                outname = self.pdb_name+"_decoy_"+repr(i)+".pdb"
                worker = Process(name = "decoy_"+repr(i), target=self._run_mover, args=(mover, outname))
                workers.append(worker)
            
            total_allowed_jobs = self.output_class.processors.get()
            print "Total allowed jobs: "+repr(total_allowed_jobs)
            total_running_jobs = 0
            job_complete=False
            
            #Check if any PDB's already exist:
            for worker in workers:
                if os.path.exists(self.pdb_name+"_"+worker.name+".pdb"):
                    workers.pop(workers.index(worker))
                    
            #Run the protocol
            while not job_complete:
                print "1"
                print "Total running jobs: "+repr(total_running_jobs)
                print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                
                time.sleep(5)
                for worker in workers:
                    if worker.is_alive():
                        pass
                        #print "Worker is alive"
                        #total_running_jobs+=1; #Increment total_running_jobs
                    elif os.path.exists(self.pdb_name+"_"+worker.name+".pdb"):
                        print "%s.exitcode = %s" %(worker.name, worker.exitcode)
                        
                        workers.pop(workers.index(worker)); #If the job is done, pop it.
                        total_running_jobs-=1
                        
                if len(workers)==0:
                    job_complete=True
                    break
                
                print "2"
                print "Total running jobs: "+repr(total_running_jobs)
                print "Total workers waiting: "+repr(len(workers)-total_running_jobs)
                
                if total_running_jobs<total_allowed_jobs:
                    for worker in workers:
                        if ((not worker.is_alive())and (not os.path.exists(self.pdb_name+"_"+worker.name+".pdb"))):
                            print "Starting Worker"
                            try:
                                worker.start()
                            except AssertionError:
                                continue
                            total_running_jobs+=1
                            if total_running_jobs>=total_allowed_jobs: break
            
                if total_running_jobs==0:
                    job_complete=True
                        
                
        else:
            for i in range(1, self.output_class.rounds.get()+1):
                print "Round "+repr(i)
                mover.apply(self.pose)
                print "Start: "+ repr(start_energy_score)+"\n"        
                print "End: "+ repr(self.score_class.score(self.pose))
        
        self.output_class.terminal_output.set(0); #Reset output to textbox
        
        print "NOTE: If > 1 decoy has been created, original decoy is still loaded. "
        print "Job Complete."
        
        return
    
        
        
    def _run_mover(self, mover, outputname):
        """
        Used for multiprocessing.  
        """
        p = Pose(); #Copy it so that each process is working on a different pose object.
        p.assign(self.pose)
        print outputname
        start = self.score_class.score(p)
        for x in range(1, self.output_class.rounds.get()+1):
            #print "Round "+repr(x)
            mover.apply(p)
        p.dump_pdb(outputname)
        print "Start: " +repr(start)
        print "End: " +repr(self.score_class.score(p))
        if self.output_class.decoys.get()==1:
            self.pose.assign(p)
            
            