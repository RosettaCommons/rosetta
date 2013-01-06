#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/qsub_script.py
## @brief  Script to run JD1 applications on a cluster using QSUB.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import os
import sys
import re
from optparse import OptionParser, IndentedHelpFormatter

"""
Script that Coincides/Works with with RosettaProtocolBuilder and QsubClusterSetup.
To use outside of these GUI's:
1) First line of config: Full path to executable to run. Specify after # symbol.
2) -database flag.  Full path to rosetta_database
3) Any command-line options.  Checkout optparse for a full description.
"""



def setup_paths(qsubtemp, tempscripts, jobname):
    """
    Setup qsubtemp and tempscript paths.  Including ones for the job.
    """
    
    print "Setting up paths in "+qsubtemp+" and "+tempscripts
    if not os.path.exists(qsubtemp):
            os.mkdir(qsubtemp)
    if not os.path.exists(tempscripts):
            os.mkdir(tempscripts)
    if not os.path.exists(tempscripts+'/'+jobname):
            os.mkdir(tempscripts+'/'+jobname)
    if not os.path.exists(qsubtemp+'/'+jobname):
            os.mkdir(qsubtemp+'/'+jobname)
    
    print "Paths Setup."
    
def location():
    """
    Allows the script to be self-aware of it's path.
    So that it can be imported from anywhere.
    """

    p = os.path.abspath(__file__)
    pathSP = os.path.split(p)
    return pathSP[0]
    
def remove_from_array(string, array):
    ind = array.index(string)
    array.pop(ind)
    return array
    
def write_jobscript(outpath, scriptPath, jobname, i, config):
    """
    Writes a Shell script for each job to be used Qsub.
    """
    
    JOBSCRIPT = open(scriptPath, 'w')
    JOBSCRIPT.write('#PBS -l nodes=1:ppn=1\n')
    JOBSCRIPT.write('mkdir /scratch')
    JOBSCRIPT.write('mkdir /scratch/'+jobname+'_'+repr(i)+'\n')
    JOBSCRIPT.write('cd /scratch/'+jobname+'_'+repr(i)+'\n')
    prot = " ".join(config)
    JOBSCRIPT.write(prot+'\n')
    JOBSCRIPT.write('mkdir '+outpath+'\n')
    JOBSCRIPT.write('mkdir '+outpath+'/raw'+'\n')
    JOBSCRIPT.write('cp -r /scratch/'+jobname+'_'+repr(i)+' '+outpath+'/raw'+'\n')
    JOBSCRIPT.write('rm -r /scratch/'+jobname+'_'+repr(i)+'\n')
    JOBSCRIPT.close()
    
def main(args):
    """
    OVERVIEW: Main function.  Parses from OptParse for backwards compatability.  Reads a config file, usually specified by the @ symbol.
    Creates a shell script for each job that will run.  Shell script will make a directory called /scratch (and /scratch/jobname).  Cd into it, and run Rosetta.
    After the run, It will copy all files into whatever is specified in -out:path (or ./RESULTS if not given), then delete /scratch/jobname
    NOTE: Will Create a /RAW directory for each job that needs to run.  This is because of JD1.
    NOTE: Use functions in main pyrosetta_toolkit GUI to combine these files, or your own script.
    """
    
    pwd = location()
    parser = OptionParser()
    parser.add_option("--config", "-c",
        default="./temp.config",
        help = "Path to the config file that you want to run on the cluster (file used by rosetta @ command)"
    )
    parser.add_option("--jobs","-j",
        type="int",
        help = "Number of Jobs to create "
    )
    parser.add_option("--stru",
        type="int",
        default=1,
        help = "Number of structure per job"
    )
    parser.add_option("--qsub",
        default="usr/local/bin/qsub",
        help = "Path to qsub executable"
    )
    parser.add_option("--qsubtemp",
        default = pwd+"/qsubtemp",
        help = "Path to where qsub should write all of its output.  If not given will make a directory in current directory"
    )
    parser.add_option("--tempscripts",
        default = pwd+"/tempscripts",
        help = "Path to where the sh scripts should go that begin each run on each processor.  if not given will make a directroy in current directory",
    )
    parser.add_option("--queue","-q",
        default="",
        help="Specific Queue to run qsub on.  Default is empty."
    )
    parser.add_option("--jobname",
        default="RosettaJob",
        help = "Name for the Job.  Used to name scripts and directories."
    )
    parser.add_option("--offset", "-o",
        type = "int",
        default = 10,
        help = "Offset to use for Seed of Job.  Will increment this with increasing jobs."
    )
    parser.add_option(
        "--pdb_list",
        default="",
        help = "List of PDBs to repeat on.  Replaces in:file:pdb.  (Very) Useful for testing protocols."
    )
    parser.add_option(
        "--debug",
        default="store_false",
        action="store_true",
        help = "Test the program.  Do not use os.system to kick off the job."
        
    )
    (options, args) = parser.parse_args(args=args[1:])

    if not os.path.exists(options.config):
        sys.exit("Config path wonky. "+options.config+" Please correct.")
    setup_paths(options.qsubtemp, options.tempscripts, options.jobname)
    FILE = open(options.config, 'r')
    config = []
    for line in FILE:
        line = line.rstrip()
        config.append(line)
    FILE.close()

    #Check to make sure executable and database option is set.
    exec_database=""; #First line of file.  Has executable and -database option present.
    if not re.search('#', config[0]):

        sys.exit("First line of file must have full path to executable after #\n" \
                 +"RosettaProtocolBuilder adds this line.  Please add the line (either with # or without) and run the program again.")
        
    else:
        if re.search("#", config[0]):
            exec_database = config[0]
            exec_database = exec_database.replace("#","")
            config[0]=exec_database
            #config.pop(0)
    
    
    #Remove and parse anything we need from the config file.
    flags_to_remove = ['#', '-out:nstruct', '-constant_seed', '-jran', '-nstruct']
    outpath = None;
    for option in config:
        if re.search('-out:path', option):
            outpath = option.split()[0]
            ind = config.index(option)
            config.pop(ind)
            if not os.path.exists(output):os.mkdir(outpath)
            continue
        #Places to get info back out  
        #Remove anything in the list from our config array.
        for flag in flags_to_remove:
            if re.search(flag, option):
                ind = config.index(option)
                config.pop(ind)
                continue
    if not outpath:
        print "No outpath specified in Config File.  Using current directory/RESULTS/jobname."
        outpath = pwd+'/RESULTS/'+options.jobname
        if not os.path.exists(pwd+'/RESULTS'):
            os.mkdir(pwd+'/RESULTS')
        if not os.path.exists(pwd+'/RESULTS/'+options.jobname):
            os.mkdir(pwd+'/RESULTS/'+options.jobname)
    #-mpi_tracer_to_file [a_file_stem]
    
    
    
    #Main Writing of Shell scripts.
    print "Removed Flags.  Created Results directory."
    config.append('-out:nstruct '+repr(options.stru))
    
    pdbList = []
    offset = options.offset
    if options.pdb_list:
        FILE=open(options.pdb_list, mode='r')
        pdb_count=0
        pdbList = []
        for line in FILE:
            line = line.rstrip()
            pdbList.append(line)
            pdb_count+=1
        FILE.close()
                #FILE.write("rootjob = jobname\n")
        x=0
        rootjob = options.jobname
        for pdbPATH in pdbList:
            name = os.path.split(pdbPATH)[1].split('.')[0]
            pdb_jobname = name+rootjob
            new_outpath = outpath+'/'+pdb_jobname
            if not os.path.exists(new_outpath):os.mkdir(new_outpath)
            
            for i in range(1, options.jobs+1):
                base_config = list(config)
                scriptname =options.tempscripts+'/'+rootjob+'/'+rootjob+'_'+pdb_jobname+'_'+repr(i)+'.in'
                jran = str(1000000+offset)
                suff = str(i)
                base_config.append('-in:file:s '+pdbPATH+' -run:constant_seed -run:jran '+jran)
                write_jobscript(new_outpath, scriptname, pdb_jobname, x, base_config)
                print 'Kicking Job number '+repr(x)+'_'+name
                offset = offset+10
                run = options.qsub+'-d '+options.qsubtemp+'/'+rootjob+' -N '+pdb_jobname+' -V -q '+options.queue+' '+scriptname
                if not options.debug:
                    os.system(run)
                x+=1
        
    else:

        for i in range(1, options.jobs+1):
            base_config = list(config); #Create NEW list, NOT a reference to the old one.
            print config
            scriptname = options.tempscripts+'/'+options.jobname+'/'+options.jobname+'_'+repr(i)+'.in'
            jran = str(1000000+offset)
            base_config.append('-constant_seed -jran '+jran+"\n")
            
            write_jobscript(outpath, scriptname, options.jobname, i, base_config)
            
            print 'Kicking Job number '+repr(i)
            offset = offset+10
            run = options.qsub+'-d '+options.qsubtemp+'/'+options.jobname+' -N '+options.jobname+' -V -q '+options.queue+' '+scriptname
            if not options.debug:
                os.system(run)

if __name__ == '__main__': main(sys.argv)