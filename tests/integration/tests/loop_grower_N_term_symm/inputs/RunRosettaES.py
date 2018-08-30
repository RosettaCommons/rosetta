#!/usr/bin/python
import sys
import argparse
import os
import multiprocessing
import time
import glob
import fileinput
import re
import shutil
import subprocess

def main():
    args = parseargs()
    check_inputs(args)
    if len(args.relaxpdbs) != 0:
        relaxpdbs(args)
    if args.dumpbeam == None and len(args.lpsfiles) == 0:
        run_segment(args)
    if args.dumpbeam != None:
        dump_beam(args,args.dumpbeam)
    if len(args.lpsfiles) != 0:
        run_monte_carlo_assembly(args)
        
def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',default='input.pdb', help='The input pdb')
    parser.add_argument('-x','--xml',default='RosettaES.xml',help='The xml to run')
    parser.add_argument('-d','--density',default='density.mrc',help='The electron density map')
    parser.add_argument('-f','--fasta',default='input.fasta',help='The sequence in fasta format')
    parser.add_argument('-name','--name',help='The name of the job')
    parser.add_argument('-l','--loop',type=int,default=0,help='Which loop are you operating on. Default will leave xml unchanged')
    parser.add_argument('-c','--cores',type=int,default=1,help='How many cores you have available per node')
    parser.add_argument('-t','--taboo',type=int,default=1,help='Do you want to use taboo filtering when available')
    parser.add_argument('-compare','--compare',type=int,default=1,help='Do you want to run the comparator when all jobs are finished?')
    parser.add_argument('-docker','--docker',type=int,default=1,help='Set to 1 if these arguments are being sent to the docker image or being run inside it')
    parser.add_argument('-rp','--relaxpdbs',nargs="+",default=[],help='The pdbs to relax')
    parser.add_argument('-rlxs','--relax_script',default='runrelax.sh',help='The script to call relax')
    parser.add_argument('-lps','--lpsfiles',nargs="+",default=[],help='The loop partial solution files used by the Monte Carlo Assembly')
    parser.add_argument('-rs','--runscript',default='runES.sh',help='The script to call RosettaES')
    parser.add_argument('-rd','--rosettadir',default='/Rosetta/',help='The Rosetta directory containing a bin/ and database/ directory')
    parser.add_argument('-db','--dumpbeam',help='The beam file you want to dump the pdbs of')
    parser.add_argument('-s','--setup',type=int,default=1,help='Do you want to pick fragments before running')
    args = parser.parse_args()
    return args

def check_inputs(args):
    is_missing = False
    if not os.path.isfile(args.fasta):
        print 'fasta file',args.fasta,'is missing'
        is_missing = True
    if not os.path.isfile(args.pdb):
        print 'pdb file',args.pdbfile,'is misssing'
        is_missing = True
    if not os.path.isfile(args.density):
        print 'density file',args.density,'is misssing'
        is_missing = True
    if not os.path.isfile(args.xml):
        print 'xml file',args.xml,'is missing'
        is_missing = True
    fragfiles = find_fragfiles(args.xml)
    for ffile in fragfiles:
        if not os.path.isfile(ffile):
            print 'fragment file',ffile,'is missing'
            is_missing = True
    if is_missing:
        exit()

def run_segment(args):
    if not os.path.isdir(args.name):
        os.mkdir(args.name)
    files = glob.glob('*')
    for cpfile in files:
        if not os.path.isdir(cpfile):
            shutil.copy(cpfile,args.name)
    os.chdir(args.name)
    growloop(args)

#pick the fragments for RosettaES
def setup(args):
    command = ('sh %sRosettaESsetup.sh %s %s'%(args.runscripts,args.pdb,args.fasta))
    subprocess.call(command, shell=True)
    if os.path.isfile('grower_prep.txt'):
        os.remove('grower_prep.txt')

#creates the command line argument and runs the individual ES job
def run_single_job(args,usebeam,beamfile,steps,count,usetaboo,taboobeams):
    command = [args.runscript,args.xml,args.pdb,args.density,usebeam,beamfile,steps,count,usetaboo,taboobeams]
    command = " ".join(str(x) for x in command)
    command = 'sh '+command
    if usebeam == 1:
        with open(beamfile,'r') as bfile:
            if len(bfile.readlines()) == 0:
                print 'beam file is empty exiting protcol'
                exit(0)
    subprocess.call(command, shell=True)

#setups the xml with correct fasta and is used to turn on or off monte carlo assembly and dumping of the structures
def update_xml(args,dumpbeam,read_from_file,finalbeam=1):
    for line in fileinput.input(args.xml, inplace=1):
        if args.loop != 0:
            line = re.sub("looporder=\S+", 'looporder="%s"'%args.loop, line)
        line = re.sub("dumpbeam=\S+", 'dumpbeam="%s"'%dumpbeam, line)
        line = re.sub("fasta=\S+", 'fasta="%s"'%args.fasta, line)
        line = re.sub("read_from_file=\S+", 'read_from_file="%s"'%read_from_file, line)
        line = re.sub("dumpfinalbeam=\S+", 'dumpfinalbeam="%s"'%finalbeam, line)
        sys.stdout.write(line)

#parse the xml to check for the fragment files. This is used to ensure they exist before launching jobs
def find_fragfiles(xml):
    fragfiles = []
    with open(xml,'r') as xmlfile:
        for line in xmlfile:
            if 'fragfile=' in line:
                frags = re.search('fragfile="(.*)"',line).group(1)
                fragfiles.append(frags)
    return fragfiles

#the main workhorse of RosettaES manages the parallelization and distribution of all jobs.
def growloop(args):
    if os.path.isfile('finished.txt'):
        print 'found finished'
        return
    if args.taboo == 1:
        if not os.path.isdir('taboo'):
            os.mkdir('taboo')
    update_xml(args,0,0)
    counter = get_counter()

    
    cleanprevious()
    
    #this first run generates doesn't need to be parallel and it generates the files required for the first round
    ready_for_parallel = False
    print 'ready for parallel'
    if counter == 1:
        beamfile = 'NA'
        usebeam = False
        filterbeam = False;
        filterbeams = 'NA'
        step = 1
        run_single_job(args,usebeam,beamfile,step,counter,filterbeam,filterbeams)
    else:
        splitjobs(counter,args.cores)
        ready_for_parallel = True

    done = False
    if os.path.isfile('finished.txt'):
        print 'found finished'
        done = True
    while done == False:
        begintime = time.time()
        if ready_for_parallel == False:
            combine(counter)
            beamcount = secondfilter(args,counter)
            splitjobs(counter,args.cores)
        ready_for_parallel = False
        runparalleljobs(args,counter)
        if args.taboo == 1:
            addtaboo(counter)
        counter+=1
        if(os.path.isfile("./finished.txt")):
           done = True
           print "found finished"
    if args.taboo == 1:
        addtaboo(counter)

    ##generate output pdbs
    beamfile = 'taboo/beam%s.txt'%counter
    dump_beam(args,beamfile)

    #args.relaxpdbs = glob.glob('final*.pdb')
    #relaxpdbs(args)
    #rename_pdbs_from_score('ESmodel','score.sc')

    #cleanup
    extra_files = glob.glob('beam*')
    extra_files.append('finished.txt')
    extra_files.append('for_taboo.txt')
    extra_files += glob.glob('final*.pdb')
    for extra_file in extra_files:
        if os.path.isfile(extra_file):
            os.remove(extra_file)

#turn on dump beam in the xml run the job then turn it back off.
def dump_beam(args,beamfile):
    update_xml(args,1,0)
    run_single_job(args,1,beamfile,0,0,0,0)
    update_xml(args,0,0)
   
#this file gets the current beam count by reading the beam_0.txt file (it no longer reads the count file which isn't necessary
def get_counter():
    counter = 1
    if os.path.exists('beam_0.txt'):
       with open('beam_0.txt') as f:
           counter = int(f.readline().split()[-2])
    return counter

#remove all beams
def cleanprevious():
    for jobfiles in glob.glob("beam_*_*.txt"):
        os.remove(jobfiles)
    for jobfiles in glob.glob('beam_*.*.txt'):
        os.remove(jobfiles)

def run(args,beamfile,jobcounter):
    previouscount = int(beamfile.split('_')[1])+1
    previousbeams = 'taboo/beam%s.txt'%(previouscount)
    filterbeam = 0
    if args.taboo == 1:
        if checktaboo(previouscount):
            filterbeam = 1
    steps = 1
    run_single_job(args,1,beamfile,steps,jobcounter,filterbeam,previousbeams)

#merge the results from each core
def combine(counter):
    combinedbeams = open("beam_%d.txt" % (counter), "w")
    for beamfile in glob.glob('beam_%d.*.txt'%counter):
        singleroundbeam = open(beamfile, "r")
        combinedbeams.write(singleroundbeam.read())
        singleroundbeam.close()
    combinedbeams.close()
    cleanprevious()

#filter the combined results
def secondfilter(args,counter): 
    beamfile = ("beam_%d.txt" % (counter))
    filterbeam = 0
    filterbeamfile = 'none'
    if args.taboo == 1:
        if checktaboo(counter):
            filterbeam = 1
            filterbeamfile = 'taboo/beam'+str(counter)+'.txt'
    run_single_job(args,1,beamfile,0,0,filterbeam,filterbeamfile)
    storedbeam = open("beam_0.txt", "r")
    beamcount = 0
    for line in storedbeam.readlines():
        if len(line.split()) == 5: 
            beamcount+=1
    if beamcount == 0:
        print 'error there are no solutions in the target beamfile exiting the script'
        exit()
    storedbeam.close()
    print beamcount
    return beamcount

#split the beamfiles up to be run on all the cores
def splitjobs(counter,cores):
    storedbeam = open("beam_0.txt", "r")
    divisioncount = 1
    partialset = open("beam_0.txt", "r")
    newbeam_linesize = 5
    for line in storedbeam.readlines():
        if len(line.split()) == newbeam_linesize:
            partialset.close()
            partialset = open("beam_%d_%d.txt" % (counter,divisioncount), "a")
            divisioncount+=1
        if divisioncount > cores:
            divisioncount = 1
        partialset.write(line)
    storedbeam.close()    
    partialset.close()   

#run the loop grower in parallel
def runparalleljobs(args,counter):
    jobcount = 1
    filecount = len(glob.glob('beam_%d_*.txt'%counter))
    while jobcount <= filecount:
          if __name__ == '__main__':
             jobs = []
             for i in range(1, args.cores+1):
                 p = multiprocessing.Process(target=run, args=(args,"beam_%d_%d.txt" % (counter,i),i))
                 jobs.append(p)
                 p.start() 
                 jobcount = jobcount + 1
                 if jobcount > filecount:
                     break
             for x in jobs:
                 x.join()

def checktaboo(count):
    return os.path.isfile('/taboo/beam%s.txt'%count)

#Check all the beams in the currentfile and in the taboo file. Add all unique beams to the taboo file.
def addtaboo(count):
    beamlist = []
    firstline = ''
    with open('beam_0.txt','r') as beamfile:
        beam = []
        for line in beamfile:
            if firstline == '':
                firstline = line
            if line == firstline:
                if len(beam) != 0:
                    beamlist.append(beam)
                    beam = []
            beam.append(line)
        beamlist.append(beam)
    
    taboobeamlist = []
    firstline = ''
    if checktaboo(count):
        with open('taboo/beam%s.txt'%count) as taboobeams:
            beam = []
            for line in taboobeams:
                if firstline == '':
                    firstline = line
                if line == firstline:
                    if len(beam) != 0:
                        taboobeamlist.append(beam)
                        beam = []
                beam.append(line)
            taboobeamlist.append(beam)
    uniquebeams = []
    for beam in beamlist:
        if beam not in taboobeamlist:
            uniquebeams.append(beam)
    with open('for_taboo.txt','w') as fortaboo:
        for beam in uniquebeams:
            for line in beam:
               fortaboo.write(line)
    os.system('cat for_taboo.txt >> taboo/beam%s.txt'%count)

def count_beams(beamfile):
    firstline = ''
    beamcount = 0
    with open(beamfile,'r') as bfile:
        for line in bfile:
            if firstline == '':
                firstline = line
            if firstline == line:
                beamcount+=1
    return beamcount

#Relax Functions
def runrelax(inputs):
    args = inputs[0]
    pdb = inputs[1]
    relaxscript = inputs[2]
    command = ('sh runrelax.sh %s %s'%(pdb,args.density))
    subprocess.call(command, shell=True)

def relaxpdbs(args):
    #make sure the score file is clean before we relax
    if os.path.isfile('score.sc'):
        os.system('rm score.sc')
    if __name__ == '__main__':
        jobs = []
        p = multiprocessing.Pool(args.cores)
        total_jobs = len(args.relaxpdbs)
        jobcount = 0
        while jobcount < total_jobs:
            pdb = args.relaxpdbs[jobcount]
            jobs.append((args,pdb,args.relax_script))
            jobcount+=1
        p.map(runrelax, jobs)

#renaming functionss
def rename_pdbs_from_score(prefix,scorefile):
    tag_to_score = []
    with open('score.sc','r') as sf:
        for line in sf:
            if line.startswith("SCORE:"):
                line = line.split()
                if isfloat(line[1]):
                    tag = line[-1]
                    score = float(line[1])
                    tag_to_score.append((score,tag))
    tag_to_score = sorted(tag_to_score)
    it = 1
    for pair in tag_to_score:
        mvcommand = 'mv %s.pdb %s_%s.pdb'%(pair[1],prefix,it,)
        print mvcommand
        os.system(mvcommand)
        it+=1

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

#Monte Carlo Assembly Function
def run_monte_carlo_assembly(args):
    combine_lps(args)
    update_xml(args,0,1)
    run_single_job(args,0,'na',0,0,0,'na')
    update_xml(args,0,0)
    #needsanother = needsmorerounds()
    #if needsanother == False:
    #    args.relaxpdbs = glob.glob('aftercomparator*')
    #    relaxpdbs(args)
    #rename_pdbs_from_score('MCAmodel','score.sc')

def combine_lps(args):
    with open("lpsfile.txt", "w") as combinedlps:
        loopcount = len(args.lpsfiles)
        combinedlps.write(str(loopcount)+ "\n")
        sorted_files = sort_lpsfiles(args)
        for lpspair in sorted_files:
            beamcount = 0;
            lpsfile = open(lpspair[1],"r")
            for line in lpsfile.readlines():
                if len(line.split()) == 1:
                    beamcount = beamcount + 1
            beamcount = beamcount/2        
            strbeamcount = str(beamcount)
            combinedlps.write(strbeamcount + "\n")
            lpsfile.seek(0)
            combinedlps.write(lpsfile.read())

def sort_lpsfiles(args):
    sorted_lps = []
    for lpsfile in args.lpsfiles:
        with open(lpsfile,'r') as lps:
            lower = int(lps.readline().split()[0])
            sorted_lps.append((lower,lpsfile))
    sorted_lps = sorted(sorted_lps)
    return sorted_lps

def needsmorerounds():
    needs_another = True
    if os.path.isfile('recommendation.txt'):
        with open('recommendation.txt') as recommendation:
            for line in recommendation:
                if '''Another round probably isn't necessary''' in line:
                    needs_another = False
    return needs_another

main()
