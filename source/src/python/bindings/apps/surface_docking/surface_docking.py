#!/usr/bin/python2.6
# Main protocol for PyRosettaSurface 
# Emily Koo

"""
Usage: surface_docking.py @flags
flags is a text file containing the arguments (one on each line. see file in demo for example)
To see arguments, see surface_docking.py --help
"""
# Import python modules
import os
from sys import argv
import argparse
import random
import datetime

# Import Rosetta modules
from rosetta import *
from rosetta.protocols.rigid import *
import rosetta.core.scoring.solid_surface
import rosetta.core.scoring.constraints
import rosetta.core.io.raw_data

# Import custom modules
from surf_param import *
from append_data import *
from constraints import *
from movers import *

# Script options
parser = argparse.ArgumentParser(description='Surface Docking Script',fromfile_prefix_chars='@')

parser.add_argument('-s', '--start', action="store", dest="s", metavar='INPUT_FILE', required=True, help='PDB file INPUT_FILE containing protein and surface')
parser.add_argument('--database', action="store",  dest="database", metavar='DATABASE_FILE', required=True, help='Rosetta database path in PyRosetta')
parser.add_argument('--f3', '--frag3', action="store", dest="f3", metavar='FRAG3_FILE', required=False, default='', help='Use 3-mer fragment file FRAG3_FILE') 
parser.add_argument('--f9', '--frag9', action="store", dest="f9", metavar='FRAG9_FILE', required=False, default='', help='Use 9-mer fragment file FRAG9_FILE')
parser.add_argument('-c', '--constraints', action='store_true', required=False, default=False, help='Use ssNMR constraints')
parser.add_argument('-d', '--disulf', action="store", dest="d", metavar='DISULF', required=False, default='', help='Use disulfide constraints file DISULF')
parser.add_argument('--nosmallshear', action='store_true', required=False, default=False, help='Disable small/shear refinements completely (e.g. for rigid-body docking ')
parser.add_argument('-n', '--nstruct', action="store", dest="n", metavar='N', required=True, type=int, help='Generate N decoys')

args = parser.parse_args(['@flags'])
print args
input = args.s
path = args.database
constraints = args.constraints
decoy_num = int(args.n)
frag3 = args.f3
frag9 = args.f9
disulf_file = args.d
nosmallshear_ref = args.nosmallshear

if path.startswith("$"):
    newpath = os.path.expandvars(path)
elif path.startswith("~"):
    newpath = os.path.expanduser(path)
else:
    newpath = path

# Rosetta options
opts = ['app','-database',newpath,'-ex1','-ex2aro','-multiple_processes_writing_to_one_directory',\
        '-mute','core','-mute','protocols','-mute','basic']

print opts
Rargs = rosetta.utility.vector1_string()
Rargs.extend(opts)
rosetta.core.init(Rargs)

#========================== FUNCTIONS ===================================
def input_checker(input, constraints, frag3, frag9, disulf_file):
    # Input Error Checking
    if not exists(input):
        print "Input file does not exist."
        sys.exit()
        
    if constraints:
        cst = 'ls *.cst > csts'
        os.system(cst)

        if os.stat('csts').st_size == 0:
            print "Constraint files do not exist."
            rm = 'rm csts'
            os.system(rm)
            sys.exit()

    if frag3 == "" and frag9 != "" or frag3 != "" and frag9 == "":
        print "You must have both fragment files or no fragment files."
        sys.exit()

    elif frag3 != "" and frag9 != "":
        if not exists(frag3) or not exists(frag9):
            print "One or more fragment files do not exist."
            sys.exit()
            
    if disulf_file != "":
        if not exists(disulf_file):
            print "Disulfide file does not exist."
            sys.exit()
        
def run(state):
# Repeat full-atom relax outer-cycle times while
# ramping fa_rep up and ramping max_angle down

    fa_relax.state = state
    
    if fa_relax.outer_cycles > 1:
        rep = 0.02
        rep_inc = (0.44-rep)/(fa_relax.outer_cycles-1)
    else:
        rep = 0.44
        rep_inc = 0
    
    for outer_cycle in range(1, fa_relax.outer_cycles + 1):

        fa_relax.outer_cycle = outer_cycle
        fa_relax.set_max_angle(30/(outer_cycle))
        fa_relax.score_high.set_weight(fa_rep, rep)
        rep += rep_inc

        fa_relax.apply(pose)
        
def append_data(name, pose, scorefxn, PDB, state):
# Append additional data to the end of each PDB file

    append_scores(name, pose, scorefxn)
    append_hbonds(name, pose, scorefxn)

    apply_constraints(pose, load_constraints(PDB)[state])
    append_constraints(name, pose, load_constraints(PDB)[state], scorefxn)
    
#======================= USER SETTINGS =======================================

sol_cycles = randint(1,5)
sol_outer_cycles = 5
sol_inner_cycles = 5
ads_cycles = 5 - sol_cycles 
ads_outer_cycles = 5
ads_inner_cycles = 5

#========================= Pose and Variables Setup ==========================

input_checker(input, constraints, frag3, frag9, disulf_file)

print "Initializing poses..."
start_pose = Pose()
PDB = input[:-4]
pose_from_pdb(start_pose, input)
if disulf_file != '':
    make_ads_disulf(start_pose, disulf_file)

# Default score functions
hackelec = 1.0
score_high = create_score_function("score12")
score_high.set_weight(hack_elec, hackelec)
score_high.set_weight(fa_pair, 0.0) # fa_pair duplicates what hack_elec scores
score_high.set_weight(atom_pair_constraint, 1.0)
score_high.set_weight(dihedral_constraint, 1.0)
        
score_pack = create_score_function("score12")
score_pack.set_weight(hack_elec, hackelec)
score_pack.set_weight(fa_rep, 0.44)
score_pack.set_weight(fa_pair, 0.0) # fa_pair duplicates what hack_elec scores
score_pack.set_weight(atom_pair_constraint, 1.0)
score_pack.set_weight(dihedral_constraint, 1.0)
    
std_scorefxn = create_score_function("score12")
std_scorefxn.set_weight(hack_elec, hackelec)
std_scorefxn.set_weight(fa_rep, 0.44)
std_scorefxn.set_weight(fa_pair, 0.0) # fa_pair duplicates what hack_elec scores
std_scorefxn.set_weight(atom_pair_constraint, 1.0)
std_scorefxn.set_weight(dihedral_constraint, 1.0)

switch_low = SwitchResidueTypeSetMover('centroid')
switch_high = SwitchResidueTypeSetMover('fa_standard')
        
#=========================Job Distributor===============================
# Default procedure with fragment files:
# 1) Ab initio (Solution state)
# 2) Centroid relation (Solution state)
# 3) Full atom relaxation (Solution state)
# 4) Full atom relaxation (Adsorbed state)

# Add folder name to output files
cwd1 = os.getcwd()
cwd2 = cwd1.split('/')
cwdn = len(cwd2)
dir = cwd2[cwdn - 1]
print 'Current folder: ', dir

jd = PyJobDistributor("AdsState_" + PDB + "_" + dir, decoy_num, score_high)
jd1 = PyJobDistributor("SolState_" + PDB + "_" + dir, decoy_num, score_high)

while (jd.job_complete == False):
    #print "#######################################################"
    #print "################# GENERATING DECOY", jd1.current_num,"##################"
    #print "#######################################################"

    x = jd1.current_name
    y = jd.current_name

    time_start = datetime.datetime.now()
    
    pose = Pose()

    if frag3 != '' and frag9 != '':
        combined_pose = Pose()
        combined_pose.assign(start_pose)
        pose.assign(start_pose.split_by_chain()[2])
     
        switch_low.apply(pose)
        # Disulfide file with just protein atoms
        apply_disulf(pose, load_disulf(disulf_file))

        # 1) Abinitio
        print "-------------Starting ab initio protocol-----------"
        abinitio = Abinitio(PDB, frag3, frag9)
        abinitio.apply(pose)
        print "-------------- End ab initio protocol--------------"
        
        # 2) Centroid relaxation
        print "--------Starting centroid relax protocol-----------"
        cen_relax = CentroidRelax()
        cen_relax.apply(pose)
        print "---------End centroid relax protocol--------------"
        
        switch_high.apply(pose)

        # Combine protein pose with surface
        combined_pose.copy_segment(pose.total_residue(), pose, combined_pose.num_jump() + 1, 1)
        pose.assign(combined_pose)
    
    else:
        pose.assign(start_pose)
        switch_high.apply(pose)

    # Disulfide file including surface atoms
    apply_disulf(pose, load_disulf(disulf_file+"_ads"))

    print ">> Starting full atom relaxation..."
    # Initial settings
    fa_relax = FullAtomRelax(score_high, score_pack, std_scorefxn, nosmallshear_ref) 
    fa_relax.set_params(pose)   
    fa_relax.loadSurf(input)
    fa_relax.constraints = constraints

    # Full repack before starting
    fa_relax.fullRepack(pose) 
    
    print "-------------- Starting solution state refinement --------------"
    # Solution state settings
    fa_relax.name = x    

    if constraints:
        fa_relax.loadConstraints(PDB)    
        fa_relax.applyConstraints(pose, "sol")   

    fa_relax.ref_cycles = sol_cycles
    fa_relax.outer_cycles = sol_outer_cycles
    fa_relax.inner_cycles = sol_inner_cycles
    
    # 3) Solution state refinement
    for cycle in range(1, sol_cycles + 1):
        fa_relax.curr_cycle = cycle
        print "Solution state cycle ", cycle
        run("sol")
 
    # Full repack at the end        
    fa_relax.fullRepack(pose) 
    
    # Get solution state pose from FullAtomRelax protocol
    sol_pose = Pose()
    sol_pose.assign(pose)
    print "-------------- End of solution state refinement --------------"

    # Slide protein into contact with surface
    fa_relax._slideProt(pose)
            
    print "-------------- Starting adsorbed state refinement --------------"
    # Adsorbed state settings
    fa_relax.name = y

    if constraints:
        fa_relax.applyConstraints(pose, "ads")  
        
    fa_relax.ref_cycles = ads_cycles
    fa_relax.outer_cycles = ads_outer_cycles 
    fa_relax.inner_cycles = ads_inner_cycles
    
    # 4) Adsorbed state refinement
    for cycle in range(1, ads_cycles + 1):
        fa_relax.curr_cycle = cycle
        print "Adsorbed cycle ", cycle
        run("ads")
    
    # Full repack at the end        
    fa_relax.fullRepack(pose) 
    
    # 2 extra refinement cycles at rep = 0.44 for unbiased prediction
    if fa_relax.constraints == 'False':
        fa_relax.score_high.set_weight(fa_rep, 0.44)
        fa_relax.set_max_angle(6)
        
        # Starting at 1 would lead to random orientation
        for outer_cycle in range(2, 4):
            fa_relax.outer_cycle = outer_cycle
            fa_relax.apply(pose)
       
        fa_relax.fullRepack(pose) 

    print "-------------- End of adsorbed state refinement --------------"
    
    # Re-save it so it's not overwritten by next job
    curr_sol = str(x)
    curr_ads = str(y)

    # Output PDBs
    jd1.output_decoy(sol_pose)
    jd.output_decoy(pose)

    append_surf_sol = "grep SURFA " + input + " >> " + curr_sol 
    append_surf_ads = "grep SURFA " + input + " >> " + curr_ads
    
    os.system(append_surf_sol)
    os.system(append_surf_ads)

    append_data(curr_sol, sol_pose, score_high, PDB, 0)
    append_data(curr_ads, pose, score_high, PDB, 1)

    time_end = datetime.datetime.now()
    time_diff = time_end - time_start
    #file = open(curr_ads + ".energies.txt", 'a')
    #file.write("*\t"+str(time_diff))
    #file.close()
