# Constraint-related functions
# Emily Koo

# Import python modules
import os
from os.path import exists

# Import Rosetta modules
from rosetta import *

def load_constraints(PDB):
    # Get cst files
    sol = 'ls '+PDB+'_sol*.cst > sol_cst'
    ads = 'ls '+PDB+'_ads*.cst > ads_cst'
    
    os.system(sol)
    os.system(ads)
    
    sol_file = []
    ads_file = []
    if exists("sol_cst"):
        sol_cst = open('sol_cst','r')
        for line in sol_cst.readlines():
            if line.strip() is not "":
                sol_file.append(line)
    if exists("ads_cst"):
        ads_cst = open('ads_cst','r')
        for line in ads_cst.readlines():
            if line.strip() is not "":
                ads_file.append(line)

    rm = 'rm sol_cst ads_cst'
    os.system(rm)
    
    return (sol_file, ads_file)    

def apply_constraints(pose, cst):

    pose.remove_constraints()
    constr = rosetta.core.scoring.constraints
    set_constraints = ConstraintSetMover()

    if cst is not None:
        cst_set = constr.ConstraintSet()
    
        for file in cst:
            constr.ConstraintIO.read_constraints(file.strip(), cst_set, pose)
            
        pose.constraint_set(cst_set)
    
def load_disulf(disulf_file):

    disulf = core.io.raw_data.DisulfideFile(disulf_file)
    return disulf
        
def apply_disulf(pose, disulf):

    if disulf is not None:
        disulf.read_in_and_set_disulfides(pose)
        
def make_ads_disulf(pose, disulf):
    
    ads = str(disulf) + "_ads"
    if not exists(ads):
    
        surface_res = pose.num_jump()
        sol_file = open(disulf, 'r')
        ads_file = open(ads, 'a')
        
        for line in sol_file:
            newline = ""
            for res in line.split():
                newres = surface_res + int(res)
                newline = newline + str(newres) + " "
            ads_file.write(newline + "\n")

        sol_file.close()
        ads_file.close()
        
