from argparse import ArgumentParser
import os
from os import system
from os.path import basename, exists, splitext
from sys import exit, stderr, stdout
from itertools import repeat
import multiprocessing as mp
import glob
import copy
import shutil
import math

def collect_xyz(location):
    cwd = os.getcwd()
    data = []
    results = glob.glob("{}/*.out".format(location))
    os.chdir(location)
    for fn in results:
        basename, extension = os.path.splitext(fn)
        results_name = os.path.split(basename)[-1]
        os.system("column_extract.pl {}.out -regex SCORE -col X -col Y -col Z -col description> tmp.txt".format(results_name))
        data_in = open("tmp.txt")
        data_in.readline()
        for line in data_in:
            x = float(line.split(" ")[0])
            y = float(line.split(" ")[1])
            z = float(line.split(" ")[2])
            descr = line.split(" ")[3].strip()
            data.append([x,y,z,descr])
    os.chdir(cwd)
    #print(data)
    return(data)
def append_to_pdb(data,input_pdb,output_pdb):
    pdb_in = open(input_pdb)
    pdb_out = open(output_pdb,"w")
    for line in pdb_in:
        pdb_out.write(line)
    ct=1
    for f in data:
        x = f[0]
        y = f[1]
        z = f[2]
        pdb_out.write("ATOM   4716 2HE2 H2O C {: 4d}   {: 4.3f} {: 4.3f}  {: 4.3f}  1.00  0.00           H\n".format(ct,x,y,z))
        ct+=1

def get_dist(pt1,pt2):
    return(math.sqrt((pt1[0]-pt2[0])**2+(pt1[1]-pt2[1])**2+(pt1[2]-pt2[2])**2))

def select_area(data,pts,dist):
    limited_data = []
    for pt_data in data:
        min_dist = 10
        added=False
        for pt_selection in pts:
            dist = get_dist(pt_data,pt_selection)
            print(dist)
            if(dist<min_dist and not added):
                print("here")
                limited_data.append(pt_data)
                added=True
    return(limited_data)

def output_lego_names(data):
    limited_data = []
    for pt_data in data:
        print("{}".format(pt_data[3]))



pts = [[-30.3,23.11,-67.800],[-31.900,20.17,-62.5],[-26.100,34.050,-78.800]]
data = collect_xyz("results")
dist = 10 #in ang
limited_data = select_area(data,pts,10)
print(len(data))
print(len(limited_data))
append_to_pdb(limited_data,"ANK12_refine_33_AB_0001.pdb","output.pdb")
output_lego_names(limited_data)

