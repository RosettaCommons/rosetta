#!/usr/bin/env python

# Edited by Emily Koo to adjust color bar scale automatically

import os
from IO import load
from sys import argv, exit
from Equations import vector, length
from PDBInfo import *
from numpy import zeros, arange

def ContactList():
   PDBs = load('PDBs')
   MinimumDistances = []
   for Decoy in PDBs:
      resetlists()
      Prot = partner(load(Decoy))[0]
      Surf = partner(load(Decoy))[1]

      for Residue in ResidueList(Decoy):
         Distance = []
         for Atom in GetResidueAtoms(Decoy, Residue):
            for SurfAtom in Surf:
               Distance.append(length(vector(Atom.XYZ, SurfAtom.XYZ)))

         MinDist =  min(Distance),Residue
         MinimumDistances.append(MinDist)

   return MinimumDistances

def MAIN():

    if os.stat('PDBs').st_size == 0:
        cmd = 'rm PDBs'
        os.system(cmd)

        exit()

    X = 0
    Y = 1
    total_res = len(ResidueList(load('PDBs')[0]))+1
    Resolution = [total_res, 300]
    Histogram = zeros(Resolution)
    Points = []
    for Line in ContactList():
        Contact = float(Line[1]), Line[0]
        Points.append(Contact)

    start_res = (Points[0])[0] - 1

    for Line in Points:
        Residue  = Line[X] - start_res
        Distance = Line[Y]
        Histogram[Residue,Distance] += 1

    Index = [arange(1,total_res,1)+1,arange(0,30,1)]

    FileOut = 'SurfaceContactMap.dat'
    o = open(FileOut, 'w')
    for x in enumerate(Index[X]):
        for y in enumerate(Index[Y]):
           o.write('%1i %1s %1i %1s %1i \n' %(x[1],' ', y[1], ' ', Histogram[x[1],y[1]]))

    # Generate GNU script to make the plot
    get_files = os.listdir('.')
    num_pdb = 0
    for files in get_files:
        if files.endswith(".pdb"):
            num_pdb += 1

    x=0.8
    range=len(ResidueList(load('PDBs')[0]))+0.5
    GNUPlotFile = 'SurfaceContactMap.gnuplot'
    O = open(GNUPlotFile, 'w')
    O.write('%0s \n' %('set view map'))
    O.write('%0s \n' %('set terminal png size 800, 800 nocrop enhanced font "VeraBd, 15"'))
    O.write('%0s \n' %('set output "SurfaceContactMap.png"'))
    O.write('%0s \n' %('set size square'))
    O.write('%0s %0.1f %0s \n' %('set xrange [0.5:',range,']'))
    O.write('%0s %0.1f %0s \n' %('set yrange [0.5:',30,']'))
    O.write('%0s \n' %('set palette defined (0 "#FFFFFF", 1 "#CCFFFF", ' + str(num_pdb) + ' "#0000FF")'))
    O.write('%0s \n' %('unset key'))
    O.write('%0s \n' %('set cbrange [0:' + str(num_pdb) + ']'))
    O.write('%0s \n' %('set cbtics (0, ' + str(num_pdb/4) + ', ' + str(num_pdb/2) + ', ' + str(3*num_pdb/4) + ', ' + str(num_pdb) + ') nomirror scale 0'))
    O.write('%0s \n' %('set colorbox user origin 0.875, 0.1625 size 0.025, 0.6925'))
    O.write('%0s \n' %('set xtic out nomirror'))
    O.write('%0s \n' %('set ytic out nomirror'))
    O.write('%0s \n' %('set mxtics'))
    O.write('%0s \n' %('set mytics 5'))
    O.write('%0s \n' %('set cblabel "Contact Frequency" font "VeraBd, 15"'))
    O.write('%0s \n' %('set border lw 3'))

    protein = partner(load(load('PDBs')[0]))[0]
    for Atom in protein:
        if Atom.AtmTyp == 'CA':
            # Gamma symbol in gnuplot
            if Atom.ResTyp == 'Z':
                output_res = "{/Symbol g}"
            else:
                output_res = AminoAcidTable(Atom.ResTyp)
            O.write('%0s %1.1f %0s %1.1f %0s \n' %("set label " '"'+output_res+'"' " at ",x,",",30.5, ' font "VeraBd, 9"'  ))
            x+=1

    O.write('%0s \n' %('set xlabel "Residue Number" font "VeraBd, 20"'))
    O.write('%0s \n' %('set ylabel "Distance (Angstroms)" font "VeraBd, 20"'))
    #O.write('%0s \n' %('set title "Osteocalcin/HAp [100] (Adsorbed)" font "VeraBd, 20"'))
    O.write('%0s \n' %('splot "SurfaceContactMap.dat" u 1:2:3 w image'))

    cmd = 'rm PDBs'
    os.system(cmd)

print "Generating surface contact map"
cmd = 'ls Ads*.pdb > PDBs'
os.system(cmd)
MAIN()

print "Surface contact map DONE"
