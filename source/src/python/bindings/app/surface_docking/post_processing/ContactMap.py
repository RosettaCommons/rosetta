#!/usr/bin/env python

# Edited by Emily Koo to create difference plot if native PDB is included in dir

from numpy import zeros, arange
from sys import argv
from PDBInfo import *
from IO import load
import os
from os.path import exists

def ContactMap(Con, native, num_pdb):

    X = 0
    Y = 1


    total_res = len(ResidueList(load('PDBsInDir')[0])) + 1
    Resolution=[total_res, total_res]
    Histogram = zeros(Resolution)

    Points = []; Points2 = []

    for Line in open(Con):
        Contact = [int(Line.split()[2]), int(Line.split()[5])]
        Points.append(Contact)

    start_res = (Points[0])[0] - 1

    for Line in Points:
        Res1 = int(Line[X]) - start_res
        Res2 = int(Line[Y]) - start_res
        Histogram[Res1,Res2] += 1

    if exists(native):
        for Line in open(native, 'r'):
            Contact = [int(Line.split()[2]), int(Line.split()[5])]
            Points2.append(Contact)

        for Line in Points2:
            Res1 = int(Line[X]) - start_res
            Res2 = int(Line[Y]) - start_res
            Histogram[Res1,Res2] -= num_pdb

    Index = [arange(1,total_res,1),arange(1,total_res,1)]

    FileOut = 'ContactMap.dat'
    o = open(FileOut, 'w')
    for x in enumerate(Index[X]):
        for y in enumerate(Index[Y]):
            o.write('%1i %1s %1i %1s %1i \n' %(x[1],' ', y[1], ' ', Histogram[x[1],y[1]]))

    # Generate GNU script to make the plot
    x1=0.8
    y2=1
    x2=y1=total_res+0.5
    Conversion=(total_res-1)+(total_res-1)/50 + 1.0

    GNUPlotFile = 'ContactMap.gnuplot'
    O = open(GNUPlotFile, 'w')
    O.write('%0s \n' %('set view map'))
    O.write('%0s \n' %('set terminal png size 800, 800 nocrop truecolor enhanced font "VeraBd, 15"'))
    O.write('%0s \n' %('set output "ContactMap.png"'))
    O.write('%0s \n' %('set size 1,1.1'))
    O.write('%0s %1.1f %0s \n' %('set xrange [0.5:',y1-1.0,']'))
    O.write('%0s %1.1f %0s \n' %('set yrange [0.5:',y1-1.0,']'))
    O.write('%0s \n' %('unset key'))

    if exists(native):
        O.write('%0s \n' %('set palette defined (0 "#FF0000", 1 "#FFFFFF", 2 "#0000FF")'))
        #O.write('%0s \n' %('set palette defined (0 "#FF0000", 1 "#0000FF")'))
        O.write('%0s \n' %('set cbrange [-' + str(num_pdb) + ':' + str(num_pdb) + ']'))
        O.write('%0s \n' %('set cbtics (-' + str(num_pdb) + ', 0, ' + str(num_pdb) + ') nomirror scale 0'))
    else:
        O.write('%0s \n' %('set palette defined (0 "#FFFFFF", 1 "#CCFFFF", ' + str(num_pdb) + ' "#0000FF")'))
        O.write('%0s \n' %('set cbrange [0:' + str(num_pdb) + ']'))
        O.write('%0s \n' %('set cbtics (0, ' + str(num_pdb/4) + ', ' + str(num_pdb/2) + ', ' + str(3*num_pdb/4) + ', ' + str(num_pdb) + ') nomirror scale 0'))

    O.write('%0s \n' %('set colorbox user origin 0.875, 0.1375 size 0.025, 0.785'))
    O.write('%0s \n' %('set xtic out nomirror'))
    O.write('%0s \n' %('set ytic out nomirror'))
    O.write('%0s \n' %('set mxtics'))
    O.write('%0s \n' %('set mytics'))
    #O.write('%0s \n' %('set cblabel "Contact Frequency" -2,0 font "VeraBd, 15"'))
    O.write('%0s \n' %('set cblabel "Contact Frequency" font "VeraBd, 15"'))
    O.write('%0s \n' %('set border lw 3'))

    for Atom in protein:
        if Atom.AtmTyp == 'CA':
            # Gamma symbol in gnuplot
            if Atom.ResTyp == 'Z':
                output_res = "{/Symbol g}"
            else:
                output_res = AminoAcidTable(Atom.ResTyp)

            O.write('%0s %1.1f %0s %1.1f %0s \n' %("set label " '"'+output_res+'"' " at ",x1,",",Conversion, ' font "VeraBd, 9"'  ))
            O.write('%0s %1.1f %0s %1.1f %0s \n' %("set label " '"'+output_res+'"' " at ",Conversion,",",y2, ' font "VeraBd, 9"'  ))
            x1+=1
            y2+=1


    O.write('%0s \n' %('set xlabel "Residue" font "VeraBd, 20"'))
    O.write('%0s \n' %('set ylabel "Residue" font "VeraBd, 20"'))
    O.write('%0s \n' %('set bmargin 0'))
    #O.write('%0s \n' %('set title "Osteocalcin/HAp [100] (Solution)" font "VeraBd, 20"'))
    O.write('%0s \n' %('splot "ContactMap.dat" u 1:2:3 w image'))

Con = argv[1]
native = argv[2]
num_pdb = 0
for files in open('PDBsInDir', 'r'):
    if files.startswith('Ads'):
        num_pdb += 1

ContactMap(Con, native, num_pdb)
