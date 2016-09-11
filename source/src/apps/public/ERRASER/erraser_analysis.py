#!/usr/bin/env phenix.python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')
imp.load_source('measure_params', file_path + '/measure_params.py')

from measure_params import compute_torsion
from erraser_util import *

pdb1 = abspath( sys.argv[1] )
pdb2 = abspath( sys.argv[2] )

subprocess_call( 'phenix.clashscore %s > start.clash' % pdb1 )
subprocess_call( 'phenix.clashscore %s > ERRASER.clash' % pdb2 )

clash1 = 0
clash2 = 0
for line in open('start.clash') :
    if "clashscore =" in line :
        clash1 = float( line.split() [-1] )
        break

for line in open('ERRASER.clash') :
    if "clashscore =" in line :
        clash2 = float( line.split() [-1] )
        break

START_rna_validate_data = phenix_rna_validate( pdb1 )
pucker1 = START_rna_validate_data['pucker']
bond1   = START_rna_validate_data['bond']
angle1  = START_rna_validate_data['angle']
suite1  = START_rna_validate_data['suite']

ERRASER_rna_validate_data = phenix_rna_validate( pdb2 )
pucker2 = ERRASER_rna_validate_data['pucker']
bond2   = ERRASER_rna_validate_data['bond']
angle2  = ERRASER_rna_validate_data['angle']
suite2  = ERRASER_rna_validate_data['suite']

suite_outlier1 = sum( ['!!' in suite[3] for suite in suite1] )
suite_outlier2 = sum( ['!!' in suite[3] for suite in suite2] )

print "Geometric Changes Introduced by ERRASER\n"
print "====Overall Statistics=================================="
print "                       Start         ERRASER            "
print "Clashscore           %7.2f         %7.2f" % ( clash1, clash2 )
print "Suite  Outlier       %7d         %7d" % ( suite_outlier1, suite_outlier2 )
print "Pucker Outlier       %7d         %7d" % ( len(pucker1), len(pucker2) )
print "Bond   Outlier       %7d         %7d" % ( len(bond1), len(bond2) )
print "Angle  Outlier       %7d         %7d" % ( len(angle1), len(angle2) )
print "\n"
#################################
print "====Suite Changes======================================="

lines = []
for i, j in zip(suite1, suite2) :
    if i[:3] == j[:3] and i[3] != j[3] :
        res = '%4s%2s%4s ' % tuple(i[:3])
        lines.append( res + '               ' + i[3] + '              ' + j[3] )
if lines != [] :
    print "'!!' stands for an outlier suite.\n"
    print "    Residue            Start         ERRASER            "
    for line in lines :
        print line
else :
    print 'None'
print "\n"

print "====Remaining Outlier Suites============================"
if suite_outlier2 != 0 :
    print "'!!' stands for an outlier suite.\n"
    print "    Residue            Start         ERRASER            "
    for i, j in zip(suite1, suite2) :
        if j[3] == '!!' :
            res = '%4s%2s%4s ' % tuple(i[:3])
            print res, '             ', i[3], '            ', j[3]
else :
    print 'None'
print "\n"

#################################
print "====Pucker Outliers====================================="
if len(pucker1) + len(pucker2) != 0 :
    print "Format: (delta torsion)/(!! or OK).\n"
    print "    Residue            Start         ERRASER            "
    for line1 in pucker1 :
        res1 = '%4s%2s%4s ' % tuple(line1[:3])
        delta1 = float(line1[3])
        info1 = '%8.1f/!!' % delta1
        info2 = '         OK'
        for line2 in pucker2 :
            res2 = '%4s%2s%4s ' % tuple(line2[:3])
            if res2 == res1 :
                delta2 = float(line2[3])
                info2 = '%8.1f/!!' % delta2
                pucker2.remove(line2)
                break
        print res1, '     %s' % info1, '    %s' % info2

    for line2 in pucker2 :
        res2 = '%4s%2s%4s ' % tuple(line2[:3])
        delta2 = float(line2[3])
        info2 = '%8.1f/!!' % delta2
        info1 = '         OK'
        print res2, '     %s' % info1, '    %s' % info2
else :
    print 'None'
print '\n'
#################################
temp = []
for i in bond1 :
    line_split = i
    res = '%4s%2s%4s ' % tuple(line_split[:3])
    atoms = line_split [3] + '-' + line_split[4]
    atoms = atoms.replace(' ', '')
    sigma = float(line_split[5])
    temp.append( [res, atoms, sigma] )
bond1 = temp[:]

temp = []
for i in bond2 :
    line_split = i
    res = '%4s%2s%4s ' % tuple(line_split[:3])
    atoms = line_split [3] + '-' + line_split[4]
    atoms = atoms.replace(' ', '')
    sigma = float(line_split[5])
    temp.append( [res, atoms, sigma] )
bond2 = temp[:]

print "====Bond Length Outliers================================="
if len(bond1) + len(bond2) != 0 :
    print "The values are the number of standard deviation away from ideal value.\n"
    print "    Residue            Bond       Start     ERRASER      "
    for i in bond1 :
        info2 = 'OK'
        for j in bond2 :
            if i[0] == j[0] and i[1] == j[1] :
                info2 = '%.2f/!!' % j[2]
                bond2.remove(j)
                break
        print i[0], '%15s' % i[1], '%8.2f/!!' % i[2], '%11s' % info2
    for i in bond2 :
        info1 = 'OK'
        print i[0], '%15s' % i[1], '%11s' % info1, '%8.2f/!!' % i[2]
else :
    print 'None'
print '\n'
#################################
temp = []
for i in angle1 :
    line_split = i
    res = '%4s%2s%4s ' % tuple(line_split[:3])
    atoms = line_split [3] + '-' + line_split[4] + '-' + line_split[5]
    atoms = atoms.replace(' ', '')
    sigma = float(line_split[6])  
    temp.append( [res, atoms, sigma] )
angle1 = temp[:]

temp = []
for i in angle2 :
    line_split = i
    res = '%4s%2s%4s ' % tuple(line_split[:3])
    atoms = line_split [3] + '-' + line_split[4] + '-' + line_split[5]
    atoms = atoms.replace(' ', '')
    sigma = float(line_split[6])
    temp.append( [res, atoms, sigma] )
angle2 = temp[:]

print "====Bond Angle Outliers==================================="
if len(angle1) + len(angle2) != 0 :
    print "The values are the number of standard deviation away from ideal value.\n"
    print "    Residue           Angle       Start     ERRASER      "
    for i in angle1 :
        info2 = 'OK'
        for j in angle2 :
            if i[0] == j[0] and i[1] == j[1] :
                info2 = '%.2f/!!' % j[2]
                angle2.remove(j)
                break
        print i[0], '%15s' % i[1], '%8.2f/!!' % i[2], '%11s' % info2
    for i in angle2 :
        info1 = 'OK'
        print i[0], '%15s' % i[1], '%11s' % info1, '%8.2f/!!' % i[2]
else :
    print 'None'
print '\n'
#################################

####################################
def find_chi_angle_std_pdb( input_pdb ) :

    def find_coord( coords_list, atm_name ) :
        for i in coords_list :
            if atm_name == i[0] :
                return i[1]
        return []

    data_list = []
    current_res = ''
    res_coords = []
    for line in open(input_pdb) :
        if len(line) > 6 and line[0:4] == 'ATOM' :
            res = line[16:27]
            if current_res != res :
                if current_res != '' :
                    data_list.append( [current_res, res_coords] )
                res_coords = []
                current_res = res
            atm_name = line[12:16].replace(' ', '')
            coord = []
            coord.append( float( line[30:38] ) )
            coord.append( float( line[38:46] ) )
            coord.append( float( line[46:54] ) )
            res_coords.append( [atm_name, coord] )

    chi_list = []
    purines = ['   A', ' ADE', ' A  ', '  rA', '   G', ' GUA', ' G  ', '  rG']
    pyrimidines = ['   U', ' URI', ' U  ', '  rU', '   C', ' CYT', ' C  ', '  rC']
    for res_coords in data_list :
        res_name = res_coords[0][0:4]
        if res_name in purines :
            atom1 = find_coord(res_coords[1], "O4'")
            atom2 = find_coord(res_coords[1], "C1'")
            atom3 = find_coord(res_coords[1], "N9")
            atom4 = find_coord(res_coords[1], "C4")
        elif res_name in pyrimidines :
            atom1 = find_coord(res_coords[1], "O4'")
            atom2 = find_coord(res_coords[1], "C1'")
            atom3 = find_coord(res_coords[1], "N1")
            atom4 = find_coord(res_coords[1], "C2")
        else :
            continue
        if atom1 == [] or atom2 == [] or atom3 == [] or atom4 == [] :
            continue
        chi_list.append( [res_coords[0], compute_torsion(atom1, atom2, atom3, atom4)] )
    return chi_list

def syn_anti_check( chi ) :
    if chi <= 140 and chi > -40 :
        return " syn"
    else :
        return "anti"
####################################

chi_list1 = find_chi_angle_std_pdb( pdb1 )
chi_list2 = find_chi_angle_std_pdb( pdb2 )
print "====Glycosidic Torsion Changes=========================="
lines = []
for i, j in zip(chi_list1, chi_list2) :
    if i[0] == j[0] :
        diff = abs(i[1] - j[1])
        if diff > 180 :
            diff = 360 - diff
        if diff > 90 :
            lines.append( i[0] + '      %6.1f/%s' % (i[1], syn_anti_check(i[1])) + '     %6.1f/%s' % (j[1], syn_anti_check(j[1])) )
if lines != [] :
    print "Format: (chi torsion)/(syn-anti annotation)            "
    print "Chi intervals for syn: (-40, 140]. anti: (-180, -40] and (140, 180].\n"
    print "    Residue            Start         ERRASER            "
    for line in lines :
        print line
else :
    print 'None'
print '\n'
#################################

print "========================================================\n"

remove( 'start.clash' )
remove( 'ERRASER.clash' )
