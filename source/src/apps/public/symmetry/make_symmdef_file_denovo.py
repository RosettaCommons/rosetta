#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import string
import sys
from sys import argv,stdout
from os.path import exists

length = -1
if argv.count('-symm_type') and argv.count('-nsub'):
	pos = argv.index('-symm_type')
	symm_type_name = argv[pos+1]
	pos = argv.index('-nsub')
	nsub = argv[pos+1]
else:
	print ""
	print "usage: make_symmdef_file_denovo.py [options] "
	print ""
	print "example: make_symmdef_file_denovo.py -symm_type cn -nsub 2"
	print "example: make_symmdef_file_denovo.py -symm_type dn -nsub 4"
	print "example: make_symmdef_file_denovo.py -symm_type dn -nsub 4 -slide_type RANDOM -slide_criteria_type CEN_DOCK_SCORE"
	print "example: make_symmdef_file_denovo.py -symm_type cn -nsub 24 -subsystem"
	print ""
	print "common options:"
	print	"	-symm_type (cn|dn)                                          : The type of symmetry. Currently cyclic or dihedral symmetry"
	print "	-nsub <integer>                                             : number of subunits"
	print "	-subsystem                                                  : Simulate a smallers subsystem. For cn 3 consecutive subunits are represented."
	print "	                                                              For dn 6 subunits are represented 3 in the upper ring and 3 in the lower ring"
	print "	-slide_type (RANDOM|SEQUENTIAL|ORDERED_SEQUENTIAL)          : Defines how a multidimensional slide should be performed."
	print "	                                                              For dn symmetry there are two sliding directions. A slide "
	print "	                                                              can be done by randomly selecting a slide direction for each"
	print "	                                                              slide step (RANDOM), randomly deciding on which direction should"
	print "	                                                              be slided first but always sequentially go through both (SEQUENTIAL),"
	print "	                                                              or define the order yourself (ORDERED_SEQUENTIAL). RANDOM default"
	print "	-slide_criteria_type (CEN_DOCK_SCORE|FA_REP_SCORE|CONTACTS) : Defines what the criteria is for abandaning a slide. Either"
	print "	                                                              the CEN_DOCK_SCORE, FA_REP_SCORE or the number of contacts. CEN_DOCK_SCORE deafult"
	print "	-slide_criteria_val <string|float>                          : Sets the actual value when a slide move is abandoned given the criteria type. By default"
	print "	                                                              set to AUTOMATIC, which means that ROSETTA figures it out by itself. Really"
	print "	                                                              only useful for the CONTACTS type"
	sys.exit()

# Check additional arguments
#initialize to empty values

subsystem_size=0
slide_type=""
slide_criteria_type=""
slide_order=""
slide_criteria_val=""

if argv.count('-slide_type'):
	pos = argv.index('-slide_type')
	slide_type=argv[pos+1]
	if slide_type != "RANDOM" and slide_type != "ORDERED" and slide_type != "ORDERED_SEQUENTIAL":
		print "ERROR: Unknown slide_type " + slide_type
		sys.exit()
if argv.count('-slide_criteria_type'):
	pos = argv.index('-slide_criteria_type')
	slide_criteria_type=argv[pos+1]
	if slide_criteria_type != "CEN_DOCK_SCORE" != slide_criteria_type != "FA_REP_SCORE" and slide_criteria_type != "CONTACTS":
		print "ERROR: Unknown slide_criteria_type " + slide_criteria_type
		sys.exit()
if argv.count('-slide_criteria_val'):
	pos = argv.index('-slide_criteria_type')
	slide_criteria_val=argv[pos+1]

subunits=int(nsub)

# Store the interfaces in this list. The element order specifies
# the interface identity (first, second etc) and the value in the
# list specifies the number of such interfaces
interfaces = []

# Divide it into two cases: subsystem or no subsystem
if not argv.count('-subsystem') or \
	argv.count('-subsystem') and  symm_type_name == "cn" and subunits < 4 \
	or argv.count('-subsystem') and  symm_type_name == "dn" and subunits < 8:
	if symm_type_name == "cn":
		if subunits == 2:
			interfaces.append(1)
		if subunits == 3:
			interfaces.append(3)
		if subunits == 4:
			interfaces.append(4)
			interfaces.append(2)
		if subunits == 5:
			interfaces.append(5)
			interfaces.append(5)
		if subunits == 6:
			interfaces.append(6)
			interfaces.append(6)
			interfaces.append(3)
	#	For large rings assume the subunit 1 and subunit >4 do
	#	not interact
		if subunits > 6:
			interfaces.append(subunits)
			interfaces.append(subunits)
			interfaces.append(subunits)

	# Now define dn symmetry. An interface value of 0 is ignored
	if symm_type_name == "dn":
		if subunits == 4:
			interfaces.append(2)
			interfaces.append(2)
			interfaces.append(2)
		if subunits == 6:
			interfaces.append(6)
			interfaces.append(0)
			interfaces.append(3)
			interfaces.append(3)
			interfaces.append(3)
		if subunits == 8:
			interfaces.append(8)
			interfaces.append(4)
			interfaces.append(0)
			interfaces.append(4)
			interfaces.append(4)
			interfaces.append(4)
			interfaces.append(4)
		if subunits == 10:
			interfaces.append(10)
			interfaces.append(10)
			interfaces.append(0)
			interfaces.append(0)
			interfaces.append(5)
			interfaces.append(5)
			interfaces.append(5)
			interfaces.append(5)
			interfaces.append(5)
		if subunits == 12:
			interfaces.append(12)
			interfaces.append(12)
			interfaces.append(6)
			interfaces.append(0)
			interfaces.append(0)
			interfaces.append(6)
			interfaces.append(6)
			interfaces.append(6)
			interfaces.append(6)
			interfaces.append(6)
			interfaces.append(6)
		if subunits > 12:
			interfaces.append(subunits)
			interfaces.append(subunits)
			interfaces.append(subunits)
			for num in range(0,subunits/2-4):
				interfaces.append(0)
			for num in range(0,subunits/2):
				interfaces.append(subunits/2)
		if subunits < 4:
			print "ERROR: the smallest complex with dihedral symmetry is 4!"
			sys.exit()
		if subunits%2 != 0:
			print "ERROR: there must be an even number of subunits for complexes with dihedral symmetry!"
			sys.exit()
	# Find out have many different types of interfaces we have
	num_interfaces=0
	for num in interfaces:
		if num !=0:
			num_interfaces=num_interfaces+1

	# Print the symmetry definition to std out
	if symm_type_name == "cn":
		print "symmetry_name " + symm_type_name[:1] + nsub
	if symm_type_name == "dn":
		 print "symmetry_name " + symm_type_name[:1] + str(subunits/2)
	print "subunits " + nsub
	print "recenter"
	print "number_of_interfaces ", num_interfaces

	# Print the E line using the information stored in interfaces
	E_string = "E = "+ nsub + "*VRT0001"
	interface_num=0
	for interface in interfaces:
		interface_num=interface_num+1
	# disregard interfaces with values == 0
		if interface !=0:
			E_string = E_string + " + " + str(interface) + "*(VRT0001:VRT" "%04d"%(interface_num+1) + ")"
	print  E_string
	# default anchor position at com
	print "anchor_residue COM"
	# define the transforms
	print "virtual_transforms_start"

	# Two cases, cn and dn
	if symm_type_name == "cn":
		print "start -1,0,0 0,1,0 0,0,0"
		print "rot Rz", nsub
		print "virtual_transforms_stop"

	# Define the jumps in the system
		for sub in range(subunits-1):
			jump_num=sub+1
			vrt_start_num=sub+1
			vrt_stop_num=sub+2
			jump_start='VRT' + str(vrt_start_num).zfill(4)
			jump_stop='VRT' + str(vrt_stop_num).zfill(4)
			print 'connect_virtual JUMP' + str(jump_num) + ' ' + jump_start + ' ' + jump_stop
	# Set the dof line. Use default values. These have to be altered by hand
		print "set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)"
		#print "slide_type ORDERED_SEQUENTIAL"
		#print "slide_order BASEJUMP"

	if symm_type_name == "dn":
	# Define the z cyclical symmetry rotation
		rotation_z=float(2*360.000/(float(subunits)))
		print "start -1,0,0 0,1,0 0,0,0"
	# First make the upper ring
		for num in range(1,subunits/2):
			print "rot Rz_angle", rotation_z
	# A twofold rotation to make the second ring
		print "rot Rx_angle 180.0"
	# And now the second ring
		for num in range(1,subunits/2):
			print "rot Rz_angle", rotation_z
		print "virtual_transforms_stop"
		for sub in range(subunits-1):
			jump_num=sub+1
			vrt_start_num=sub+1
			vrt_stop_num=sub+2
			jump_start='VRT' + str(vrt_start_num).zfill(4)
			jump_stop='VRT' + str(vrt_stop_num).zfill(4)
			print 'connect_virtual JUMP' + str(jump_num) + ' ' + jump_start + ' ' + jump_stop
	# Define the dof line, default values. Have to change manually.
		print "set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)"
	#Define the jump that relates one ring to the other along z
		z_jump=subunits/2
		z_random_rotation=float(360.000/float(subunits))
		print "set_dof JUMP" + str(z_jump) + ' z(50) angle_z(0:' + str(z_random_rotation) + ')'
	# Print the slide information if flags are given. Otherwise default values will be used by ROSETTA
		if slide_type != "":
			print "slide_type", slide_type
		if slide_criteria_type != "":
			print "slide_criteria_type", slide_criteria_type
		if slide_criteria_val != "":
			print "slide_criteria_val", slide_criteria_val
	# Slide order is printed in a default order. User have to manually change the order if he desires.
		if argv.count('-slide_order'):
	# This only make sense if ORDERED_SEQUENCE is used.
			if slide_type == "ORDERED_SEQUENTIAL":
				print "slide_order BASEJUMP" + ' JUMP' + str(z_jump)
# subsystem case. We cover this special case in a hard coded manner. Only one type of setup.
else:
	# Print the symmetry definition to std out
	if symm_type_name=="cn":
		print "symmetry_name " + symm_type_name[:1] + nsub
		print "subunits 3"
		print "recenter"
		print "number_of_interfaces 1"
		# Print the E line using the information stored in interfaces
		E_string = "E = "+ nsub + "*VRT0001"
		interface_num=0
		E_string = E_string + " + " + nsub + "*(VRT0001:VRT0002)"
		print  E_string

	elif symm_type_name=="dn":
		print "symmetry_name " + symm_type_name[:1] + str(subunits/2)
		print "subunits 6"
		print "number_of_interfaces 1"
		print "recenter"
		# Print the E line using the information stored in interfaces
		E_string = "E = "+ nsub + "*VRT0001"
		# disregard interfaces with values == 0
		E_string = E_string + " + " + nsub + "*(VRT0001:VRT0002) + " + str(subunits/2) + "*(VRT0001:VRT0004) + " + str(subunits/2) + "*(VRT0001:VRT0005) + " + str(subunits/2) + "*(VRT0001:VRT0006)"
		print  E_string
	# default anchor position at com
	print "anchor_residue COM"
	# define the transforms
	print "virtual_transforms_start"

	# Two cases, cn and dn
	if symm_type_name == "cn":
		rot_z1=float(360.000/float(subunits))
		rot_z2=360.000 -2*rot_z1
		print "start -1,0,0 0,1,0 0,0,0"
		print "rot Rz_angle", rot_z1
		print "rot Rz_angle", rot_z2
		print "virtual_transforms_stop"

	# Define the jumps in the system
		print 'connect_virtual JUMP1 VRT0001 VRT0002'
		print 'connect_virtual JUMP2 VRT0002 VRT0003'
	# Set the dof line. Use default values. These have to be altered by hand
		print "set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)"
		#print "slide_type ORDERED_SEQUENTIAL"
		#print "slide_order BASEJUMP"

	if symm_type_name == "dn":
	# Define the z cyclical symmetry rotation
		rotation_z=float(2*360.000/float(subunits))
		rotation_z2= 360.000 - 2*rotation_z
		print "start -1,0,0 0,1,0 0,0,0"
	# First make the upper ring
		print "rot Rz_angle", rotation_z
		print "rot Rz_angle", rotation_z2
	# A twofold rotation to make the second ring
		print "rot Rx_angle 180.0"
	# And now the second ring
		print "rot Rz_angle", rotation_z
		print "rot Rz_angle", rotation_z
		print "virtual_transforms_stop"
		print 'connect_virtual JUMP1 VRT0001 VRT0002'
		print 'connect_virtual JUMP2 VRT0002 VRT0003'
		print 'connect_virtual JUMP3 VRT0003 VRT0004'
		print 'connect_virtual JUMP4 VRT0004 VRT0005'
		print 'connect_virtual JUMP5 VRT0005 VRT0006'
	# Define the dof line, default values. Have to change manually.
		print "set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)"
	#Define the jump that relates one ring to the other along z
		z_jump=subunits/2
		z_random_rotation=float(360.000/float(subunits))
		print "set_dof JUMP3" + ' z(50) angle_z(0:' + str(z_random_rotation) + ')'
	# Print the slide information if flags are given. Otherwise default values will be used by ROSETTA
		if slide_type != "":
			print "slide_type", slide_type
		if slide_criteria_type != "":
			print "slide_criteria_type", slide_criteria_type
		if slide_criteria_val != "":
			print "slide_criteria_val", slide_criteria_val
	# Slide order is printed in a default order. User have to manually change the order if he desires.
		if argv.count('-slide_order'):
	# This only make sense if ORDERED_SEQUENCE is used.
			if slide_type == "ORDERED_SEQUENTIAL":
				print "slide_order BASEJUMP" + ' JUMP' + str(z_jump)
