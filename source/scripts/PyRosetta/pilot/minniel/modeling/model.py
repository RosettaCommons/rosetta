#!/usr/bin/env python
from os import system,popen
import string
from sys import argv
import math
from math import sin,cos
from xyzMath import *
from lex import *
from rosetta import *
#init(extra_options = "-extra_res_fa ./HZN.params")
#init()
import rosetta.protocols.rigid
import numpy as np

# function for assigning grid points
def assign_grid_points(start, max, step):
	return [i for i in np.arange(start, max + step, step)]

def helix_vector(pose, start_pos, end_pos):
	# average the first 4 CA, and the last 4 CA to draw a vector through the first helix
	helix_start_xyz = numeric.xyzVector_Real(0,0,0)
	for i in range (0, 4):
		seqpos = i + start_pos
		helix_start_xyz = helix_start_xyz + pose.residue(seqpos).xyz("CA")
	helix_start_xyz = helix_start_xyz / 4.

	helix_end_xyz = numeric.xyzVector_Real(0,0,0)
	for i in range (0, 4):
		seqpos = end_pos - 1
		helix_end_xyz = helix_end_xyz + pose.residue(seqpos).xyz("CA")
	helix_end_xyz = helix_end_xyz / 4.

	helix_middle_xyz = (helix_start_xyz + helix_end_xyz) / 2.
	return helix_start_xyz, helix_middle_xyz, helix_end_xyz

def get_grid_values(fine=False, sample_values=[]):
    grid_sizes = [ ]
    if(not fine):
        # R is distance of long helix to the Z-axis
        R_samples = assign_grid_points( -6, 6, 2 )
        grid_sizes.append( len(R_samples) )

        # X is distance the long helix is translated up and down the X-axis
        X_samples = assign_grid_points( -4, 4, 2 )
        grid_sizes.append( len(X_samples) )

        # angX is a twist angle around the X axis (degrees)
        angX_samples = assign_grid_points( -6, -3, 1 )
        grid_sizes.append( len(angX_samples) )

        # angY is a twist angle around the Y axis (degres)
        angY_samples = assign_grid_points( -3, 3, 1 )
        grid_sizes.append( len(angY_samples) )

        # angZ is a twist angle around the Z axis (degres)
        angZ_samples = assign_grid_points( -3, 3, 1 )
        grid_sizes.append( len(angZ_samples) )
    elif len(sample_values) == 5:
        R_samples = assign_grid_points( sample_values[1] - 3, sample_values[1] + 3, 1 )
        grid_sizes.append( len(R_samples) )

        # X is distance the long helix is translated up and down the X-axis
        X_samples = assign_grid_points( sample_values[0] - 3, sample_values[0] + 3, 1 )
        grid_sizes.append( len(X_samples) )

        # angX is a twist angle around the X axis (degrees)
        angX_samples = assign_grid_points( sample_values[2] - 2, sample_values[2] + 2, 0.5 )
        grid_sizes.append( len(angX_samples) )

        # angY is a twist angle around the Y axis (degres)
        angY_samples = assign_grid_points( sample_values[3] - 2, sample_values[3] + 2, 0.5 )
        grid_sizes.append( len(angY_samples) )

        # angZ is a twist angle around the Z axis (degres)
        angZ_samples = assign_grid_points( sample_values[4] - 2, sample_values[4] + 2, 1 )
        grid_sizes.append( len(angZ_samples) )
    else:
        raise ValueError("Cannot create fine grid without sample values")
    grid_values = [R_samples, X_samples, angX_samples, angY_samples, angZ_samples]
    return (grid_values, grid_sizes)


class Model:
	init(extra_options = "-in:auto_setup_metals")
	rosetta.basic.Tracer.super_mute(True)

	def __init__(self, inputfile, runNum, nonlocal_removed=False):
		self.inputfile  		= inputfile
		self.nonlocal_removed           = nonlocal_removed
		self.runNum 			= runNum
		self.pose			= pose_from_pdb(self.inputfile)
                self.ca_contacts                = []

	def remove_nonlocal_res(self):
		residue_count = 2
		pose = self.pose
                pose.conformation().detect_bonds()
		while residue_count < pose.total_residue() + 1:
			print "Checking residue: %d %s" % (residue_count, pose.residue(residue_count).name3())
			pose.residue(residue_count).update_connections_to_other_residue(pose.residue(residue_count - 1) )
			if pose.residue(residue_count).is_protein() and not pose.residue(residue_count).is_polymer_bonded(pose.residue(residue_count - 1)):
				pose.delete_polymer_residue(residue_count)
                                print "Deleted"
			else:
				residue_count += 1
		self.nonlocal_removed = True
		self.pose = pose

	def elongate( self, prepend_number=3, append_number=3 ):
		if(not self.nonlocal_removed):
			self.remove_nonlocal_res()
		pose = self.pose
		frag_pose = pose_from_sequence('A'*3,"fa_standard")

		#First find what the first polymer residue is
		first_polymer_res = 1
		for i in range(1, pose.total_residue()+1):
			if pose.residue(i).is_polymer():
				first_polymer_res = i
				break
		if prepend_number > 0:
			#If our sites have CA in them to start with, we'll have to first make sure we're getting the first polymer
			#residue and not a nonpolymer residue (i.e. the calcium ion) and last polymer residue
			#This will affect all of the hard-coded 1s and pose.total_residue()s (and pose.total_residue()-i)
			#It will also affect the range where we're setting the torsion angles (start loop at first_polymer_res - 1)
			remove_lower_terminus_type_from_conformation_residue(pose.conformation(), first_polymer_res)
			for i in range (0, prepend_number-1):
				pose.prepend_polymer_residue_before_seqpos(frag_pose.residue(2),first_polymer_res,True)
			pose.prepend_polymer_residue_before_seqpos(frag_pose.residue(1),first_polymer_res,True)  #last one to prepend should be n-terminal residue
			for i in range (first_polymer_res-1, prepend_number+first_polymer_res):  #make sure old n-terminus also has helical phi,psi
				pose.set_phi(i+1, -65)
				pose.set_psi(i+1, -45)
				pose.set_omega(i+1,180)

                #Now find what the last polymer residue is
		last_polymer_res=pose.total_residue()
		for i in range(pose.total_residue(), 0, -1):
			if pose.residue(i).is_polymer():
				last_polymer_res = i
				break

		if append_number > 0:
			remove_upper_terminus_type_from_conformation_residue(pose.conformation(), last_polymer_res)
			for i in range (0, append_number-1):
				pose.append_polymer_residue_after_seqpos(frag_pose.residue(2),last_polymer_res,True)
				last_polymer_res += 1
			pose.append_polymer_residue_after_seqpos(frag_pose.residue(3),last_polymer_res,True)  #last one to append should be c-terminal type
			last_polymer_res += 1
			for i in range (0, append_number+1):  #make sure old n-terminus also has helical phi,psi
				pose.set_phi(last_polymer_res-i, -65)
				pose.set_psi(last_polymer_res-i, -45)
				pose.set_omega(last_polymer_res-i,180)

		#minimize the new bits
		movemap = MoveMap()
		for i in range (first_polymer_res, first_polymer_res+prepend_number+1):
			movemap.set_bb(i,True)
		for i in range (last_polymer_res-append_number,last_polymer_res+1):
			movemap.set_bb(i,True)
		scorefxn = get_fa_scorefxn()
		minmover = MinMover(movemap, scorefxn, "linmin", 0.01, True)
		minmover.apply(pose)
		pose.dump_pdb(self.inputfile[:-3] + "elongated.pdb")
		self.inputfile = self.inputfile[:-3] + "elongated.pdb"
		self.pose = pose


	def create_start_models(self, first_helical_residue_1, last_helical_residue_1, first_helical_residue_2, last_helical_residue_2, origin_flag, fine=False, score_file_path=""):
		if(not self.nonlocal_removed):
			self.remove_nonlocal_res()

                pose = self.pose
		inputfile = self.inputfile
		runNum = self.runNum
		origin =  numeric.xyzVector_Real(0,0,0)
		residue_count = 2
		pose.conformation().detect_bonds()

		# start by aligning the first helix with the Y-axis
		first_helix_points = helix_vector(pose, first_helical_residue_1, last_helical_residue_1)
		first_helix_start_xyz = first_helix_points[0]
		first_helix_middle_xyz = first_helix_points[1]
		first_helix_end_xyz = first_helix_points[2]

		# move end of the first helix to origin
		y_axis = numeric.xyzVector_Real(0,1,0)
		z_axis = numeric.xyzVector_Real(0,0,1)
		x_axis = numeric.xyzVector_Real(1,0,0)
		identity_matrix = numeric.rotation_matrix_radians(z_axis, 0)
		if origin_flag == "-l":
			pose.apply_transform_Rx_plus_v(identity_matrix, -first_helix_end_xyz) #just performing translation
		elif origin_flag == "-h":
			pose.apply_transform_Rx_plus_v(identity_matrix, -first_helix_middle_xyz) #just performing translation
		else:
			print "Incorrect Origin Flag"
			sys.exit()

		# rotate the first helix on to the z-axis
		helix_vec = first_helix_end_xyz - first_helix_start_xyz
		angle_to_z_axis = numeric.angle_radians(helix_vec, origin, -z_axis)
		rotation_axis = helix_vec.cross(-z_axis)
		rot_matrix = numeric.rotation_matrix_radians(rotation_axis, angle_to_z_axis)
		pose.apply_transform_Rx_plus_v(rot_matrix, origin)

		# now rotate around the z-axis to point the second helix in the -x direction
		# this will set things up so that in the starting structures the second helix
		# is pointing back towards the other chain
		second_helix_points = helix_vector(pose, first_helical_residue_2, last_helical_residue_2)
		second_helix_start_xyz = second_helix_points[0]
		second_helix_end_xyz = second_helix_points[2]
		helix_vec = second_helix_end_xyz - second_helix_start_xyz
		# angle to x-axis of second helix
		angle_to_x_axis = numeric.angle_radians(helix_vec,origin,-x_axis)
		rot_matrix = numeric.rotation_matrix_radians(z_axis,angle_to_x_axis) #just performing rotation
		pose.apply_transform_Rx_plus_v(rot_matrix, origin)


		if( fine ): # if a fine sampling, parse grid data from file name of top 5
			if(os.path.isfile(score_file_path)):
				with open(score_file_path, 'r') as score_file:
					line_count = 1
					for line in score_file:
						if (line_count >= 6):#stop after top 5
							break
						else:
							line = line.strip()
							score_fields = line.split()
							fileName = score_fields[-1]
							print "Sampling around " + fileName
							fileName_fields = fileName.split('_')
							sample_values = [float(i) for i in fileName_fields[5:10]]
							print sample_values
                                                        print fileName
                                                        grid_values, grid_sizes = get_grid_values(fine, sample_values)
							self.sample_grid(grid_values, grid_sizes, fine) #apply rotations and output pdb
							line_count += 1
		else: #do a coarse sampling with preset grid values
                        grid_values, grid_sizes = get_grid_values()
			self.sample_grid(grid_values, grid_sizes, fine)

	def sample_grid(self, grid_values, grid_sizes, fine):
		#apply rotations and output pdb
                pose = self.pose
		inputfile = self.inputfile
		runNum = self.runNum

		R_samples = grid_values[0]
		X_samples = grid_values[1]
		angX_samples = grid_values[2]
		angY_samples = grid_values[3]
		angZ_samples = grid_values[4]

		y_axis = numeric.xyzVector_Real(0,1,0)
		z_axis = numeric.xyzVector_Real(0,0,1)
		x_axis = numeric.xyzVector_Real(1,0,0)
		identity_matrix = numeric.rotation_matrix_radians(z_axis, 0)
		origin =  numeric.xyzVector_Real(0,0,0)

		start_pose = Pose()
		start_pose.assign(pose)  # save starting pose
		lex = LexicographicalIterator( grid_sizes )   # need lex.py, allows multidimensional grid search
		while not lex.at_end:          # enumerates all degrees of freedom
			pose.assign(start_pose)    #retreive starting pose

			R = R_samples[ lex.pos[0] ]
			X = X_samples[ lex.pos[1] ]
			angX = angX_samples[ lex.pos[2] ]
			angY = angY_samples[ lex.pos[3] ]
			angZ = angZ_samples[ lex.pos[4] ]
			lex.increment()

                        # rotate around Y axis
			rotY = math.radians( angY )
			rot_matrix = numeric.rotation_matrix_radians(y_axis,rotY)
			pose.apply_transform_Rx_plus_v(rot_matrix, origin)

                        # rotate around X axis
			rotX = math.radians( angX )
			rot_matrix = numeric.rotation_matrix_radians(x_axis,rotX)
			pose.apply_transform_Rx_plus_v(rot_matrix, origin)

                        # rotate around Z axis
			rotZ = math.radians( angZ )
			rot_matrix = numeric.rotation_matrix_radians(z_axis,rotZ)
			pose.apply_transform_Rx_plus_v(rot_matrix, origin)

			# set R (distance from the z-axis along the y-axis)
			translateR = numeric.xyzVector_Real(R,0,0)
			pose.apply_transform_Rx_plus_v(identity_matrix,translateR)

			# set x (translate along X axis)
			translateY = numeric.xyzVector_Real(0,X,0)
			pose.apply_transform_Rx_plus_v(identity_matrix,translateY)

                        #check for N & cterm constraints
                        pose2 = Pose()
                        pose2.assign(pose)
                        rot_matrix = numeric.rotation_matrix_radians(z_axis,math.pi)
                        pose2.apply_transform_Rx_plus_v(rot_matrix, origin)
                        dist_C2N = pose.residue(1).xyz('C').distance(pose2.residue(pose2.total_residue()-1).xyz('C'))
                        if(dist_C2N > 16 or dist_C2N < 4):
                            print "Bad Distance: %d" % dist_C2N
                            continue

                        #check for/create file directory for output files
                        if(not os.path.isdir("output")):
                                os.mkdir("output")
                        file_fields = inputfile.split('/')  #file path
                        pdb_file = file_fields[-1]  #file name
                        pdb_fields = pdb_file.split('.')
                        pdb_name = '.'.join(pdb_fields[0:2]) #just the pdb code & site type/number
                        if( not os.path.isdir("output/" + pdb_name)):
                                os.mkdir("output/" + pdb_name)
                        output_path_prefix = "output/" + pdb_name + "/" + pdb_file[:-4]

			# if fine, create symmetric partner & sample rotamers of sidechains near CA
                        if ( fine and not("local" in inputfile)):
                            pose2 = Pose()
                            pose2.assign(pose)
                            rot_matrix = numeric.rotation_matrix_radians(z_axis,math.pi)
                            pose2.apply_transform_Rx_plus_v(rot_matrix, origin)
                            ca_xyz = pose.residue(pose.total_residue()).xyz(1)
                            ca_contacts = []
                            for resNum in range(1, pose.total_residue()):
                                rset = pose.residue(resNum).residue_type_set()
                                if pose2.residue(resNum).name3() == 'GLY':
                                    new_res = rosetta.core.conformation.ResidueFactory.create_residue(rset.name_map('ALA'),pose2.residue(resNum),pose2.conformation())
                                    rosetta.core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose2.residue(resNum),new_res,pose2.conformation())
                                    pose2.replace_residue(resNum,new_res,False)
                                dist_c_beta = pose2.residue(resNum).xyz('CB').distance(ca_xyz)
                                if dist_c_beta < 6.0:
                                    ca_contacts.append(resNum)

                            best_rot_dict = {}
                            print "Checking Residues for CA contacts",
                            for resNum in ca_contacts: #ignoring CA which is the last residue in pose
                                print "#",
                                min_ls_value = 2.42 # the max least squares (from 3.5 to 2.4) to be a valid CA bond 
                                #keep track of mutated residue
                                current_res_name = pose2.residue(resNum).name3()
                                if pose.residue(resNum).is_bonded(pose.residue(pose.total_residue())):
                                    continue
                                rset = pose.residue(resNum).residue_type_set()
                                for rotamer in self.get_rotamer():
                                    #mutate residue to returned rotamer
                                    if rotamer[0] != current_res_name:
                                        new_res = rosetta.core.conformation.ResidueFactory.create_residue(rset.name_map(rotamer[0]),pose2.residue(resNum),pose2.conformation())
                                        rosetta.core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose2.residue(resNum),new_res,pose2.conformation())
                                        pose2.replace_residue(resNum,new_res,False)
                                        current_res_name = rotamer[0]
                                    #rotate to returned rotamer
                                    pose2.set_chi(1,resNum,rotamer[1])
                                    pose2.set_chi(2,resNum,rotamer[2])
                                    if( rotamer[0] == 'ASP' ):
                                        oxy_xyz = (pose2.residue(resNum).xyz('OD1') + pose2.residue(resNum).xyz('OD2')) / 2
                                        dist_1 = pose2.residue(resNum).xyz('OD1').distance(ca_xyz)
                                        dist_2 = pose2.residue(resNum).xyz('OD2').distance(ca_xyz)
                                        vertex_xyz = pose2.residue(resNum).xyz('CG')
                                    else:
                                        pose2.set_chi(3,resNum,rotamer[3]) #if glu, set chi3 as well.
                                        oxy_xyz = (pose2.residue(resNum).xyz('OE1') + pose2.residue(resNum).xyz('OE2')) / 2
                                        dist_1 = pose2.residue(resNum).xyz('OE1').distance(ca_xyz)
                                        dist_2 = pose2.residue(resNum).xyz('OE2').distance(ca_xyz)
                                        vertex_xyz = pose2.residue(resNum).xyz('CD')
                                    if not (dist_1 < 3.4 and dist_2 < 3.4): #CA bonds < 3.5 Ang
                                        continue
                                    angle_to_ca = numeric.angle_degrees(ca_xyz, vertex_xyz, oxy_xyz)
                                    if not ( angle_to_ca < 45 or angle_to_ca > 315):
                                        print "\nBad Angle: %d" % angle_to_ca
                                        continue
                                    rot_least_sq = (dist_1 - 2.4)**2 + (dist_2 - 2.4)**2
                                    if( rot_least_sq < min_ls_value ):
                                        print "\nAdding Rotamer to dictionary"
                                        best_rot_dict[resNum] = rotamer
                                        min_ls_value = rot_least_sq
                            print ""
                            #output all new pdbs
                            for resNum, rotamer in best_rot_dict.iteritems():
                                output_pose = Pose()
                                output_pose.assign(pose)
                                new_res = rosetta.core.conformation.ResidueFactory.create_residue(rset.name_map(rotamer[0]),output_pose.residue(resNum),output_pose.conformation())
                                rosetta.core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(output_pose.residue(resNum),new_res, output_pose.conformation())
                                output_pose.replace_residue(resNum,new_res,False)
                                #rotate to returned rotamer
                                output_pose.set_chi(1,resNum,rotamer[1])
                                output_pose.set_chi(2,resNum,rotamer[2])
                                if( rotamer[0] == 'GLU'):
                                    output_pose.set_chi(3,resNum,rotamer[3])
                                tag = '%s_run%d_fine_%d_%s'%(output_path_prefix, runNum, resNum, rotamer[0])
                                out_file_name='%s_%.2f_%.2f_%.2f_%.2f_%.2f.%s'%(tag,X,R,angX,angY,angZ,'pdb')
                                #for run 1 output was %(tag,R,X,angX,angY,angZ,'pdb')
                                output_pose.dump_pdb(out_file_name)

                        else: #output pdbs
                            if(fine): #ca contacts are all local
                                tag = '%s_run%d_fine'%(output_path_prefix, runNum)
                            else: #coarse output
                                tag = '%s_run%d'%(output_path_prefix, runNum)
                            out_file_name='%s_%.2f_%.2f_%.2f_%.2f_%.2f.%s'%(tag,X,R,angX,angY,angZ,'pdb')
                            #for run 1 output was %(tag,R,X,angX,angY,angZ,'pdb')
                            pose.dump_pdb(out_file_name)
		self.pose = start_pose #reset pose for next iteration

	def get_rotamer(self):
            #sample Aspartate rotamers
            for chi1 in  np.arange(-82.5, -62.5, 2):
                for chi2 in np.arange(-22.5, -2.5, 5):
                    rotamer = ['ASP', chi1, chi2]
                    yield rotamer
          
            #sample Glutamate rotamers
            for chi1 in  np.arange(-185, -165, 2):
                for chi2 in np.arange(165, 185, 5):
                    for chi3 in np.arange(-20, 0, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
                    for chi3 in np.arange(-50, -30, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
            for chi1 in  np.arange(-80, -60, 2):
                for chi2 in np.arange(165, 185, 5):
                    for chi3 in np.arange(-20, 0, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
                    for chi3 in np.arange(-50, -30, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
                for chi2 in np.arange(-70, -50, 5):
                    for chi3 in np.arange(-20, 0, 2):
                       rotamer = ['GLU', chi1, chi2, chi3]
                       yield rotamer
                    for chi3 in np.arange(-50, -30, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
            for chi1 in  np.arange(-187, -167, 2):
                for chi2 in np.arange(-167, -187, 5):
                    for chi3 in np.arange(-20, 0, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
                    for chi3 in np.arange(-50, -30, 2):
                        rotamer = ['GLU', chi1, chi2, chi3]
                        yield rotamer
