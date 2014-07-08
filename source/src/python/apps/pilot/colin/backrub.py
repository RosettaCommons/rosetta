#!/usr/bin/env python

import optparse
import sys

parser = optparse.OptionParser()

parser.set_usage("%prog [options] [-- rosetta_options]")
parser.add_option("--ntrials", action="store", type="int", default=1000,
                  help="number of Monte Carlo trials to run")
parser.add_option("--sc_prob", action="store", type="float", default=0.25,
                  help="probability of making a side chain move")
parser.add_option("--sm_prob", action="store", type="float", default=0.,
                  help="probability of making a small move")
parser.add_option("--sc_prob_uniform", action="store", type="float", default=0.1,
                  help="probability of uniformly sampling chi angles")
parser.add_option("--mc_kt", action="store", type="float", default=0.6,
                  help="value of kT for Monte Carlo")
parser.add_option("--mm_bend_weight", action="store", type="float", default=1.0,
                  help="weight of mm_bend bond angle energy term")
parser.add_option("--initial_pack", action="store_true", default=False,
                  help="force a repack at the beginning regardless of whether mutations are set in the resfile")

(options, rosetta_args) = parser.parse_args()

if "--" in sys.argv:
	rosetta_args = [" ".join(sys.argv[:sys.argv.index("--")+1])] + rosetta_args
else:
	rosetta_args = [sys.argv[0]] + rosetta_args

import time
import_time_start = time.time()
print "Importing rosetta module..."

from rosetta import core, protocols
import rosetta.core.mm
import rosetta.core.pack.task.operation
import rosetta.core.scoring.methods
# *** Error while importing rosetta.core.scoring.constraints ***
#import rosetta.core.scoring.constraints
import rosetta.protocols.branch_angle
import rosetta.protocols.jobdist

import_time_end = time.time()
print "Time to import rosetta module: %.2f seconds" % (import_time_end - import_time_start)

rosetta.init(*rosetta_args)

# create a TaskFactory with the resfile
main_task_factory = core.pack.task.TaskFactory()
main_task_factory.push_back(core.pack.task.operation.InitializeFromCommandline())
if core.get_file_vector_option("packing:resfile")[1] != "resfile":
	main_task_factory.push_back(core.pack.task.operation.ReadResfile())
else:
	main_task_factory.push_back(core.pack.task.operation.RestrictToRepacking())
# C-beta atoms should not be altered during packing because branching atoms are optimized
main_task_factory.push_back(core.pack.task.operation.PreserveCBeta())

# set up the score function and add the bond angle energy term
score_fxn = core.scoring.get_score_function()
score_fxn.set_weight(core.scoring.mm_bend, options.mm_bend_weight)
energymethodoptions = score_fxn.energy_method_options()
energymethodoptions.decompose_bb_hb_into_pair_energies(True)
energymethodoptions.bond_angle_central_atoms_to_score(core.get_string_vector_option("backrub:pivot_atoms"))
# didn't port centroid constraints
#core.scoring.constraints.add_fa_constraints_from_cmdline_to_scorefxn(score_fxn)

# set up the BackrubMover
backrubmover = protocols.moves.BackrubMover()
# read known and unknown optimization parameters from the database
backrubmover.branchopt().read_database()
# TypeError: No to_python (by-value) converter found for C++ type: utility::pointer::owning_ptr<core::mm::MMBondAngleResidueTypeParamSet const>
#if energymethodoptions.bond_angle_residue_type_param_set():
#	backrubmover.branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set())

# set up the SmallMover
smallmover = protocols.moves.SmallMover();
smallmover.nmoves(1);
if options.sm_prob > 0:
	movemap = core.kinematics.MoveMap()
	movemap.init_from_file(core.get_file_option("in:file:movemap"))
	smallmover.movemap(movemap)

# set up the SidechainMover
sidechainmover = protocols.moves.SidechainMover()
sidechainmover.set_task_factory(main_task_factory)
sidechainmover.set_prob_uniform(options.sc_prob_uniform)

# set up the PackRotamersMover
packrotamersmover = protocols.moves.PackRotamersMover()
packrotamersmover.task_factory(main_task_factory)
packrotamersmover.score_function(score_fxn)

input_jobs = protocols.jobdist.load_s_and_l();

for jobnum in xrange(1, len(input_jobs)+1):

	print "Processing " + input_jobs[jobnum].input_tag() + "..."
	custom_fold_tree = False

	input_pose = core.pose.Pose()

	# didn't port centroid input
	rosetta.core.io.pdb.pose_from_pdb(input_pose, input_jobs[jobnum].input_tag())
	#custom_fold_tree = read_fold_tree_from_file(input_pose, input_jobs[jobnum].input_tag())
	#core.scoring.constraints.add_fa_constraints_from_cmdline_to_pose(input_pose)
	#input_pose.dump_pdb(input_jobs[jobnum].output_tag(0) + "_postread.pdb")

	backrubmover.clear_segments()
	backrubmover.set_input_pose(input_pose)

	print "Backrub segment lengths: " #+ option[ backrub::min_atoms ] + "-" + option[ backrub::max_atoms ] + " atoms"

	print "Backrub main chain pivot atoms: " #+ option[ backrub::pivot_atoms ].value_string()

	backrubmover.add_mainchain_segments_from_options()

	print "Backrub Segments Added: " + backrubmover.num_segments()

	print "Score After PDB Load:"
	score_fxn.show(input_pose)

	backrubmover.optimize_branch_angles(input_pose)
	#input_pose.dump_pdb(input_jobs[jobnum].output_tag(0) + "_postoptbranch.pdb")
	sidechainmover.idealize_sidechains(input_pose)
	#input_pose.dump_pdb(input_jobs[jobnum].output_tag(0) + "_postidealizesc.pdb")



