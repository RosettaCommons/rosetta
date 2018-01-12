# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Dan Mandell

from __future__ import print_function

from pyrosetta import *
from pyrosetta.rosetta import *

from pyrosetta.rosetta.protocols.loops.loop_closure.kinematic_closure import *
import pyrosetta.rosetta.protocols.loops.loop_mover.refine

from sys import exit
from random import randrange
import math

# Demonstrates loop remodeling with kinematic closure (KIC).
# Centroid stage is coded explicitly in PyRosetta, and will
# display in PyMOL if PyMOLPyRosettaServer.py in running in
# Pymol. All-atom refine stage is called through LoopMover_Refine_KIC
# and will not display in PyMOL for efficiency, but will output
# loop rms and energy to command line terminal.

MAX_KIC_BUILD_ATTEMPTS = 10000

# Kale: to reduce test time even further try options: -loops:outer_cycles 1 -loops:max_inner_cycles 1
init(extra_options='-constant_seed -run:test_cycles True')  # -loops:test_cycles is only needed for self-test, please make sure to remove it on production run
# WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

p = core.import_pose.pose_from_file( "../test/data/2cpl_min.pdb" )

starting_p = Pose()
starting_p.assign( p )

scorefxn_low  = protocols.loops.get_cen_scorefxn() #  create_score_function( 'cen_std' )
scorefxn_high = protocols.loops.get_fa_scorefxn()  #  create_score_function_ws_patch( 'standard', 'score12' )

pymol = PyMOLMover() # If Pymol server is running, centroid stage will display

loop_begin = 145
loop_end = 155
loop_cut = 150
my_loop = protocols.loops.Loop( loop_begin, loop_end, loop_cut )
my_loops = protocols.loops.Loops()
my_loops.add_loop( my_loop )
print(my_loop)

protocols.loops.set_single_loop_fold_tree( p, my_loop )

movemap = MoveMap()
movemap.set_bb_true_range(loop_begin, loop_end )
movemap.set_chi( True )

print( p.fold_tree() )

print( "setting up movers" )
# use the KinematicMover explicitly in centroid stage
kic_mover = KinematicMover()

#centroid/fullatom conversion movers
to_centroid = protocols.simple_moves.SwitchResidueTypeSetMover( 'centroid' )
to_fullatom = protocols.simple_moves.SwitchResidueTypeSetMover( 'fa_standard' )
recover_sidechains = protocols.simple_moves.ReturnSidechainMover( starting_p )

#set up sidechain packer movers
task_pack = core.pack.task.TaskFactory.create_packer_task( starting_p )
task_pack.restrict_to_repacking()
task_pack.or_include_current( True )
pack = protocols.minimization_packing.PackRotamersMover( scorefxn_high, task_pack )

#convert to centroid mode
to_centroid.apply( p )

#set up centroid stage MoveMap
mm = MoveMap()
movemap.set_bb( True )

# set up centroid stage line minimizer
tol = 0.001
min_type = "linmin"
linmin_mover = protocols.minimization_packing.MinMover( mm, scorefxn_low, min_type, tol, True )

# save starting pose
starting_p_centroid = Pose()
starting_p_centroid.assign( p )

print( "building random loop conformation with ideal bond lengths, bond angles," )
print( "and omega angles using KIC", end='' )
success = False
kic_mover.set_idealize_loop_first( True )
kic_mover.set_pivots( loop_begin, loop_cut, loop_end )
pymol.apply( p )
for i in range( MAX_KIC_BUILD_ATTEMPTS ):
  print( "\n  attempt %d..." %i, end='' )
  kic_mover.apply( p )
  pymol.apply( p )
  if kic_mover.last_move_succeeded():
    success = True
    kic_mover.set_idealize_loop_first( False )
    print( "succeeded." )
    break
if not success:
  print( "Could not complete initial KIC loop building in %d attempts. Exiting" \
      %MAX_KIC_BUILD_ATTEMPTS )
  exit()
scorefxn_low( p )
linmin_mover.apply( p )

print( "centroid stage loop KIC remodeling" )
outer_cycles = 1 # 10 # inner_cycles and outer_cycles should be much higher in production run!
inner_cycles = 1 # 30
init_temp = 2.0
final_temp = 1.0
gamma = math.pow( ( final_temp/init_temp ),( 1.0/( outer_cycles*inner_cycles ) ) )
kT = init_temp
mc = MonteCarlo( p, scorefxn_low, kT )

for i in range( 1,outer_cycles+1 ):
  mc.recover_low( p )
  scorefxn_low( p )
  for j in range( 1, inner_cycles+1 ):
    kT = kT * gamma
    mc.set_temperature( kT )
    kic_start = randrange( loop_begin, loop_end - 1 )
    kic_end = randrange( kic_start+2, loop_end+1 )
    middle_offset = ( kic_end - kic_start ) // 2
    kic_middle = kic_start + middle_offset
    kic_mover.set_pivots( kic_start, kic_middle, kic_end )
    kic_mover.set_temperature( kT )
    kic_mover.apply( p )
    linmin_mover.apply( p )
    mc.boltzmann( p )
    pymol.apply( p )
    rms = protocols.loops.loop_rmsd( p, starting_p, my_loops )
    print( "centroid stage loop rmsd to starting loop: %f" % rms )
mc.recover_low( p )

print( "high-res KIC refinement" )
to_fullatom.apply( p )
recover_sidechains.apply( p )
pack.apply( p )

loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_KIC( my_loops )
loop_refine.apply( p ) # won't show in Pymol for efficiency
pymol.apply( p ) # just show refined model

print( "outputting decoy .test_kic.pdb" )
p.dump_pdb( ".test_kic.pdb" )
