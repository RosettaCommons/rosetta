# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
from rosetta import *
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

args = [ "app",
         "-database minirosetta_database", \
         "-loops:fast" ] # reduce the number of cycles for testing purposes
init( *args )

p = Pose()
pose_from_file( p, "2cpl_min.pdb" )

starting_p = Pose()
starting_p.assign( p )

scorefxn_low = create_score_function( 'cen_std' )
scorefxn_high = create_score_function( 'standard', 'score12' )

pymol = rosetta.PyMOLMover() # If Pymol server is running, centroid stage will display

loop_begin = 145
loop_end = 155
loop_cut = 150
my_loop = Loop( loop_begin, loop_end, loop_cut )
my_loops = Loops()
my_loops.add_loop( my_loop )
print my_loop

set_single_loop_fold_tree( p, my_loop )

movemap = MoveMap()
movemap.set_bb_true_range(loop_begin, loop_end )
movemap.set_chi( True )

print p.fold_tree()

print "setting up movers"
# use the KinematicMover explicitly in centroid stage
kic_mover = KinematicMover()

#centroid/fullatom conversion movers
to_centroid = protocols.simple_moves.SwitchResidueTypeSetMover( 'centroid' )
to_fullatom = protocols.simple_moves.SwitchResidueTypeSetMover( 'fa_standard' )
recover_sidechains = protocols.simple_moves.ReturnSidechainMover( starting_p )

#set up sidechain packer movers
task_pack = TaskFactory.create_packer_task( starting_p )
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

print "building random loop conformation with ideal bond lengths, bond angles,"
print "and omega angles using KIC",
success = False
kic_mover.set_idealize_loop_first( True )
kic_mover.set_pivots( loop_begin, loop_cut, loop_end )
pymol.apply( p )
for i in range( MAX_KIC_BUILD_ATTEMPTS ):
  print "\n  attempt %d..." %i,
  kic_mover.apply( p )
  pymol.apply( p )
  if kic_mover.last_move_succeeded():
    success = True
    kic_mover.set_idealize_loop_first( False )
    print "succeeded."
    break
if not success:
  print "Could not complete initial KIC loop building in %d attempts. Exiting" \
      %MAX_KIC_BUILD_ATTEMPTS
  exit()
scorefxn_low( p )
linmin_mover.apply( p )

print "centroid stage loop KIC remodeling"
outer_cycles = 10
inner_cycles = 30
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
    middle_offset = ( kic_end - kic_start ) / 2
    kic_middle = kic_start + middle_offset
    kic_mover.set_pivots( kic_start, kic_middle, kic_end )
    kic_mover.set_temperature( kT )
    kic_mover.apply( p )
    linmin_mover.apply( p )
    mc.boltzmann( p )
    pymol.apply( p )
    rms = loop_rmsd( p, starting_p, my_loops )
    print "centroid stage loop rmsd to starting loop: %f" %rms
mc.recover_low( p )

print "high-res KIC refinement"
to_fullatom.apply( p )
recover_sidechains.apply( p )
pack.apply( p )

loop_refine = LoopMover_Refine_KIC( my_loops )
loop_refine.apply( p ) # won't show in Pymol for efficiency
pymol.apply( p ) # just show refined model

print "outputting decoy test_kic.pdb"
p.dump_pdb( "test_kic.pdb" )
