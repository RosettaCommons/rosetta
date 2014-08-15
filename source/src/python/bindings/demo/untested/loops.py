# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
import math

init()

p = Pose()
pose_from_pdb(p, "test_in.pdb")

starting_p = Pose()
starting_p.assign(p)

scorefxn_low = create_score_function('cen_std')
scorefxn_high = create_score_function_ws_patch('standard', 'score12')

loop_begin = 77
loop_end = 85
cutpoint = 81
my_loop = Loop( loop_begin, loop_end, cutpoint)
print my_loop

set_single_loop_fold_tree(p, my_loop)
#set_loop_cutpoint_in_pose_fold_tree(cutpoint, p, loop_begin, loop_end) # rename and use my_loop object as argument

movemap = MoveMap()
movemap.set_bb_true_range(loop_begin, loop_end)
movemap.set_chi( True )

print p.fold_tree()

print "setting up movers"
#backbone movers
fragset3mer = ConstantLengthFragSet(3, "test_in3_fragments")
mover_3mer = ClassicFragmentMover(fragset3mer,movemap)
ccd_closure = CCDLoopClosureMover(my_loop, movemap)

#centroid/fullatom conversion movers
to_centroid = protocols::simple_moves::SwitchResidueTypeSetMover('centroid')
to_fullatom = protocols::simple_moves::SwitchResidueTypeSetMover('fa_standard')
recover_sidechains = protocols::simple_moves::ReturnSidechainMover(starting_p)

#set up sidechain packer movers
task_pack = TaskFactory.create_packer_task(starting_p)
task_pack.restrict_to_repacking()
task_pack.or_include_current( True )
pack = protocols::simple_moves::PackRotamersMover( scorefxn_high, task_pack )

#convert to centroid mode
to_centroid.apply(p)

starting_p_centroid = Pose()
starting_p_centroid.assign(p)

print "set up job distributor"
jd = Job_dist("loop_output", 100, scorefxn_high)
jd.native_pose = starting_p

while (jd.job_complete == False):
  p.assign(starting_p_centroid)

  print "randomizing loop"
  for i in range(loop_begin, loop_end+1):
    p.set_phi(i, -180)
    p.set_psi(i, 180)

  for i in range(loop_begin, loop_end+1):
    mover_3mer.apply(p)

  print "low res loop modeling"
  outer_cycles = 10
  inner_cycles = 30
  init_temp = 2.0
  final_temp = 0.8
  gamma = math.pow((final_temp/init_temp),(1.0/(outer_cycles*inner_cycles)))
  kT = init_temp

  mc = MonteCarlo(p, scorefxn_low, kT)

  for i in range(1,outer_cycles+1):
    mc.recover_low(p)
    scorefxn_low(p)
    for j in range(1, inner_cycles+1):
      kT = kT * gamma
      mc.set_temperature(kT)
      mover_3mer.apply(p)
      ccd_closure.apply(p)
      mc.boltzmann(p)
  mc.recover_low(p)

  print "high-res refinement"

  to_fullatom.apply(p)
  recover_sidechains.apply(p)
  pack.apply(p)

  my_loops = Loops()
  my_loops.add_loop(my_loop)
  loop_refine = LoopMover_Refine_CCD(my_loops)
  loop_refine.apply(p)

  print "outputting decoy..."
  Lrms = loop_rmsd(p, starting_p, my_loops, True)
  jd.additional_decoy_info = " Lrmsd: " + str(Lrms)
  jd.output_decoy(p)

