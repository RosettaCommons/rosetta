#! /usr/bin/python
# List of commands used in PyRosetts Workshop #5

from __future__ import print_function

from rosetta import *
from pyrosetta import *
from pyrosetta.teaching import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

start = pose_from_file("../test/data/workshops/1YY8.clean.pdb")
test = Pose()
test.assign(start)

start.pdb_info().name("start")
test.pdb_info().name("test")

pmm = PyMolMover()
pmm.apply(start)
pmm.apply(test)
pmm.keep_history(True)
print( pmm )

# Small and Shear Moves
kT = 1.0
n_moves = 1
movemap = MoveMap()
movemap.set_bb(True)
small_mover = SmallMover(movemap, kT, n_moves)
shear_mover = ShearMover(movemap, kT, n_moves)

small_mover.angle_max("H", 5)
small_mover.angle_max("E", 5)
small_mover.angle_max("L", 5)

small_mover.apply(test)
print( small_mover )
print( shear_mover )

test2 = Pose()
test2.assign(start)
test2.pdb_info().name("test2")
pmm.apply(test2)

movemap.set_bb(False)
movemap.set_bb(50, True)
movemap.set_bb(51, True)
print( movemap )

small_mover.apply(test)
shear_mover.apply(test2)
pmm.apply(test)
pmm.apply(test2)

# Minimization Moves
min_mover = MinMover()

mm4060 = MoveMap()
mm4060.set_bb_true_range(40, 60)

scorefxn = create_score_function("talaris2013")

min_mover.movemap(mm4060)
min_mover.score_function(scorefxn)

# Commenting out for now because this lead to seg-fault in debug builds
# AddPyMolObserver(test2, True)

min_mover.apply(test2)
print( min_mover )

# Monte Carlo Object
mc = MonteCarlo(test, scorefxn, kT)

mc.boltzmann(test)

mc.show_scores()
mc.show_counters()
mc.show_state()

# Trial Mover
trial_mover = TrialMover(small_mover, mc)
for i in range (2): # was 10, swtching to 2 to decrese test running time in debug mode
    trial_mover.apply(test)

print( trial_mover.num_accepts() )
print( trial_mover.acceptance_rate() )
mc.show_state()

# Sequence and Repeat Movers
seq_mover = protocols.moves.SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)
print( seq_mover )

trialmover = TrialMover(seq_mover, mc)
print( trialmover )
repeat_mover = RepeatMover(trialmover, 2) # was 5, swtching to 2 to decrese test running time in debug mode
repeat_mover.apply(test)
mc.show_state()
print( repeat_mover )
# Refinement Protocol
relax = protocols.relax.ClassicRelax()
relax.set_scorefxn(scorefxn)
#relax.apply(test)  # This takes way too long for a test.
