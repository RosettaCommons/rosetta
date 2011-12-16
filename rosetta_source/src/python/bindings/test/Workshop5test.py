#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from rosetta import *
init()

start = pose_from_pdb("test/data/workshops/1YY8.clean.pdb")
test = Pose()
test.assign(start)

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

test.assign(start)
test2 = Pose()
test2.assign(start)

movemap.set_bb(False)
movemap.set_bb(50, True)
movemap.set_bb(51, True)

small_mover.apply(test)
shear_mover.apply(test2)

# Minimization Moves
min_mover = MinMover()

mm4060 = MoveMap()
mm4060.set_bb_true_range(40, 60)

scorefxn = create_score_function("standard")

min_mover.movemap(mm4060)
min_mover.score_function(scorefxn)

min_mover.apply(test2)

# Monte Carlo Object
mc = MonteCarlo(test, scorefxn, kT)

mc.boltzmann(test)

mc.show_scores()
mc.show_counters()
mc.show_state()

# Trial Mover
trial_mover = TrialMover(small_mover, mc)
trial_mover.apply(test)

print trial_mover.num_accepts()
print trial_mover.acceptance_rate()
mc.show_state()

# Sequence and Repeat Movers
seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)

repeat_mover = RepeatMover(trial_mover, 5)

# Refinement Protocol
relax = ClassicRelax()
relax.set_scorefxn(scorefxn)
#relax.apply(test)  # This takes way too long for a test.
