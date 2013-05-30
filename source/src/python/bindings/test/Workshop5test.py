#! /usr/bin/python
# List of commands used in PyRosetts Workshop #5

from rosetta import *
init()

start = pose_from_pdb("test/data/workshops/1YY8.clean.pdb")
test = Pose()
test.assign(start)

start.pdb_info().name("start")
test.pdb_info().name("test")

pmm = PyMOL_Mover()
pmm.apply(start)
pmm.apply(test)
pmm.keep_history(True)
print pmm

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
print small_mover
print shear_mover

test2 = Pose()
test2.assign(start)
test2.pdb_info().name("test2")
pmm.apply(test2)

movemap.set_bb(False)
movemap.set_bb(50, True)
movemap.set_bb(51, True)
print movemap

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
AddPyMolObserver(test2, True)
min_mover.apply(test2)
print min_mover

# Monte Carlo Object
mc = MonteCarlo(test, scorefxn, kT)

mc.boltzmann(test)

mc.show_scores()
mc.show_counters()
mc.show_state()

# Trial Mover
trial_mover = TrialMover(small_mover, mc)
for i in range (10):
    trial_mover.apply(test)

print trial_mover.num_accepts()
print trial_mover.acceptance_rate()
mc.show_state()

# Sequence and Repeat Movers
seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)
print seq_mover

trialmover = TrialMover(seq_mover, mc)
print trialmover
repeat_mover = RepeatMover(trialmover, 5)
repeat_mover.apply(test)
mc.show_state()
print repeat_mover
# Refinement Protocol
relax = ClassicRelax()
relax.set_scorefxn(scorefxn)
#relax.apply(test)  # This takes way too long for a test.
