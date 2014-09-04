from rosetta import *
from rosetta.core.chemical.carbohydrates import CarbohydrateInfo
from rosetta.protocols.simple_moves.carbohydrates import RingConformationMover

init()

pose = pose_from_pdb("maltotriose.pdb")

n_res = pose.total_residue()
first_sugar = n_res - 2
last_sugar = n_res

residues = [pose.residue(resnum) for resnum in range(first_sugar, last_sugar + 1)]

for res in residues:
    print res
    print res.carbohydrate_info()

sf = create_score_function("standard")

sf.show(pose)

mm1 = MoveMap()
mm1.set_bb_true_range(first_sugar + 1, last_sugar)

mm2 = MoveMap()
mm2.set_bb_true_range(first_sugar, last_sugar)

sm_move = SmallMover()
sm_move.movemap(mm1)
sm_move.angle_max(15.0)

ring_flip = RingConformationMover()
ring_flip.movemap(mm2)

mc = MonteCarlo(pose, sf, 0.5)

frame = 0
pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
frame += 1

for i in range(10):
    ring_flip.apply(pose)
    pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
    frame += 1

    sm_move.apply(pose)
    pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
    frame += 1

    sm_move.apply(pose)
    pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
    frame += 1

    sm_move.apply(pose)
    pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
    frame += 1

    mc.boltzmann(pose)
    sm_move.apply(pose)
    pose.dump_pdb("movie/frame" + str(frame) + ".pdb")
    frame += 1

sf.show(pose)