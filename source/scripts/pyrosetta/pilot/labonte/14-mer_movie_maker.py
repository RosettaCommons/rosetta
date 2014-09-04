from rosetta import *
from rosetta.core.chemical.carbohydrates import CarbohydrateInfo
from rosetta.protocols.simple_moves import RingConformationMover

init(extra_options="-include_sugars -read_pdb_link_records " +
                   "-enable_branching -mute core -mute protocols")

pose = pose_from_pdb("N-linked_14-mer_glycan.pdb")

n_res = pose.total_residue()
first_sugar = 6
last_sugar = n_res

sf = get_fa_scorefxn()

initial_score = sf(pose)

mm_flips = MoveMap()
mm_flips.set_bb_true_range(first_sugar, last_sugar)

mm_small_moves = MoveMap()
mm_small_moves.set_bb_true_range(first_sugar, last_sugar)
mm_small_moves.set_bb(6, False)
#mm_small_moves.set_bb(15, False)
#mm_small_moves.set_bb(18, False)

pt = standard_packer_task(pose)
pt.restrict_to_repacking()
pt.temporarily_set_pack_residue(2, False)
pt.or_include_current(True)

ring_flip = RingConformationMover()
ring_flip.movemap(mm_flips)

sm_move = SmallMover()
sm_move.movemap(mm_small_moves)
sm_move.angle_max(30.0)

pack = PackRotamersMover(sf, pt)

mc = MonteCarlo(pose, sf, 0.5)

observer = PyMOLObserver(pose, True)

frame = 0
pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
frame += 1

for i in range(20):
    ring_flip.apply(pose)
    pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
    frame += 1

    sm_move.apply(pose)
    pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
    frame += 1

    sm_move.apply(pose)
    pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
    frame += 1

    pack.apply(pose)
    pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
    frame += 1

    mc.boltzmann(pose)
    pose.dump_pdb("movie/14-mer_frame" + str(frame) + ".pdb")
    frame += 1

print 'Initial Score:', initial_score
print 'Final Score:', sf(pose)