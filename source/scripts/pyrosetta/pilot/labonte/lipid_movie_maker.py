from rosetta import *
from random import *

init(extra_options="-include_sugars -include_lipids -read_pdb_link_records " +
                   "-enable_branching -mute core -mute protocols")

FILEPATH = "lipid_movie/GalCer_frame"

pose = pose_from_pdb("../Rosetta3/Rosetta/main/tests/integration/tests/carbohydrates/input/GalCer.pdb")
pose.pdb_info().name("GalCer")

sf = get_fa_scorefxn()

initial_score = sf(pose)

mm = MoveMap()
mm.set_bb(True)
mm.set_chi(True)

pt = standard_packer_task(pose)
pt.restrict_to_repacking()
pt.temporarily_set_pack_residue(2, False)
pt.or_include_current(True)

min_move = MinMover(mm, sf, "linmin", 0.001, True)

pack = PackRotamersMover(sf, pt)

mc = MonteCarlo(pose, sf, 0.5)

pm = PyMOL_Mover()

frame = 0
pm.apply(pose)
pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
frame += 1

for i in range(20):
    resnum = 2
    while resnum == 2:
        resnum = randint(1, 3)
        if resnum != 2:
            res = pose.residue(resnum)
            torsion_num = randint(1, res.n_mainchain_atoms())
            id = TorsionID(resnum, BB, torsion_num)
            current_angle = pose.torsion(id)
            pose.set_torsion(id, current_angle + gauss(0, 30))
    pm.apply(pose)
    pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
    frame += 1

    resnum = 2
    while resnum == 2:
        resnum = randint(1, 3)
        if resnum != 2:
            res = pose.residue(resnum)
            torsion_num = randint(1, res.n_mainchain_atoms())
            id = TorsionID(resnum, BB, torsion_num)
            current_angle = pose.torsion(id)
            pose.set_torsion(id, current_angle + gauss(0, 30))
    pm.apply(pose)
    pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
    frame += 1

    pack.apply(pose)
    pm.apply(pose)
    pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
    frame += 1

    min_move.apply(pose)
    pm.apply(pose)
    pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
    frame += 1

    mc.boltzmann(pose)
    pm.apply(pose)
    pose.dump_pdb(FILEPATH + str(frame) + ".pdb")
    frame += 1

print 'Initial Score:', initial_score
print 'Final Score:', sf(pose)