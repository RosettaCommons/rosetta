#!/usr/bin/env python
# :noTabs=true:

import time, bz2

from random import random

import rosetta, rosetta.core.pose

import rosetta.core.pose.signals


def test(pose):
    pymol.update_energy = True
    for t in range(12000):
        pose.set_phi(40, pose.phi(40) + 1)
        pose_s.set_phi(4, pose_s.phi(4) + 1)

        scorefxn(pose)
        scorefxn(pose_s)

        seq.apply(pose)
        seq.apply(pose_s)

        #pymol.send_energy( pose )
        #pymol.send_energy( pose_s )

        time.sleep(.1)


        '''
        for i in range(1, pose.n_residue()+1):
            pose.set_phi(i, random()*100. )
            pose.set_psi(i, random()*100. )
            seq.apply(pose)
            #pose.dump_pdb('_%s.pdb' % i)
            time.sleep(1)
            '''


def coloring_demo(pose):
    pymol.update_energy = False
    pymol.apply(pose)
    for t in range(12000):
        N = pose.total_residue()
        for r in range(1, N+1 ):
            C = { r : rosetta.protocols.moves.XC_red, 1+N-r: rosetta.protocols.moves.XC_white}
            #print r, N, C
            pymol.send_colors(pose, C, default_color=rosetta.protocols.moves.XC_blue)

            #pymol.send_energy( pose_s )

            time.sleep(.1)


rosetta.init()

pose = rosetta.Pose();  pose.name = 'CustomNamedPose'
pose_s = rosetta.Pose()
rosetta.pose_from_pdb(pose, "../test/data/test_in.pdb")
rosetta.pose_from_pdb(pose_s, "../test/data/test_in_short.pdb")

scorefxn = rosetta.create_score_function('standard')
scorefxn(pose)

pymol = rosetta.PyMOL_Mover()

pymol.apply(pose_s)
coloring_demo(pose_s)


seq = rosetta.protocols.moves.SequenceMover()
seq.add_mover(pymol)

seq.apply(pose)
seq.apply(pose_s)

#pm.sendEnergies(pose, 'some_thing', [1,2,3])
#pm.sendEnergies(pose, 'some_thing_other', [1,2,3./7., 1/3.])


po = rosetta.PyMOL_Observer(keep_history=True)
po.add_observer(pose)


#from rosetta import *


scorefxn(pose_s)

#pymol.send_specific_energy( pose , 'total_energy')
#pymol.send_specific_energy( pose_s , 'total_energy')

