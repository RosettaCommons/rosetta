#!/usr/bin/env python
# :noTabs=true:

import rosetta

rosetta.init()

pose = rosetta.pose_from_sequence('EVAAAVAT')

pymol = rosetta.PyMOL_Mover()

pymol.apply(pose)

scorefxn = rosetta.create_score_function('standard')
scorefxn(pose)

pymol.send_energy(pose)

pymol.send_colors(pose, {}, default_color=rosetta.protocols.moves.XC_orange)

colors = { 2 : rosetta.protocols.moves.XC_red, 5 : rosetta.protocols.moves.XC_white}

pymol.send_colors(pose, colors, default_color=rosetta.protocols.moves.XC_blue)


observer = rosetta.AddPyMolObserver(pose)
pose.set_psi(3, 10)
scorefxn(pose)

rosetta.AddPyMolObserver(pose, keep_history=True)
pose.set_psi(2, 10)

pose.set_psi(4, 10)

observer.detach(pose)

