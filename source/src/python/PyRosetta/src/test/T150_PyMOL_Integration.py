#!/usr/bin/env python
# :noTabs=true:

from __future__ import print_function

import rosetta, pyrosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose = pyrosetta.pose_from_sequence('EVAAAVAT')

pymol = pyrosetta.PyMOLMover()

pymol.apply(pose)

scorefxn = pyrosetta.get_fa_scorefxn() #  rosetta.create_score_function('standard')
scorefxn(pose)

pymol.send_energy(pose)
#pymol.send_energy(pose, label=True)  DEPRECATED: needed to be ported to C++ version

#pymol.send_colors(pose, {}, default_color="orange")
pymol.send_colors(pose, rosetta.std.map_int_int(), default_color=rosetta.protocols.moves.XC_orange)

#colors = {2: "red", 5: "white"}
colors = rosetta.std.map_int_int();
colors[2] = rosetta.protocols.moves.XC_red
colors[5] = rosetta.protocols.moves.XC_white
pymol.send_colors(pose, colors, default_color=rosetta.protocols.moves.XC_blue)

#pymol.label_energy(pose, "fa_atr")  DEPRECATED: needed to be ported to C++ version

#pymol.send_hbonds(pose)  DEPRECATED: needed to be ported to C++ version
#pymol.send_ss(pose)  DEPRECATED: needed to be ported to C++ version
#pymol.send_polars(pose)  DEPRECATED: needed to be ported to C++ version

mm = pyrosetta.MoveMap()
#pymol.send_movemap(pose, mm)  DEPRECATED: needed to be ported to C++ version

#pymol.send_foldtree(pose)  DEPRECATED: needed to be ported to C++ version
#pymol.view_foldtree_diagram(pose)  DEPRECATED: needed to be ported to C++ version

#pymol.plot_graph("Line", "white", [0, 1, 2, 3, 4], [0, 2, 4, 6, 8])  DEPRECATED: needed to be ported to C++ version
#pymol.send_point("Line", "white", 5, 10)  DEPRECATED: needed to be ported to C++ version

observer = rosetta.protocols.moves.AddPyMOLObserver(pose)
pose.set_psi(3, 10)
scorefxn(pose)

rosetta.protocols.moves.AddPyMOLObserver(pose, keep_history=True)
pose.set_psi(2, 10)

pose.set_psi(4, 10)

observer.detach(pose)
