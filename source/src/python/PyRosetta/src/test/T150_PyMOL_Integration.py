#!/usr/bin/env python
# :noTabs=true:

import rosetta

rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose = rosetta.pose_from_sequence('EVAAAVAT')

pymol = rosetta.PyMOL_Mover()

pymol.apply(pose)

scorefxn = rosetta.get_fa_scorefxn() #  rosetta.create_score_function('standard')
scorefxn(pose)

pymol.send_energy(pose)
pymol.send_energy(pose, label=True)

pymol.send_colors(pose, {}, default_color="orange")
colors = {2: "red", 5: "white"}
pymol.send_colors(pose, colors, default_color="blue")

pymol.label_energy(pose, "fa_atr")

pymol.send_hbonds(pose)
pymol.send_ss(pose)
pymol.send_polars(pose)

mm = rosetta.MoveMap()
pymol.send_movemap(pose, mm)

pymol.send_foldtree(pose)
pymol.view_foldtree_diagram(pose)

pymol.plot_graph("Line", "white", [0, 1, 2, 3, 4], [0, 2, 4, 6, 8])
pymol.send_point("Line", "white", 5, 10)

observer = rosetta.AddPyMolObserver(pose)
pose.set_psi(3, 10)
scorefxn(pose)

rosetta.AddPyMolObserver(pose, keep_history=True)
pose.set_psi(2, 10)

pose.set_psi(4, 10)

observer.detach(pose)
