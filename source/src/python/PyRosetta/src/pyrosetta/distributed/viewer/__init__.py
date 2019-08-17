# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

"""
Display PackedPose or Pose objects, or .pdb files, in py3Dmol within a Jupyter notebook.

Usage:
import pyrosetta.distributed.viewer as viewer

Example Jupyter notebook commands:
--------------------------------------------------------------------------------

view = viewer.init("path/to/pdbfile.pdb")
view.show()

--------------------------------------------------------------------------------

import logging
logging.basicConfig(level=logging.WARNING)
import pyrosetta
pyrosetta.init("-mute all")

pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("5BVL")

view = viewer.init(pose, window_size=(800, 600))
view() # Equivalent to view.show()

--------------------------------------------------------------------------------

poses = [pyrosetta.toolbox.rcsb.pose_from_rcsb(id) for id in ["5BVL", "6MSR", "1QCQ"]]

view = viewer.init(poses) \
+ viewer.setStyle(colorscheme="lightgreyCarbon") \
+ viewer.setHydrogenBonds()
view()

--------------------------------------------------------------------------------

import pyrosetta.distributed.io as io
packed_pose = io.to_packed(pyrosetta.toolbox.pose_from_rcsb("2FD7"))
polar_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
    pyrosetta.rosetta.core.chemical.ResidueProperty(52)
)

view = viewer.init(packed_pose)
view.add(viewer.setStyle(radius=0.1))
view.add(viewer.setStyle(residue_selector=polar_residue_selector, colorscheme="whiteCarbon", radius=0.25, label=False))
view.add(viewer.setHydrogens(color="white", polar_only=True, radius=0.1))
view.add(viewer.setHydrogenBonds(color="black"))
view.add(viewer.setDisulfides(radius=0.1))
view()

--------------------------------------------------------------------------------

view = sum(
    [
        viewer.init(packed_pose),
        viewer.setStyle(cartoon=False, style="sphere", radius=1.5, colorscheme="darkgreyCarbon"),
        viewer.setZoom(factor=1.5)
    ]
)
view.show()

--------------------------------------------------------------------------------

pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("6MSR")
chA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
chB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")

view = sum(
    [
        viewer.init(pose),
        viewer.setStyle(cartoon_color="lightgrey", radius=0.25),
        viewer.setSurface(residue_selector=chA, colorscheme="greenCarbon", opacity=0.65, surface_type="VDW"),
        viewer.setSurface(residue_selector=chB, color="blue", opacity=1.0, surface_type="SAS"),
        viewer.setDisulfides(radius=0.25),
        viewer.setZoom(factor=1.5)
    ]
)
view()

--------------------------------------------------------------------------------

helix_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("H")
sheet_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("E")
loop_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("L")

modules = [
    viewer.setBackgroundColor(color="grey"),
    viewer.setStyle(residue_selector=helix_selector, cartoon_color="blue", label=False, radius=0),
    viewer.setStyle(residue_selector=sheet_selector, cartoon_color="red", label=False, radius=0),
    viewer.setStyle(residue_selector=loop_selector, cartoon_color="white", label=False, radius=0)
]

view = viewer.init(poses, window_size=(1200, 600), modules=modules)
view()

--------------------------------------------------------------------------------

view.reinit() # Subtract all visualization modules previously added to the Viewer
view()

--------------------------------------------------------------------------------

# View live trajectory:

pose = pyrosetta.toolbox.pose_from_rcsb("2FD7")
view = viewer.init(pose, delay=0.15) + viewer.setStyle(radius=0.1) + viewer.setDisulfides(radius=0.1)
backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()

for _ in range(100):
    backrub.apply(pose)
    minimize.apply(pose)
    view.show()

--------------------------------------------------------------------------------

# Display preset custom viewers for routine visualizations:

viewer.presets.coreBoundarySurface(poses, window_size=(800, 600), continuous_update=True)

--------------------------------------------------------------------------------

The Viewer quickly renders .pdb files, dynamically instantiating Pose objects if required
for certain visualization modules (matching the name "viewer.set*"). So when adding
visualization modules to the Viewer or using presets, passing Pose or PackedPose objects to the
Viewer is suggested for quicker rendering. If a Pose object or list, tuple, or set of Pose
objects are provided to the Viewer, the Pose(s) pointer location(s) in memory remain fixed, and so
the Viewer can dynamically update upon Pose conformational changes by calling the show() method.
The Viewer applies visualization modules in the same order they are added (from left to right),
so layering different styles (and ResidueSelectors) on top of one another becomes possible.
The user must have already initialized PyRosetta providing .params files for any ligands and 
non-canonical residues in the Pose, otherwise pyrosetta.distributed automatically initializes
PyRosetta with default options.
"""

import warnings

from pyrosetta.distributed.viewer.core import (
    expand_notebook,
    init
)
from pyrosetta.distributed.viewer.modules import *
from pyrosetta.distributed.viewer import presets


__all__ = [
    "expand_notebook",
    "init",
    "presets",
    "setBackgroundColor",
    "setDisulfides",
    "setHydrogenBonds",
    "setHydrogens",
    "setStyle",
    "setSurface",
    "setZoom",
    "setZoomTo"
]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        get_ipython().Completer.limit_to__all__ = True
    except:
        pass

expand_notebook()
