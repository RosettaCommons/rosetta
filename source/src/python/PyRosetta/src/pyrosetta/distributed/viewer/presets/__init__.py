# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

"""
Display preset custom viewers for routine visualizations.

Example Jupyter notebook commands:
--------------------------------------------------------------------------------
import pyrosetta.distributed.viewer as viewer
viewer.presets.ligandsAndMetals(poses, window_size=(800, 600), continuous_update=True)
--------------------------------------------------------------------------------

Contribute your own preset Viewer:
 1. Edit ~Rosetta/main/source/src/python/PyRosetta/src/pyrosetta/distributed/viewer/presets/__init__.py 
 2. Copy and modify the "templatePreset" function, renaming it to the name of your new preset Viewer.
 3. Add the name of your new preset Viewer to the __all__ list.
"""

import logging
import pyrosetta

from pyrosetta.distributed import viewer


__all__ = [
    "coreBoundarySurface",
    "ligandsAndMetals"
]

_logger = logging.getLogger("pyrosetta.distributed.viewer.presets")


def coreBoundarySurface(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Display core residues as 'blackCarbon' sticks, boundary residues as 'greyCarbon' sticks, and surface residues
    as 'whiteCarbon' sticks, with 'spectrum' cartoon representation, using the default arguments in
    pyrosetta.rosetta.core.select.residue_selector.LayerSelector() to select layers.

    @klimaj
    """
    core_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    core_selector.set_layers(True, False, False)
    boundary_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    boundary_selector.set_layers(False, True, False)
    surface_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    surface_selector.set_layers(False, False, True)

    view = viewer.init(packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs)
    view.add(viewer.setStyle())
    view.add(
        viewer.setStyle(residue_selector=core_selector, style="stick", colorscheme="blackCarbon", radius=0.25, label=False)
    )
    view.add(
        viewer.setStyle(residue_selector=boundary_selector, style="stick", colorscheme="greyCarbon", radius=0.25, label=False)
    )
    view.add(
        viewer.setStyle(residue_selector=surface_selector, style="stick", colorscheme="whiteCarbon", radius=0.25, label=False)
    )
    view.add(viewer.setDisulfides(radius=0.25))

    return view.show()


def ligandsAndMetals(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Display residues with ResidueProperty.LIGAND as 'brownCarbon' sticks with opaque surface,
    and ResidueProperty.METAL as 'chainHetatm' spheres, with 'spectrum' cartoon representation,
    disulfide bonds, polar hydrogens, and dashed hydrogen bonds.

    @klimaj
    """
    metals_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
         pyrosetta.rosetta.core.chemical.ResidueProperty(31)
    )
    ligands_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
         pyrosetta.rosetta.core.chemical.ResidueProperty(2)
    )

    view = viewer.init(packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs) \
        + viewer.setBackgroundColor("lightgrey") \
        + viewer.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15) \
        + viewer.setStyle(residue_selector=ligands_selector, style="stick", colorscheme="brownCarbon", radius=0.5, label=True) \
        + viewer.setStyle(residue_selector=metals_selector, style="sphere", colorscheme="chainHetatm", radius=1.5, label=True) \
        + viewer.setHydrogenBonds() \
        + viewer.setDisulfides(radius=0.15) \
        + viewer.setHydrogens(color="white", radius=0.033, polar_only=True) \
        + viewer.setSurface(residue_selector=ligands_selector, surface_type="VDW", opacity=0.5, colorscheme="brownCarbon") \
        + viewer.setZoomTo(residue_selector=ligands_selector)

    return view()


def templatePreset(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Add a description of the preset Viewer here

    @author
    """
    view = viewer.init(packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs)
    
    # Add custom Viewer commands here

    return view()
