# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.


import logging
import py3Dmol
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io

from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.pose.full_model_info import (
    get_res_num_from_pdb_info,
    get_chains_from_pdb_info
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector
)


_logger = logging.getLogger("pyrosetta.distributed.viewer")


class setBackgroundColor:
    """
    Set viewer background color with either Hexcode or standard colors.
    Default: 0xffffffff
    """
    def __init__(self, color=0xffffffff):

        self.color = color

    def apply(self, viewer, pose, pdbstring):

        viewer.setBackgroundColor(self.color)

        return viewer


class setDisulfides:
    """
    Display disulfide bonds according to pyrosetta.rosetta.core.conformation.is_disulfide_bond
    for all combinations of cysteine residues in the Pose.
    """
    def __init__(self, color="gold", radius=0.5):

        self.color = color
        self.radius = radius

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        if not pose:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        cys_res = []
        for i, aa in enumerate(pose.sequence(), start=1):
            if aa == "C":
                cys_res.append(i)
        for i in cys_res:
            for j in cys_res:
                if pyrosetta.rosetta.core.conformation.is_disulfide_bond(pose.conformation(), i, j):
                    i_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(pose.residue(i).atom_index("SG"), i)
                    )
                    j_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(pose.residue(j).atom_index("SG"), j)
                    )
                    viewer.addCylinder(
                        {
                            "radius": self.radius, "color": self.color, "fromCap": True, "toCap": True,
                            "start": {"x": i_xyz[0], "y": i_xyz[1], "z": i_xyz[2]},
                            "end": {"x": j_xyz[0], "y": j_xyz[1], "z": j_xyz[2]}
                        }
                    )
        
        return viewer


class setHydrogenBonds:
    """Show hydrogen bonds according to pose.get_hbonds()."""
    def __init__(self, color="black", dashed=True, radius=None):

        self.dashed = dashed
        self.color = color
        self.radius = radius

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        if not pose:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        hbond_set = pose.get_hbonds()
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_xyz = pose.residue(r.don_res()).xyz(r.don_hatm())
                    acc_xyz = pose.residue(r.acc_res()).xyz(r.acc_atm())
                    if self.radius:
                        if self.dashed:
                            _logger.warning(
                                " ".join("setHydrogenBonds argument 'radius' cannot be set with argument 'dashed' set to True. \
                                Setting argument 'dashed' to False.".split())
                            )
                        viewer.addCylinder(
                            {
                                "radius": self.radius, "color": self.color, "fromCap": True, "toCap": True,
                                "start": {"x": don_xyz[0], "y": don_xyz[1], "z": don_xyz[2]},
                                "end": {"x": acc_xyz[0], "y": acc_xyz[1], "z": acc_xyz[2]}
                            }
                        )
                    else:
                        viewer.addLine(
                            {
                                "dashed": self.dashed, "color": self.color,
                                "start": {"x": don_xyz[0], "y": don_xyz[1], "z": don_xyz[2]},
                                "end": {"x": acc_xyz[0], "y": acc_xyz[1], "z": acc_xyz[2]}
                            }
                        )
        
        return viewer


class setHydrogens:
    """Show all or only polar hydrogen atoms."""
    def __init__(self, color="white", radius=0.05, polar_only=False):

        self.color = color
        self.radius = radius
        self.polar_only = polar_only

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        def _addCylinder(_viewer, i_xyz, j_xyz):
            _viewer.addCylinder(
                {
                    "radius": self.radius, "color": self.color, "fromCap": True, "toCap": True,
                    "start": {"x": i_xyz[0], "y": i_xyz[1], "z": i_xyz[2]},
                    "end": {"x": j_xyz[0], "y": j_xyz[1], "z": j_xyz[2]}
                }
            )
            return _viewer

        if not pose:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        if pose.is_fullatom():
            for i in range(1, pose.total_residue() + 1):
                r = pose.residue(i)
                h_begin = r.attached_H_begin()
                h_end = r.attached_H_end()
                for h in range(1, len(h_begin) + 1):
                    i_index = h_begin[h]
                    j_index = h_end[h]
                    if all(q != 0 for q in [i_index, j_index]):
                        i_xyz = r.atom(h).xyz()
                        for j in range(i_index, j_index + 1):
                            if self.polar_only:
                                if r.atom_is_polar_hydrogen(j):
                                    j_xyz = r.atom(j).xyz()
                                    viewer = _addCylinder(viewer, i_xyz, j_xyz)
                            else:
                                j_xyz = r.atom(j).xyz()
                                viewer = _addCylinder(viewer, i_xyz, j_xyz)

        return viewer


class setStyle:
    """
    Show and color cartoon, and/or show heavy atoms with style 'line', 'cross', 'stick', or 'sphere'
    with provided 'color' and 'radius'. If 'residue_selector' argument is provided, apply styles only
    to selected residues. If the 'command' argument is provided, override all other arguments and 
    pass py3Dmol.view.setStyle commands directly to the Viewer object.

    Example:
    
    view = viewer.init(poses) \
        + viewer.setStyle(
            command=
                {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}}
        )
    """
    def __init__(self, residue_selector=None,
                 cartoon=True, cartoon_color="spectrum",
                 style="stick", colorscheme="blackCarbon", radius="0.1",
                 label=True, label_fontsize=12, label_background=False, label_fontcolor="black",
                 command=None):

        _valid_styles = ["line", "cross", "stick", "sphere"]
        if not any(style == s for s in _valid_styles): 
            raise "setStyle argument 'style' must be either: {0}".format(", ".join(_valid_styles))

        if residue_selector:
            if not isinstance(residue_selector, ResidueSelector):
                raise ViewerInputError(residue_selector)

        self.residue_selector = residue_selector
        self.cartoon = cartoon
        self.cartoon_color = cartoon_color
        self.style = style
        self.colorscheme = colorscheme
        if float(radius) == 0.:
            radius = 1e-10
        self.radius = radius
        self.label = label
        self.label_fontsize = label_fontsize
        self.label_background = label_background
        self.label_fontcolor = label_fontcolor
        self.command = command

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        if self.command:
            if isinstance(cmd, tuple):
                viewer.setStyle(*self.command)
            elif isinstance(cmd, dict):
                viewer.setStyle(self.command)
            else:
                raise ValueError("setStyle argument 'command' should be an instance of tuple or dict.")
        else: 
            if self.residue_selector:
                if not pose:
                    pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)
                
                resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)
                
                if (not resi) and (not chain):
                    pass
                else:
                    if self.cartoon:
                        viewer.setStyle(
                            {"resi": resi, "chain": chain}, {"cartoon": {"color": self.cartoon_color},
                            self.style: {"colorscheme": self.colorscheme, "radius": self.radius}}
                        )
                    else:
                        viewer.setStyle(
                            {"resi": resi, "chain": chain},
                            {self.style: {"colorscheme": self.colorscheme, "radius": self.radius}}
                        )
                    if self.label:
                        viewer.addResLabels(
                            {"resi": resi, "chain": chain},
                                {"fontSize": self.label_fontsize,
                                 "showBackground": self.label_background,
                                 "fontColor": self.label_fontcolor}
                        )
            else:
                if self.cartoon:
                    viewer.setStyle(
                        {
                            "cartoon": {"color": self.cartoon_color},
                            self.style: {"colorscheme": self.colorscheme, "radius": self.radius}
                        }
                    )
                else:
                    viewer.setStyle(
                        {
                            self.style: {"colorscheme": self.colorscheme, "radius": self.radius}
                        }
                    )

        return viewer


class setSurface:
    """
    Show surface with argument 'opacity' for residues in 'residue_selector' within the pose.
    py3Dmol supports the following 'surface_type' arguments:

    "VDW": Van der Waals surface
    "MS": Molecular surface
    "SES": Solvent excluded surface
    "SAS": Solvent accessible surface
    """
    def __init__(self, residue_selector=None, surface_type="VDW", opacity=0.5, color=None, colorscheme=None):

        if not residue_selector:
            residue_selector = TrueResidueSelector()
        elif not isinstance(residue_selector, ResidueSelector):
            raise ViewerInputError(residue_selector)
        if not any(surface_type == s for s in ["VDW", "MS", "SES", "SAS"]):
            raise ValueError("Input surface_type argument must be one of the strings: 'VDW', 'MS', 'SES', 'SAS'"
            )
        
        _surface_types_dict = {
            "VDW": py3Dmol.VDW,
            "MS": py3Dmol.MS,
            "SES": py3Dmol.SES,
            "SAS": py3Dmol.SAS,
        }

        self.residue_selector = residue_selector
        self.surface_type = _surface_types_dict[surface_type]
        self.opacity = opacity
        self.color = color
        self.colorscheme = colorscheme

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        if not pose:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)
            
        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            if self.colorscheme:
                viewer.addSurface(
                    self.surface_type, 
                    {"opacity": self.opacity, "colorscheme": self.colorscheme},
                    {"resi": resi, "chain": chain}
                )
            elif self.color:
                viewer.addSurface(
                    self.surface_type, 
                    {"opacity": self.opacity, "color": self.color},
                    {"resi": resi, "chain": chain}
                )
            else:
                viewer.addSurface(
                    self.surface_type, 
                    {"opacity": self.opacity},
                    {"resi": resi, "chain": chain}
                )

        return viewer


class setZoom:
    """
    Set zoom magnification factor. Values >1 zoom in, <1 zoom out.
    Default: 2
    """
    def __init__(self, factor=2):

        self.factor = factor
    
    def apply(self, viewer, pose, pdbstring):

        viewer.zoom(self.factor)

        return viewer


class setZoomTo:
    """Zoom into a ResidueSelector."""
    def __init__(self, residue_selector=None):

        if not residue_selector:
            residue_selector = TrueResidueSelector()
        elif not isinstance(residue_selector, ResidueSelector):
            raise ViewerInputError(residue_selector)
        
        self.residue_selector = residue_selector

    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring):

        if not pose:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)
        
        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            viewer.zoomTo({"resi": resi, "chain": chain})

        return viewer


class ViewerInputError(Exception):
    """Exception raised for errors with the input argument 'residue_selector'."""
    def __init__(self, obj):

        super().__init__(
            " ".join(
                "Input 'residue_selector' argument should be an instance of \
                pyrosetta.rosetta.core.select.residue_selector.ResidueSelector. \
                Input argument 'residue_selector' was invoked with: {0}".format(obj).split()
            )
        )


def _pose_to_residue_chain_tuples(pose, residue_selector, logger=_logger):
    """
    Given a Pose object and ResidueSelector object, return a tuple of lists containing
    PDB residue numbers and chain IDs for the selection.
    """

    pdb_numbering = list(zip(get_res_num_from_pdb_info(pose), get_chains_from_pdb_info(pose)))
    residues_from_subset = list(get_residues_from_subset(residue_selector.apply(pose)))
    residue_chain_tuples = [pdb_numbering[i - 1] for i in residues_from_subset]
    
    if len(residue_chain_tuples) == 0:
        logger.info("ResidueSelector {0} is empty and did not select any residues!".format(residue_selector))
        return [], []
    else:
        return map(list, zip(*residue_chain_tuples))

def _pdbstring_to_pose(pdbstring, class_name, logger=_logger):
    """Convert pdbstring to pose with logging."""
    logger.info(
        " ".join("{0} requires pyrosetta.rosetta.core.pose.Pose object but given input .pdb file. \
        Now instantiating pyrosetta.rosetta.core.pose.Pose object from input .pdb file. \
        For faster performance, either input pyrosetta.rosetta.core.pose.Pose \
        or pyrosetta.distributed.packed_pose.core.PackedPose objects to pyrosetta.distributed.viewer.init, \
        or do not add {0} objects that require a pyrosetta.rosetta.core.pose.Pose object.  \
        ".format(class_name).split())
    )
    return io.to_pose(io.pose_from_pdbstring(pdbstring))
