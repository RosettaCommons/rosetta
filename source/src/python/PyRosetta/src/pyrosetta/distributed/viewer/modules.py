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
    Set Viewer background color with either Hexcode or standard colors.
    
    Parameters
    ----------
    first : optional
        `color`

        Hexcode literal (e.g. 0xffffffff) or `str` indicating a standard color (e.g. "black").
        Default: 0xffffffff

    Returns
    -------
    A Viewer instance.
    """

    def __init__(self, color=0xffffffff):

        self.color = color

    def apply(self, viewer, pose, pdbstring):

        viewer.setBackgroundColor(self.color)

        return viewer


class setDisulfides:
    """
    Display disulfide bonds according to `pyrosetta.rosetta.core.conformation.is_disulfide_bond()`
    for all combinations of cysteine residues in each initialized `.pdb` file, `Pose` or `PackedPose` object.
    
    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "black").
        Default: "gold"

    second : optional
        `radius`

        `float` or `int` indicating the radius of the stick connecting the atoms participating in each disulfide bond.
        Default: 0.5

    Returns
    -------
    A Viewer instance.
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
    """
    Display hydrogen bonds according to `pyrosetta.rosetta.core.pose.Pose.get_hbonds()`
    in each initialized `.pdb` file, `Pose` or `PackedPose` object.
    
    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "yellow").
        Default: "black"

    second : optional
        `dashed`

        `True` or `False` to show hydrogen bonds as dashed lines.
        If `True`, then option `radius` must be `None`.
        If `False`, then option `radius` must be specified.
        Default: True

    third : optional
        `radius`

        `float` or `int` indicating the radius of the solid (non-dashed) stick connecting the atoms participating
        in each hydrogen bond. If set, this automatically sets the option `dashed` to `False`.
        Default: None

    Returns
    -------
    A Viewer instance.
    """
    def __init__(self, color="black", dashed=True, radius=None):

        self.color = color
        self.dashed = dashed
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
    """
    Show all or only polar hydrogen atoms in each initialized `.pdb` file, `Pose` or `PackedPose` object.
    
    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "grey").
        Default: "white"

    second : optional
        `radius`

        `float` or `int` indicating the radius of the hydrogen atom stick represnetations.
        Default: 0.05

    third : optional
        `polar_only`

        `True` or `False`. `True` to show only polar hydrogen atoms, and `False` to show all hydrogen atoms.
        Default: False

    Returns
    -------
    A Viewer instance.
    """
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
    Show and color cartoon, and/or show heavy atoms with provided style, color and radius for each initialized
    `.pdb` file, `Pose` or `PackedPose` object. If the `residue_selector` argument is provided, apply styles
    only to the selected residues. If the `command` argument is provided, override all other arguments and pass
    `py3Dmol.view.setStyle()` commands to the Viewer.
    
    Parameters
    ----------
    first : optional
        `residue_selector`
        
        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` on which to apply the style(s).
        Default: None

    second : optional
        `cartoon`
        
        `True` or `False` to show cartoon representation.
        Default: True

    third : optional
        `cartoon_color`
        
        Hexcode literal (e.g. 0xAF10AB) or `str` indicating a standard color (e.g. "grey") for the cartoon representation.
        If "spectrum", apply reversed color gradient based on residue numbers. The option `cartoon` must also be set to `True`.
        Default: "spectrum"
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorSpec

    fourth : optional
        `style`

        `str` indicating a representation style of heavy atoms, choosing from either "line", "cross", "stick", or "sphere".
        Default: "stick"

    fifth : optional
        `colorscheme`
        
        `str` indicating the color scheme for heavy atoms represented by the `style` option. Options include:
            A lower-case standard color followed by "Carbon" (e.g. "orangeCarbon")
            "ssPyMOL": PyMol secondary colorscheme
            "ssJmol": Jmol secondary colorscheme
            "Jmol": Jmol primary colorscheme
            "default": default colorscheme
            "amino": amino acid colorscheme
            "shapely": shapely protien colorscheme
            "nucleic": nucleic acid colorscheme
            "chain": standard chain colorscheme
            "chainHetatm": chain Hetatm colorscheme
        Default: "blackCarbon"
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec
        
    sixth : optional
        `radius`

        `float` or `int` indicating the radius of the heavy atoms represented by the `style` option.
        Default: 0.1

    seventh : optional
        `label`

        `True` or `False` to show labels next to residues selected by the `residue_selector` option. 
        Default: True

    eighth : optional
        `label_fontsize`

        `int` or `float` indicating the font size of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: 12

    ninth : optional
        `label_background`

        `True` or `False` to show the background of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: False

    tenth : optional
        `label_fontcolor`

        `str` indicating a standard font color (e.g. "grey") for label text next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: "black"

    eleventh : optional
        `command`
        
        `dict` or `tuple` of `dict`s of `py3Dmol.view.setStyle()` commands. If specified, this option overrides all other options.
        Default: None
        Example:
            command = {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "greyCarbon", "radius": 0.15}}
            view = viewer.init(poses) + viewer.setStyle(command=command)
            view.show()

    Returns
    -------
    A Viewer instance.
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
            if isinstance(self.command, tuple):
                viewer.setStyle(*self.command)
            elif isinstance(self.command, dict):
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
    Show the specified surface for each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `residue_selector`
        
        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` to select residues
        on which to apply the surface.
        Default: None

    second : optional
        `surface_type`

        `str` indicating surface type to be displayed. py3Dmol supports the following options:
            "VDW": Van der Waals surface
            "MS": Molecular surface
            "SES": Solvent excluded surface
            "SAS": Solvent accessible surface
        Default: "VDW"

    third : optional
        `opacity`

        `float` or `int` between 0 and 1 for opacity of the displayed surface.
        Default: 0.5

    fourth : optional
        `color`

        `str` indicating a standard color (e.g. "grey") of the surface to be displayed. 
        Either `color` or `colorscheme` may be specified, where `colorscheme` overrides `color`.
        Default: None

    fifth : optional
        `colorscheme`

        `str` indicating the color scheme of the surface to be displayed.
        Either `color` or `colorscheme` may be specified, where `colorscheme` overrides `color`.
        Options include:
            A lower-case standard color followed by "Carbon" (e.g. "yellowCarbon")
            "ssPyMOL": PyMol secondary colorscheme
            "ssJmol": Jmol secondary colorscheme
            "Jmol": Jmol primary colorscheme
            "default": default colorscheme
            "amino": amino acid colorscheme
            "shapely": shapely protien colorscheme
            "nucleic": nucleic acid colorscheme
            "chain": standard chain colorscheme
            "chainHetatm": chain Hetatm colorscheme
        Default: None
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec
        
    Returns
    -------
    A Viewer instance.
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
    Set the zoom magnification factor of each initialized `.pdb` file, `Pose` or `PackedPose` object.
    Values >1 zoom in, and values <1 zoom out.
    
    Parameters
    ----------
    first : optional
        `factor`

        `float` or `int` indicating the zoom magnification factor.
        Default: 2

    Returns
    -------
    A Viewer instance.
    """
    def __init__(self, factor=2):

        self.factor = factor
    
    def apply(self, viewer, pose, pdbstring):

        viewer.zoom(self.factor)

        return viewer


class setZoomTo:
    """
    Zoom to a `ResidueSelector` in each initialized `.pdb` file, `Pose` or `PackedPose` object.
    
    Parameters
    ----------
    first : optional
        `residue_selector`

        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` into which to zoom.
        Default: None

    Returns
    -------
    A Viewer instance.
    """
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
    """Exception raised for errors with the input argument `residue_selector`."""
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
    Given a `Pose` object and `ResidueSelector` object, return a `tuple` of `list`s containing
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
    """Convert pdbstring to a `Pose` with logging."""
    logger.info(
        " ".join("{0} requires pyrosetta.rosetta.core.pose.Pose object but given input .pdb file. \
        Now instantiating pyrosetta.rosetta.core.pose.Pose object from input .pdb file. \
        For faster performance, either input pyrosetta.rosetta.core.pose.Pose \
        or pyrosetta.distributed.packed_pose.core.PackedPose objects to pyrosetta.distributed.viewer.init, \
        or do not add {0} objects that require a pyrosetta.rosetta.core.pose.Pose object.  \
        ".format(class_name).split())
    )
    return io.to_pose(io.pose_from_pdbstring(pdbstring))
