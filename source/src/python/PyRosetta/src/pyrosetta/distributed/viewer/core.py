# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.


import functools
import logging
import os
import pyrosetta
import pyrosetta.distributed.io as io
import sys
import time

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose


_logger = logging.getLogger("pyrosetta.distributed.viewer")

try:
    from IPython.core.display import display, HTML
    from IPython.display import clear_output
except ImportError:
    _logger.error("IPython.core.display or IPython.display module cannot be imported.")

try:
    import numpy
    import py3Dmol
    from ipywidgets import interact, IntSlider
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.viewer' requires the third-party packages " \
        + "'numpy', 'py3Dmol', and 'ipywidgets' as dependencies!\n" \
        + "Please install these packages into your python environment. " \
        + "For installation instructions, visit:\n" \
        + "https://pypi.org/project/numpy/\n" \
        + "https://pypi.org/project/py3Dmol/\n" \
        + "https://ipywidgets.readthedocs.io/en/latest/user_install.html"
    )
    raise


class Viewer:
 
    def __init__(self, poses, pdbstrings, window_size, modules, delay, continuous_update, *args, **kwargs):

        self.poses = poses
        self.pdbstrings = pdbstrings
        self.window_size = window_size
        self.modules = modules
        self.delay = delay
        self.continuous_update = continuous_update
        self._toggle_window(self.window_size)
        self._toggle_scrolling()

    def __add__(self, other, *args, **kwargs):
        
        self.modules += [other]
        
        return Viewer(
            poses=self.poses,
            pdbstrings=self.pdbstrings,
            window_size=self.window_size,
            modules=self.modules,
            delay=self.delay,
            continuous_update=self.continuous_update,
            *args, **kwargs
        )

    def __radd__(self, other):

        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __call__(self):

        return self.show()

    def _clear_output(self):
        
        try:
            _logger.debug("IPython.display clearing Jupyter notebook cell output.")
            clear_output(wait=True)
        except NameError as e:
            _logger.debug(e)

    def _toggle_scrolling(self):
        
        try:
            _logger.debug("IPython.core.display toggling scrolling in Jupyter notebook cell.")
            display(HTML("<script>$('.output_scroll').removeClass('output_scroll')</script>"))
        except NameError as e:
            _logger.debug(e)

    def _toggle_window(self, _window_size):
        
        try:
            _logger.debug("IPython.core.display toggling cell window area in Jupyter notebook.")
            HTML("""<style>
                    .output_wrapper, .output {
                        height:auto !important;
                        max-height:%ipx;
                    }
                    .output_scroll {
                        box-shadow:none !important;
                        webkit-box-shadow:none !important;
                    }
                    </style>
            """ % numpy.ceil(_window_size[1])
            )
        except NameError as e:
            _logger.debug(e)

    def add(self, other):
        """Add a module to the Viewer instance."""
        return self.__add__(other)

    def reinit(self):
        """Subtract all modules from the Viewer instance."""
        self.modules = []

    def reset(self):
        """Delete Viewer instance attributes."""
        self.poses = None
        self.pdbstrings = None
        self.window_size = None
        self.modules = None
        self.delay = None
        self.continuous_update = None

    def show(self):
        """Display Viewer in Jupyter notebook."""
        def view(i=0):
            
            _viewer = py3Dmol.view(*self.window_size)
            _pose = self.poses[i]
            _pdbstring = self.pdbstrings[i]

            if _pose:
                _viewer.addModels(io.to_pdbstring(_pose), "pdb")
            else:
                _viewer.addModels(_pdbstring, "pdb")
            _viewer.zoomTo()

            for module in self.modules:
                _viewer = module.apply(viewer=_viewer, pose=_pose, pdbstring=_pdbstring)

            self._clear_output()
            
            if _pose:
                _logger.debug("Decoy {0}: {1}".format(i, _pose.pdb_info().name()))

            return _viewer.show()
        
        time.sleep(self.delay)

        num_decoys = len(self.pdbstrings)
        if num_decoys > 1:
            s_widget = IntSlider(
                min=0, max=num_decoys - 1, description="Decoys", continuous_update=self.continuous_update
            )
            widget = interact(view, i=s_widget)
        else:
            widget = view()

        return widget


class ViewerInputError(Exception):
    """Exception raised for errors with the input argument 'packed_and_poses_and_pdbs'."""
    def __init__(self, obj):

        super().__init__(
            " ".join(
                "Input argument 'packed_and_poses_and_pdbs' should be an instance of \
                pyrosetta.rosetta.core.pose.Pose, pyrosetta.distributed.packed_pose.core.PackedPose, \
                or a valid path string to a .pdb file, or a list, set, or tuple of these objects. \
                Input argument 'packed_and_poses_and_pdbs' was invoked with: {0}".format(obj).split()
            )
        )


def init(packed_and_poses_and_pdbs=None, window_size=None, modules=None,
         delay=None, continuous_update=None, *args, **kwargs):
    """
    Initialize the Viewer object.

    Parameters
    ----------
    first : required
        `packed_and_poses_and_pdbs`

        `PackedPose`, `Pose`, or `str` of a valid path to a .pdb file, or a `list`, `set`, or `tuple` of these objects.

    second : optional
        `window_size`

        `list` or `tuple` of `int` or `float` values for the (width, height) dimensions of the displayed window screen size.
        Default: (1200, 800)

    third : optional
        `modules`
        
        `list` of instantiated visualization modules to run upon changing amongst `packed_and_poses_and_pdbs` objects
        with the slider, matching the namespace pyrosetta.distributed.viewer.set*
        Default: []

    fourth : optional
        `delay`
        
        `float` time delay in seconds before rendering the Viewer in a Jupyter notebook, which is useful to prevent
        overburdening the Jupyter notebook client if `for` looping over quick modifications to a `Pose`.
        Default: 0.25

    fifth : optional
        `continuous_update`
        
        `True` or `False`. When using the interactive slider widget, `False` restricts rendering to mouse release events.
        Default: False

    Returns
    -------
    A Viewer instance.
    """

    _default_window_size = (1200, 800)
    _default_modules = []
    _default_delay = 0.25 # seconds
    _default_continuous_update = False

    @functools.singledispatch
    def to_pose(obj):
        raise ViewerInputError(obj)

    to_pose.register(type(None), lambda obj: None)
    to_pose.register(PackedPose, lambda obj: io.to_pose(obj))
    to_pose.register(Pose, lambda obj: obj)
    to_pose.register(str, lambda obj: None)
    
    @functools.singledispatch
    def to_pdbstring(obj):
        raise ViewerInputError(obj)

    to_pdbstring.register(type(None))
    def _(obj):
        raise ViewerInputError(obj)

    to_pdbstring.register(PackedPose, lambda obj: io.to_pdbstring(obj))
    to_pdbstring.register(Pose, lambda obj: io.to_pdbstring(obj))
        
    @to_pdbstring.register(str)
    def _(obj):
        if not os.path.isfile(obj):
            raise ViewerInputError(obj)
        else:
            with open(obj, "r") as f:
                return f.read()

    if isinstance(packed_and_poses_and_pdbs, (list, set, tuple)):
        poses, pdbstrings = map(
            list, zip(*[(to_pose(p), to_pdbstring(p)) for p in packed_and_poses_and_pdbs])
        )
    else:
        poses = [to_pose(packed_and_poses_and_pdbs)]
        pdbstrings = [to_pdbstring(packed_and_poses_and_pdbs)]

    @functools.singledispatch
    def to_window_size(obj):
        _logger.warning(
            "Input argument 'window_size' cannot be parsed. Setting 'window_size' to default."
        )
        return _default_window_size
    
    to_window_size.register(type(None), lambda obj: _default_window_size)

    @to_window_size.register(tuple)
    @to_window_size.register(list)
    def _(obj):
        assert len(obj) == 2, "Input argument 'window_size' must be a list or tuple of length 2."
        return obj

    window_size = to_window_size(window_size)

    if not modules:
        modules = _default_modules
    assert isinstance(modules, list), "Input argument 'modules' should be an instance of list."

    @functools.singledispatch
    def to_delay(obj):
        _logger.warning(
            "Input argument 'delay' should be an instance of float. Setting 'delay' to default."
        )
        return _default_delay

    to_delay.register(type(None), lambda obj: _default_delay)
    to_delay.register(int, lambda obj: float(obj))

    @to_delay.register(str)
    def _(obj):
        try:
            _delay = float(obj)
        except ValueError:
            _logger.warning(
                "Input argument 'delay' cannot be parsed. Setting 'delay' to default."
            )
            _delay = _default_delay
        return _delay

    delay = to_delay(delay)

    if not continuous_update:
        continuous_update = _default_continuous_update
    assert type(continuous_update) == bool, "Input argument 'continuous_update' must be boolean."

    return Viewer(poses=poses, pdbstrings=pdbstrings, window_size=window_size, modules=modules,
        delay=delay, continuous_update=continuous_update, *args, **kwargs)

def expand_notebook():
    """Expand Jupyter notebook cell width."""
    try:
        _logger.debug("IPython.core.display expanding Jupyter notebook cell width.")
        display(HTML("<style>.container { width:100% !important; }</style>"))
    except NameError:
        _logger.exception("IPython.core.display module not imported.")
