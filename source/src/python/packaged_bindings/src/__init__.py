# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

from __future__ import absolute_import

###############################################################################
# Imports.

import logging
logger = logging.getLogger("rosetta")

# Rosetta
import warnings
warnings.filterwarnings("ignore", "to-Python converter for .+ already registered; second conversion method ignored.", RuntimeWarning, "^rosetta\\.")

#from rosetta.core.import_pose import pose_from_file
#from rosetta.core.io.pdb import dump_pdb

# Version and build information for module consumers
from .version import commit as __version__
def get_include():
    """Return the directory that contains the Rosetta header files.

    Extension modules that need to compile against librosetta should use this
    function to locate the appropriate include directory.

    Example:
        When using ``distutils``, for example in ``setup.py``.

        import rosetta
        ...
        Extension('extension_name', ...
                include_dirs=[rosetta.get_include()])
        ...
    """
    from os import path
    return os.path.abspath(path.join(__path__[0], "include"))

from .initialization import init

###############################################################################
# Modifications to Rosetta.
def pose_from_sequence(seq, res_type="fa_standard", auto_termini=True):
    """
    Returns a pose generated from a single-letter sequence of amino acid
    residues in <seq> using the <res_type> ResidueType and creates N- and C-
    termini if <auto_termini> is set to True.

    Unlike make_pose_from_sequence(), this method generates a default PDBInfo
    and sets all torsion angles to 180 degrees.

    Example:
        pose = pose_from_sequence("THANKSEVAN")
    See also:
        Pose
        make_pose_from_sequence()
        pose_from_file()
        pose_from_rcsb()
    """

    from .core.pose import Pose, PDBInfo
    from .core.pose import make_pose_from_sequence

    pose = Pose()
    make_pose_from_sequence(pose, seq, res_type, auto_termini)
    #print 'Setting phi, psi, omega...'
    for i in range(0, pose.total_residue()):
        pose.set_phi(i + 1, 180)
        pose.set_psi(i + 1, 180)
        pose.set_omega(i + 1, 180)
    #print 'Attaching PDBInfo...'
    # Empty PDBInfo (rosetta.core.pose.PDBInfo()) is not correct here;
    # we have to reserve space for atoms....
    pose.pdb_info(PDBInfo(pose))
    pose.pdb_info().name(seq[:8])
    #print pose
    return pose
