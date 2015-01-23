# :noTabs=true:

"""
(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under
(c) license.
(c) The Rosetta software is developed by the contributing members of the
(c) Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org.
(c) Questions about this can be addressed to University of Washington UW
(c) TechTransfer, email: license@u.washington.edu.
"""
from __future__ import absolute_import

###############################################################################
# Imports.
# Standard library.
import os
import sys
import platform
import os.path

import logging
logger = logging.getLogger("rosetta")

# Rosetta
import warnings
warnings.filterwarnings("ignore", "to-Python converter for .+ already registered; second conversion method ignored.", RuntimeWarning, "^rosetta\\.")

# Double-checked right order...
try:
    __import__("rosetta.utility")
    __import__("rosetta.utility.excn")
    __import__("rosetta.utility.file")
except ImportError:
    pass

try:
    __import__("rosetta.numeric")
except ImportError:
    pass

try:
    __import__("rosetta.basic")
    __import__("rosetta.basic.datacache")
    __import__("rosetta.basic.resource_manager")
except ImportError:
    pass

try:
    __import__("rosetta.core")
    __import__("rosetta.core.graph")
    __import__("rosetta.core.chemical")
    __import__("rosetta.core.chemical.orbitals")
    __import__("rosetta.core.scoring")
    __import__("rosetta.core.scoring.methods")
    __import__("rosetta.core.scoring.constraints")
    __import__("rosetta.core.scoring.etable")
    __import__("rosetta.core.kinematics")

    __import__("rosetta.core.io.silent")
    __import__("rosetta.core.pose")

    __import__("rosetta.core.conformation")
    __import__("rosetta.core.id")
    __import__("rosetta.core.io")
    __import__("rosetta.core.io.pdb")
    __import__("rosetta.core.fragment")

    __import__("rosetta.core.pack")
    __import__("rosetta.core.pack.task")
    __import__("rosetta.core.scoring.hbonds")

    __import__("rosetta.core.pose.signals")
except ImportError:
    pass

try:
    __import__("rosetta.protocols")
    __import__("rosetta.protocols.moves")
    __import__("rosetta.protocols.jd2")
    __import__("rosetta.protocols.jd2.archive")
    __import__("rosetta.protocols.canonical_sampling")
    __import__("rosetta.protocols.simple_moves")
    __import__("rosetta.protocols.jumping")
    __import__("rosetta.protocols.abinitio")

    __import__("rosetta.protocols.filters")
    __import__("rosetta.protocols.docking")
    __import__("rosetta.protocols.init")

    __import__("rosetta.protocols.loops")
    __import__("rosetta.protocols.wum")
    __import__("rosetta.protocols.relax")
except ImportError:
    pass


from rosetta.core.import_pose import pose_from_pdb
from rosetta.core.io.pdb import dump_pdb
from rosetta.core.pose import make_pose_from_sequence

import rosetta.version as version

import rosetta.logging_support as logging_support

###############################################################################
# Constants and globals

__version__ =  version.commit_id + ':' + version.commit

# Create global '_PLATFORM' that will hold info of current system.
if sys.platform.startswith("linux"):
    _PLATFORM = "linux"  # can be linux1, linux2, etc.
elif sys.platform == "darwin":
    _PLATFORM = "macos"
elif sys.platform == "cygwin":
    _PLATFORM = "cygwin"
elif sys.platform == 'win32':
    _PLATFORM = 'windows'
else:
    _PLATFORM = "_unknown_"

#PlatformBits = platform.architecture()[0][:2]  # unused?

_python_py_exit_callback = None


###############################################################################
#Exception handling.
import rosetta.utility as utility
class PyRosettaException(Exception):
    #def __init__(self): pass

    def __str__(self):
        return 'PyRosettaException'


class PythonPyExitCallback(utility.PyExitCallback):
    def exit_callback(self):
        raise PyRosettaException()

    def __init__(self):
        utility.PyExitCallback.__init__(self)

###############################################################################
#
def rosetta_database_from_env():
    """Read rosetta database directory from environment or standard install locations.

    Returns database path if found, else None."""

    # Figure out database dir....
    if 'PYROSETTA_DATABASE' in os.environ:
        database = os.path.abspath(os.environ['PYROSETTA_DATABASE'])
        if os.path.isdir(database):
            logger.info('PYROSETTA_DATABASE environment variable was set to: %s; using it....', database)
            return database
        else:
            logger.warning('Invalid PYROSETTA_DATABASE environment variable was specified: %s', database)

    database_names = ["rosetta_database", "database"]

    for database_name in database_names:
        candidate_paths = []

        #Current directory database
        candidate_paths.append(database_name)

        #Package directory database
        candidate_paths.append(os.path.join(os.path.dirname(__file__), "..", database_name))

        #Home directory database
        candidate_paths.append(os.path.join(os.environ['HOME'], database_name))

        #Cygwin root install
        if sys.platform == "cygwin":
            candidate_paths.append(os.path.join('/', database_name))

        # Mac /usr/lib database install
        candidate_paths.append(os.path.join('rosetta', database_name))

        for candidate in candidate_paths:
            if os.path.isdir(candidate):
                database = os.path.abspath(candidate)
                logger.info('Found rosetta database at: %s; using it....', database)
                return database

    # No database found.
    return None

# rosetta.init()
def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=True):
    """Initialize Rosetta.  Includes core data and global options.

    options string with default Rosetta command-line options args.
            (default: '-ex1 -ex2aro')
    kargs -
        extra_options - Extra command line options to pass rosetta init.
                        (default None)
        set_logging_handler - Route rosetta tracing through logging logger 'rosetta'.
                        (default True)

    Examples:
        init()                     # uses default flags
        init(extra_options='-pH')  # adds flags to supplement the default
        init('-pH -database /home/me/pyrosetta/rosetta_database')  # overrides default flags - be sure to include the dB last
    """

    logging_support.initialize_logging()

    global _python_py_exit_callback
    _python_py_exit_callback = PythonPyExitCallback()
    utility.PyExitCallback.set_PyExitCallBack(_python_py_exit_callback)

    if set_logging_handler: logging_support.set_logging_handler()

    args = ['PyRosetta'] + options.split() + extra_options.split()

    # Attempt to resolve database location from environment if not present, else fallback
    # to rosetta's standard resolution
    if not "-database" in args:
        database = rosetta_database_from_env()
        if database is not None: args.extend(["-database", database])

    v = utility.vector1_string()
    v.extend(args)

    logger.info("Version: %s", __version__)

    try:
        from .protocols.init import init
    except ImportError:
        from .core.init import init

    init(v)


# MPI version of init function, use it instead of init(...)
def mpi_init(*args, **kargs):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    kargs['extra_options'] = kargs.get('extra_options', '') + ' -seed_offset %s' % (rank*10000)

    init(*args, **kargs)


def MPIJobDistributor(njobs, fun):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    myjobs = []

    if rank == 0:
        jobs = range(njobs)
        jobs.extend( [None]*(size - njobs % size) )
        n = len(jobs)/size
        for i in range(size):
            queue = []  # list of jobs for individual cpu
            for j in range(n):
                queue.append(jobs[j*size+i])

            if( i == 0 ):
                myjobs = queue
            else:
                # now sending the queue to the process
                logger.info('Sending %s to node %s' % (queue, i) )
                comm.send(queue, dest=i)
    else:
        # getting decoy lists
        myjobs = comm.recv(source=0)

    logger.info('Node %s, got queue:%s' % (rank, myjobs) )

    for j in myjobs:
        if j is not None: fun(j)

###############################################################################
# Modifications to Rosetta.
# Add iter property to Pose.

# Vector compatibility.
def _extendfunc(vec, othervec):
    for i in othervec:
        vec.append(i)


def _add_extend(vectype):
    vectype.extend = _extendfunc


for k, v in utility.__dict__.items():
    if k.startswith("vector1_"):
      _add_extend(v)


def new_vector1_init(self, arg1=None, arg2=False):
    self.__old_init()
    if hasattr(arg1, "__iter__"):
        self.extend(arg1)
    elif isinstance(arg1, type(1)):
        for i in xrange(arg1):
            self.append(arg2)


def replace_init(cls, init):
  cls.__old_init = cls.__init__
  cls.__init__ = init


def Vector1(list_in):
    """Creates a Vector1 object, deducing type from the given list."""

    if all([isinstance(x, bool) for x in list_in]):
        t = utility.vector1_bool
    elif all([isinstance(x, int) for x in list_in]):
        t = utility.vector1_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = utility.vector1_double
    elif all([isinstance(x, str) for x in list_in]):
        t = utility.vector1_string
    elif all([isinstance(x, core.id.AtomID) for x in list_in]):
        t = utility.vector1_AtomID
    else:
        raise Exception('Vector1: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    v = t()
    for i in list_in:
        v.append(i)
    return v


def Set(list_in):
    """Creates a Vector1 object, deducing type from the given list."""
    if all([isinstance(x, int) for x in list_in]):
        t = utility.set_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = utility.set_double
    elif all([isinstance(x, str) for x in list_in]):
        t = utility.set_string
    else:
        raise Exception('Set: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    s = t()
    for i in list_in: s.add(i)
    return s


# New methods.
def generate_nonstandard_residue_set(params_list):
    """
    Generates a ResidueTypeSet from a list of .params filenames.

    .params files must be generated beforehand. Typically, one would obtain a
    molfile (.mdl) generated from the xyz coordinates of a residue, small
    molecule, or ion.  The script molfile_to_params.py can be used to convert
    to a Rosetta-readable .params file.  It can be found in the /test/tools
    folder of your PyRosetta installation or downloaded from the Rosetta
    Commons.

    Example:
        params = ["penicillin.params", "amoxicillin.params"]
        type_set = generate_nonstandard_residue_set(params)
        pose = pose_from_pdb(type_set, "TEM-1_with_substrates.pdb")
    See also:
        ResidueTypeSet
        Vector1()
        pose_from_pdb()
    """
    res_set = ChemicalManager.get_instance().nonconst_residue_type_set("fa_standard")
    # res_set.read_files(Vector1(params_list),
    #                    ChemicalManager.get_instance().atom_type_set("fa_standard"),
    #                    ChemicalManager.get_instance().element_set('default'),
    #                    ChemicalManager.get_instance().mm_atom_type_set("fa_standard"),
    #                    ChemicalManager.get_instance().orbital_type_set("fa_standard"),)
    res_set.read_files(Vector1(params_list))
    return res_set


def standard_task_factory():
	tf = TaskFactory()
	tf.push_back(InitializeFromCommandline())
	tf.push_back(NoRepackDisulfides())
	return tf


def standard_packer_task(pose):
	tf = standard_task_factory()
	task = tf.create_task_and_apply_taskoperations(pose)
	return task

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
        pose_from_pdb()
        pose_from_rcsb()
    """
    from rosetta.core.pose import Pose
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
    pose.pdb_info(rosetta.core.pose.PDBInfo(pose))
    pose.pdb_info().name(seq[:8])
    #print pose
    return pose


# By Michael Pacella
def etable_atom_pair_energies(atom1, atom2, sfxn):
    """
    Usage: lj_atr, lj_rep, solv=etable_atom_pair_energies(atom1, atom2, sfxn)
	Description: given a pair of atoms and scorefunction, use the precomputed
	'etable' to return LJ attractive, LJ repulsive, and LK solvation energies
    """
    score_manager = core.scoring.ScoringManager.get_instance()
    etable_ptr = score_manager.etable(
                                    sfxn.energy_method_options().etable_type())
    etable=etable_ptr.get()
    etable_energy = core.scoring.etable.AnalyticEtableEnergy(etable,
                                                  sfxn.energy_method_options())

	# Construct AtomPairEnergy container to hold computed energies.
    ape = core.scoring.etable.AtomPairEnergy()

	# Set all energies in the AtomPairEnergy to zero prior to calculation.
    ape.attractive, ape.bead_bead_interaction, ape.repulsive, ape.solvation = \
                                                             0.0, 0.0, 0.0, 0.0

	# Calculate the distance squared and set it in the AtomPairEnergy.
    ape.distance_squared = atom1.xyz().distance_squared(atom2.xyz())

	# Evaluate energies from pre-calculated etable, using a weight of 1.0
	# in order to match the raw energies from eval_ci_2b.
    etable_energy.atom_pair_energy(atom1, atom2, 1.0, ape)

	# Calculate atom-atom scores.
    lj_atr = ape.attractive
    lj_rep = ape.repulsive
    solv = ape.solvation

    return lj_atr, lj_rep, solv
