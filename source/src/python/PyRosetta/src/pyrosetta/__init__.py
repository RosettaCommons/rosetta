from __future__ import absolute_import
# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

###############################################################################
# Imports.
# Standard library.
import os, sys, platform, os.path, json, random, gzip

import pyrosetta.rosetta as rosetta
import pyrosetta.method_bindings

import warnings
import logging
logger = logging.getLogger("rosetta")

import pyrosetta.logging_support as logging_support

# PyRosetta-3 comapatability
# WARNING WARNING WARNING: do not add anything extra imports/names here! If you feel strongly that something needs to be added please contact author first!
from pyrosetta.rosetta.core.kinematics import FoldTree, MoveMap
from pyrosetta.rosetta.core.io.pdb import dump_pdb
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring import ScoreFunction

from pyrosetta.rosetta.protocols.moves import PyMOLMover, SequenceMover, RepeatMover, TrialMover, MonteCarlo
from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn

from pyrosetta.io import *

get_score_function = rosetta.core.scoring.get_score_function
create_score_function = rosetta.core.scoring.ScoreFunctionFactory.create_score_function

rosetta.utility.vector1_string = rosetta.utility.vector1_std_string


###############################################################################
# Constants and globals

# FIXME: create 'version' struct in utility instead
rosetta_version =  rosetta.utility.Version.package() + ' ' + rosetta.utility.Version.version() + rosetta.utility.Version.commit()

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
class PyRosettaException(Exception):
    def __str__(self):
        return 'PyRosettaException'


class PythonPyExitCallback(rosetta.utility.py.PyExitCallback):
    def __init__(self):
        rosetta.utility.py.PyExitCallback.__init__(self)

    def exit_callback(self):
        raise PyRosettaException()


###############################################################################
#
def _rosetta_database_from_env():
    """Read rosetta database directory from environment or standard install locations.

    Database resolution proceeds by first searching the current installation for a 'database' or 'rosetta_database'
    path. If not found the search then continues to the users's home dir, cygwin, and osx standard installation
    locations.

    Returns database path if found, else None."""

    # Figure out database dir....
    if 'PYROSETTA_DATABASE' in os.environ:
        database = os.path.abspath(os.environ['PYROSETTA_DATABASE'])
        if os.path.isdir(database):
            logger.info('PYROSETTA_DATABASE environment variable was set to: %s; using it....', database)
            return database
        else:
            logger.warning('Invalid PYROSETTA_DATABASE environment variable was specified: %s', database)

    candidate_paths = []

    database_names = ["rosetta_database", "database"]
    for database_name in database_names:
        #Package directory database
        candidate_paths.append(os.path.join(os.path.dirname(__file__), database_name))
        candidate_paths.append(os.path.join(os.path.dirname(__file__), "..", database_name))

    for database_name in database_names:
        #Current directory database
        candidate_paths.append(database_name)

        #Home directory database
        if 'HOME' in os.environ:
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


# this variable and function are deprecated becuase they are not used -- JRJ 2017.17.01
# use instead rosetta.basic.database.full_name("path/to/db/file")
#_ROSETTA_DATABASE_PATH_ = None
#def get_rosetta_database_path(): return _ROSETTA_DATABASE_PATH_

def _is_interactive():
    """Determine if in an interactive context.

    See: https://stackoverflow.com/questions/2356399/tell-if-python-is-in-interactive-mode
    """

    import __main__ as main
    return not hasattr(main, '__file__')

def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=None, notebook=None):
    """Initialize Rosetta.  Includes core data and global options.

    options string with default Rosetta command-line options args.
            (default: '-ex1 -ex2aro')
    kargs -
        extra_options - Extra command line options to pass rosetta init.
                        (default None)
        set_logging_handler - Route rosetta tracing through logging logger 'rosetta':
            None - Set handler if interactive, otherwise not.
            False - Write logs via c++-level filehandles.
            "interactive" - Register python log handling and make visible if not.
            "logging" - Register python log handling, do not update logging config.
            True - Register python log handling, make visible if logging isn't configured.

    Examples:
        init()                     # uses default flags
        init(extra_options='-pH')  # adds flags to supplement the default
        init('-pH -database /home/me/pyrosetta/rosetta_database')  # overrides default flags - be sure to include the dB last
    """

    if set_logging_handler is None and _is_interactive():
        set_logging_handler = "interactive"
    elif notebook is not None:
        warnings.warn(
            "pyrosetta.init 'notebook' argument is deprecated and may be removed in 2018. "
            "See set_logging_handler='interactive'.",
            stacklevel=2
        )
        set_logging_handler = "interactive"


    assert set_logging_handler in (None, True, False, "interactive", "logging")

    logging_support.maybe_initialize_handler(set_logging_handler)
    if (set_logging_handler):
        logging_support.set_logging_sink()

    args = ['PyRosetta'] + options.split() + extra_options.split()

    # Attempt to resolve database location from environment if not present, else fallback
    # to rosetta's standard resolution
    if not "-database" in args:
        database = _rosetta_database_from_env()
        if database is not None: args.extend(["-database", database])

    # this variable is deprecated becuase it is not used -- JRJ 2017.17.01
    # use instead rosetta.basic.database.full_name("path/to/db/file")
    #_ROSETTA_DATABASE_PATH_ = database

    v = rosetta.utility.vector1_string()
    v.extend(args)

    logger.info(version())
    rosetta.protocols.init.init(v)


def version():
    return "PyRosetta-4 2017 [Rosetta " + rosetta_version + ' ' + rosetta.utility.Version.date() + \
           "] retrieved from: " + rosetta.utility.Version.url() + \
           "\n(C) Copyright Rosetta Commons Member Institutions." + \
           "\nCreated in JHU by Sergey Lyskov and PyRosetta Team.\n"


# if config['core']:
#     Pose.__iter__ = _Pose_residue_iterator
#     '''
#     # Add get() method to Pose (for compatibility with Rosetta Mover containers).
#     def _get(self):
#         """
#         Returns the Pose itself.
#         The entire purpose of this method is so that custom Movers in PyRosetta do
#         not crash when apply(pose) is called from within Mover containers (such as
#         TrialMover, RepeatMover, etc.)  When these C++ Movers call the apply()
#         method of the custom PyRosetta Mover, they send a PoseAP, which Python can-
#         not handle without converting to a raw pointer with get().  The addition of
#         this get() method to Pose allows for one to always use get() in the apply()
#         method of the custom Mover without concern for whether the apply() was
#         called from within Python code or C++ code.
#         """
#         return self
#     '''
#     Pose.get = lambda x: x

# Vector compatibility: Adding 'extend' to all utility.vector* functions
def _vector_extend_func(vec, othervec):
    for i in othervec: vec.append(i)
for k, vectype in rosetta.utility.__dict__.items():
    if k.startswith("vector1_") or k.startswith("vector0_") or k.startswith("vectorL_"): vectype.extend = _vector_extend_func


# def new_vector1_init(self, arg1=None, arg2=False):
#     self.__old_init()
#     if hasattr(arg1, "__iter__"):
#         self.extend(arg1)
#     elif isinstance(arg1, type(1)):
#         for i in xrange(arg1):
#             self.append(arg2)


# def replace_init(cls, init):
#   cls.__old_init = cls.__init__
#   cls.__init__ = init


def Vector1(list_in):
    """Creates a Vector1 object, deducing type from the given list."""

    if all([isinstance(x, bool) for x in list_in]):
        t = rosetta.utility.vector1_bool
    elif all([isinstance(x, int) for x in list_in]):
        t = rosetta.utility.vector1_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = rosetta.utility.vector1_double
    elif all([isinstance(x, str) for x in list_in]):
        t = rosetta.utility.vector1_string
    elif all([isinstance(x, rosetta.core.id.AtomID) for x in list_in]):
        t = rosetta.utility.vector1_AtomID
    else:
        raise Exception('Vector1: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    v = t()
    for i in list_in:
        v.append(i)
    return v


def Set(list_in):
    """Creates a std::set object, deducing type from the given list."""
    if all([isinstance(x, int) for x in list_in]):
        t = rosetta.utility.set_int
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = rosetta.utility.set_double
    elif all([isinstance(x, str) for x in list_in]):
        t = rosetta.utility.set_string
    else:
        raise Exception('Set: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    s = t()
    for i in list_in: s.add(i)
    return s


# New methods.
def generate_nonstandard_residue_set(pose, params_list):
    """
    Places the ResidueTypes corresponding to a list of .params filenames into a given pose

    .params files must be generated beforehand. Typically, one would obtain a
    molfile (.mdl) generated from the xyz coordinates of a residue, small
    molecule, or ion.  The script molfile_to_params.py can be used to convert
    to a Rosetta-readable .params file.  It can be found in the /test/tools
    folder of your PyRosetta installation or downloaded from the Rosetta
    Commons.

    Example:
        params = ["penicillin.params", "amoxicillin.params"]
        pose = Pose()
        generate_nonstandard_residue_set(pose, params)
        pose_from_file(pose, "TEM-1_with_substrates.pdb")
    See also:
        ResidueTypeSet
        Vector1()
        pose_from_file()
    """
    res_set = pose.conformation().modifiable_residue_type_set_for_conf()
    res_set.read_files_for_base_residue_types(Vector1(params_list))
    pose.conformation().reset_residue_type_set_for_conf( res_set )
    return pose.residue_type_set_for_pose()

def standard_task_factory():
        tf = rosetta.core.pack.task.TaskFactory()
        tf.push_back(rosetta.core.pack.task.operation.InitializeFromCommandline())
        #tf.push_back(rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(rosetta.core.pack.task.operation.NoRepackDisulfides())
        return tf


def standard_packer_task(pose):
        tf = standard_task_factory()
        task = tf.create_task_and_apply_taskoperations(pose)
        return task


# def add_extra_options():
#     rosetta.protocols.abinitio.AbrelaxApplication.register_options()
#     rosetta.protocols.abinitio.IterativeAbrelax.register_options()
#     rosetta.protocols.abinitio.register_options_broker()


# By Michael Pacella
def etable_atom_pair_energies(atom1, atom2, sfxn):
    """
    Usage: lj_atr, lj_rep, solv=etable_atom_pair_energies(atom1, atom2, sfxn)
        Description: given a pair of atoms and scorefunction, use the precomputed
        'etable' to return LJ attractive, LJ repulsive, and LK solvation energies
    """
    score_manager = rosetta.core.scoring.ScoringManager.get_instance()
    etable_ptr = score_manager.etable( sfxn.energy_method_options().etable_type() )
    etable = etable_ptr.lock()
    etable_energy = rosetta.core.scoring.etable.AnalyticEtableEnergy(etable,
                                                  sfxn.energy_method_options())

        # Construct AtomPairEnergy container to hold computed energies.
    ape = rosetta.core.scoring.etable.AtomPairEnergy()

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

def output_scorefile(pose, pdb_name, current_name, scorefilepath, \
                 scorefxn, nstruct, native_pose=None, additional_decoy_info=None):
    """
    Moved from PyJobDistributor (Jared Adolf-Bryfogle)
    Creates a scorefile if none exists, or appends the current one.
    Calculates and writes CA_rmsd if native pose is given,
    as well as any additional decoy info
    """
    if not os.path.exists(scorefilepath):
        with open(scorefilepath, 'w') as f:
            f.write("pdb name: " + pdb_name + "     nstruct: " +
                    str(nstruct) + '\n')

    score = scorefxn(pose)	 # Calculates total score.
    score_line = pose.energies().total_energies().weighted_string_of(scorefxn.weights())
    output_line = "filename: " + current_name + " total_score: " + str(round(score, 2))

    # Calculates rmsd if native pose is defined.
    if native_pose:
        rmsd = rosetta.core.scoring.CA_rmsd(native_pose, pose)
        output_line = output_line + " rmsd: " + str(round(rmsd, 2))

    with open(scorefilepath, 'a') as f:
        if additional_decoy_info:
            f.write(output_line + ' ' + score_line + ' '+additional_decoy_info + '\n')
        else:
            f.write(output_line + ' ' + score_line + '\n')


class PyJobDistributor:
    def __init__(self, pdb_name, nstruct, scorefxn, compress=False):
        self.pdb_name = pdb_name
        self.nstruct = nstruct
        self.compress = compress

        self.current_id = None
        self.current_name = None	      # Current decoy name
        self.scorefxn = scorefxn	      # Used for final score calculation
        self.native_pose = None		      # Used for rmsd calculation
        self.additional_decoy_info = None     # Used for any additional decoy information you want stored

        self.sequence = list(range(nstruct));  random.shuffle(self.sequence)
        self.start_decoy()		      # Initializes the job distributor


    @property
    def job_complete(self): return len(self.sequence) == 0

    def start_decoy(self):
        while(self.sequence):
            self.current_id = self.sequence[0]

            self.current_name = self.pdb_name + '_' + str(self.current_id) + ('.pdb.gz' if self.compress else '.pdb')
            self.current_in_progress_name = self.current_name + '.in_progress'

            if (not os.path.isfile(self.current_name))  and  (not os.path.isfile(self.current_in_progress_name)):
                with open(self.current_in_progress_name, 'w') as f: f.write("This decoy is in progress.")
                print( 'Working on decoy: {}'.format(self.current_name) )
                break


    def output_decoy(self, pose):
        if os.path.isfile(self.current_name): # decoy is already exist, probably written to other process -> moving to next decoy if any

            if os.path.isfile(self.current_in_progress_name): os.remove(self.current_in_progress_name)

            if not self.job_complete:
                self.start_decoy()
                self.output_decoy(self, pose)

        else:
            if self.compress:
                s = pyrosetta.rosetta.std.ostringstream()
                pose.dump_pdb(s)

                z = pyrosetta.rosetta.utility.io.ozstream(self.current_name)
                z.write(s.str(), len(s.str() ) )
                del z

                #s = pyrosetta.rosetta.std.ostringstream()
                #pose.dump_pdb(s)
                #with gzip.open(self.current_name, 'wb') as f: f.write( s.str() )

            else:
                dump_pdb(pose, self.current_name)

            score_tag = '.fasc' if pose.is_fullatom() else '.sc'

            scorefile = self.pdb_name + score_tag
            output_scorefile(pose, self.pdb_name, self.current_name, scorefile, self.scorefxn,
                             self.nstruct, self.native_pose, self.additional_decoy_info)

            self.sequence.remove(self.current_id)

            if os.path.isfile(self.current_in_progress_name): os.remove(self.current_in_progress_name)

            self.start_decoy()


###############################################################################
# Decorator generation for custom PyRosetta energy methods.
_mem_EnergyMethods_ = []
_mem_EnergyCreators_ = []


class CD:
    '''Class to represent named tuples.'''
    def __init__(self, **entries):
        self.__dict__.update(entries)

    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'):
                r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2] + '|'


_ScoreTypesRegistryByType_ = [
    CD(base=rosetta.core.scoring.methods.ContextIndependentTwoBodyEnergy,
       first=rosetta.core.scoring.PyRosettaTwoBodyContextIndepenedentEnergy_first,
       last=rosetta.core.scoring.PyRosettaTwoBodyContextIndepenedentEnergy_last,
       methods={}),
    CD(base=rosetta.core.scoring.methods.ContextDependentTwoBodyEnergy,
       first=rosetta.core.scoring.PyRosettaTwoBodyContextDependentEnergy_first,
       last=rosetta.core.scoring.PyRosettaTwoBodyContextDependentEnergy_last,
       methods={}),
    CD(base=None,
       first=rosetta.core.scoring.PyRosettaEnergy_first,
       last=rosetta.core.scoring.PyRosettaEnergy_last,
       methods={}),
]

ScoreTypesRegistry = {}


def defineEnergyMethodCreator(class_, scoreType):
    class Abstract_EnergyMethodCreator(
                             rosetta.core.scoring.methods.EnergyMethodCreator):
        def __init__(self):
            rosetta.core.scoring.methods.EnergyMethodCreator.__init__(self)

        def create_energy_method(self, energy_method_options):
            e = self.EnergyMethodClass()
            _mem_EnergyMethods_.append(e)
            return e

        def score_types_for_method(self):
            sts = rosetta.utility.vector1_core_scoring_ScoreType()
            sts.append(self.scoreType)
            return sts

    class_name = class_.__name__ + '_Creator'
    new_class = type(class_name, (Abstract_EnergyMethodCreator,),
                     {'EnergyMethodClass': class_,
                      'scoreType': rosetta.core.scoring.ScoreType(scoreType)})
    #globals()[class_name ] = new_class

    return new_class


class EnergyMethod:
    """
    Decorator function for custom EnergyMethods in PyRosetta.
    """
    def __init__(self, scoreName=None, scoreType=None, version=1):
        self.scoreName = scoreName
        self.scoreType = scoreType
        self.version = version

    def __call__(self, original_class):
        self.scoreName = self.scoreName or original_class.__name__
        # Try to automatically determine first avaliable scoreType.
        if not self.scoreType:
            for s in _ScoreTypesRegistryByType_:
                if not s.base or issubclass(original_class, s.base):
                    self.scoreType = max(s.methods.keys() or [int(s.first) - 1]) + 1
                    if self.scoreType > int(s.last):
                        err_msg = 'Cannot find free ScoreType to create %s! (looking in range [%s, %s])' % (self.scoreName, s.first, s.last)
                        raise Exception(err_msg)
                    s.methods[self.scoreType] = self.scoreName
                    ScoreTypesRegistry[self.scoreType] = self.scoreName
                    break

        def _clone(self):
            _mem_EnergyMethods_.append( self.__class__() )
            return _mem_EnergyMethods_[-1]
            #return type(self)()

        def _f_version(self):
            return self.version

        def _indicate_required_context_graphs(self, v):
            pass

        creator = defineEnergyMethodCreator(original_class, self.scoreType)

        if 'clone' not in original_class.__dict__:
            original_class.clone = _clone
        if 'version' not in original_class.__dict__:
            original_class.version = _f_version
        if 'indicate_required_context_graphs' not in original_class.__dict__:
            original_class.indicate_required_context_graphs = _indicate_required_context_graphs

        original_class.creator = creator
        original_class.scoreType = rosetta.core.scoring.ScoreType(self.scoreType)

        _mem_EnergyCreators_.append( creator() )
        rosetta.core.scoring.methods.PyEnergyMethodRegistrator(_mem_EnergyCreators_[-1])

        return original_class

def rrange(n):
    """Return an iterator from 1 to n inclusive"""
    # think of this as "Rosetta range"
    assert(type(n) == int)
    return range(1, n + 1)
