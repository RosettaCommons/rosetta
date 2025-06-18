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
import os, sys, shlex

import pyrosetta.rosetta as rosetta

# for backward compatibility
rosetta.utility.vector1_string = rosetta.utility.vector1_std_string

rosetta.std.list_std_shared_ptr_core_pack_task_operation_TaskOperation_std_allocator_std_shared_ptr_core_pack_task_operation_TaskOperation_t = rosetta.std.list_std_shared_ptr_core_pack_task_operation_TaskOperation_t
rosetta.std.vector_std_vector_double_std_allocator_double_t = rosetta.std.vector_std_vector_double_t
rosetta.std.set_std_string_std_less_std_string_std_allocator_std_string_t = rosetta.std.set_std_string_t

import pyrosetta.bindings
import pyrosetta.protocols

import warnings
import logging
logger = logging.getLogger("pyrosetta.rosetta")

import pyrosetta.logging_support as logging_support

# this try/except block is due to decorator-module requirement of the Distributed framework
try:
    from pyrosetta.distributed.utility.log import LoggingContext
except:
    pass

from pyrosetta.toolbox import etable_atom_pair_energies, PyJobDistributor

# PyRosetta-3 comapatability
# WARNING WARNING WARNING: do not add anything extra imports/names here! If you feel strongly that something needs to be added please contact author first!
from pyrosetta.rosetta.core.kinematics import FoldTree, MoveMap
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring import ScoreFunction

from pyrosetta.rosetta.protocols.moves import PyMOLMover, SequenceMover, RepeatMover, TrialMover, MonteCarlo
from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn

from pyrosetta.io import (
    pose_from_pdb,
    pose_from_file,
    poses_from_files,
    pose_from_sequence,
    poses_from_sequences,
    poses_from_silent,
    poses_from_multimodel_pdb,
    poses_to_silent,
    dump_file,
    dump_scored_pdb,
    dump_pdb,
    dump_multimodel_pdb,
    dump_cif,
    dump_mmtf,
    Pose,
)

from pyrosetta.rosetta.core.scoring import get_score_function
create_score_function = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function

###############################################################################
# Exception handling.
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

def _is_interactive():
    """Determine if in an interactive context.

    See: https://stackoverflow.com/questions/2356399/tell-if-python-is-in-interactive-mode
    """

    import __main__ as main
    return not hasattr(main, '__file__')

def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=None, notebook=None, silent=False):
    """Initialize Rosetta.  Includes core data and global options.

    options string with default Rosetta command-line options args.
            (default: '-ex1 -ex2aro')
    kwargs -
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

    if not isinstance(options, (str, list)) or not isinstance(extra_options, (str, list)):
        raise Exception("Command line rosetta arguments must be passed as list or string!")

    if (set_logging_handler is None) and _is_interactive():
        set_logging_handler = "interactive"
    elif notebook is not None:
        warnings.warn(
            "pyrosetta.init 'notebook' argument is deprecated and may be removed in a future release. "
            "See set_logging_handler='interactive'.",
            stacklevel=2
        )
        set_logging_handler = "interactive"


    assert set_logging_handler in (None, True, False, "interactive", "logging"), \
        "pyrosetta.init 'set_logging_handler' argument must be either: None, True, False, 'interactive', or 'logging'."

    logging_support.maybe_initialize_handler(set_logging_handler)
    if (set_logging_handler):
        logging_support.set_logging_sink()

    if isinstance(options, str):
        options = shlex.split(options)
    if isinstance(extra_options, str):
        extra_options = shlex.split(extra_options)

    args = ['PyRosetta'] + options + extra_options

    # Attempt to resolve database location from environment if not present, else fallback
    # to rosetta's standard resolution
    if not "-database" in args:
        database = _rosetta_database_from_env()
        if database is not None: args.extend(["-database", database])

    v = rosetta.utility.vector1_string()
    v.extend(args)

    if not silent:
        print( version() )
        logger.info( version() )
    else:
        logger.debug( version() )
    rosetta.protocols.init.init(v)
    pyrosetta.protocols.h5_fragment_store_provider.init_H5FragmentStoreProvider()
    pyrosetta.protocols.h5_structure_store_provider.init_H5StructureStoreProvider()

# FIXME: create 'version' struct in utility instead
def _version_string():
    version, commit = rosetta.utility.Version.version(), rosetta.utility.Version.commit()
    version = version.split(".")
    if commit.startswith(version[-1]):
        version.pop()
    version.append(commit)
    return rosetta.utility.Version.package() + " " + ".".join(version)


def version():
    return (
        '┌───────────────────────────────────────────────────────────────────────────────┐\n'
        '│                                  PyRosetta-4                                  │\n'
        '│               Created in JHU by Sergey Lyskov and PyRosetta Team              │\n'
        '│               (C) Copyright Rosetta Commons Member Institutions               │\n'
        '│                                                                               │\n'
        '│ NOTE: USE OF PyRosetta FOR COMMERCIAL PURPOSES REQUIRES PURCHASE OF A LICENSE │\n'
        '│          See LICENSE.PyRosetta.md or email license@uw.edu for details         │\n'
        '└───────────────────────────────────────────────────────────────────────────────┘\n'
    ) + \
    'PyRosetta-4 ' + rosetta.utility.Version.date().split('-').pop(0) + \
    ' [Rosetta ' + _version_string() + ' ' + rosetta.utility.Version.date() + \
    '] retrieved from: ' + rosetta.utility.Version.url()


###############################################################################
# Vector compatibility: Adding 'extend' to all utility.vector* functions
def _vector_extend_func(vec, othervec):
    for i in othervec: vec.append(i)
for k, vectype in rosetta.utility.__dict__.items():
    if k.startswith("vector1_") or k.startswith("vector0_") or k.startswith("vectorL_"): vectype.extend = _vector_extend_func


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
        t = rosetta.utility.vector1_core_id_AtomID
    else:
        raise Exception('Vector1: attemting to create vector of unknown type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    v = t()
    for i in list_in:
        v.append(i)
    return v


def Set(list_in):
    """Creates a std::set object, deducing type from the given list."""
    if all([isinstance(x, int) for x in list_in]):
        t = rosetta.std.set_int_t
    elif all([isinstance(x, float) or isinstance(x, int) for x in list_in]):
        t = rosetta.std.set_double_t
    elif all([isinstance(x, str) for x in list_in]):
        t = rosetta.std.set_std_string_t
    else:
        raise Exception('Set: attemting to create vector of unknow type ' +
                        'or mixed type vector init_list = ' + str(list_in))

    s = t()
    for i in list_in: s.add(i)
    return s

###############################################################################
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


###############################################################################
# Decorator generation for custom PyRosetta energy methods.
_mem_EnergyMethods_ = []
_mem_EnergyCreators_ = []

from collections import namedtuple
CD = namedtuple("CD", "base first last methods")

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
