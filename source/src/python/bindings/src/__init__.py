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

###############################################################################
# Imports.
# Standard library.
import os, sys, platform, os.path, json

import logging
logger = logging.getLogger("rosetta")

# Rosetta
import warnings
warnings.filterwarnings("ignore", "to-Python converter for .+ already registered; second conversion method ignored.", RuntimeWarning, "^rosetta\\.")

# bindings config
# print __name__
# print sys.modules[__name__].__file__
#config = json.load( file( os.path.join( os.path.split(__file__)[0], 'config.json' ) ) )
#config = json.load( file( os.path.join( os.path.split(__file__)[0], 'config.json' ) ) )
#config = {"low_memory_mode": False, "protocols": False, "core": False, "basic": False, "numeric": False, "utility": True, 'monolith': True}
#config = {"low_memory_mode": False, "protocols": False, "core": True, "basic": True, "numeric": True, "utility": True, 'monolith': True}
config = {"low_memory_mode": False, "protocols": True, "core": True, "basic": True, "numeric": True, "utility": True, 'monolith': True}
if '__file__' in vars():
    config_file_name = os.path.join( os.path.split(__file__)[0], 'config.json' )
    if os.path.isfile(config_file_name): config = json.load( file( config_file_name ) )

# import rosetta
# print dir(rosetta)
# print version
# import version
# import rosetta.version

# from version import *
# from rosetta.version import *

# print 'version.url', version.url
# print 'rosetta.version.url', rosetta.version.url

# if config['monolith']:
#     import _rosetta_ as rosetta
#     from _rosetta_ import *
#     #for m in ['utility', 'numeric', 'basic', 'core', 'protocols']:

#     def add_modules(module, parent_name):
#         for a in dir(module):
#             if hasattr( getattr(module, a), '__package__'):
#                 name = parent_name + '.' + a
#                 #print 'adding:', name
#                 sys.modules[name] = getattr(module, a)
#                 add_modules(getattr(module, a), name)

#     add_modules(_rosetta_, __name__) #if config[m]: sys.modules[__name__+'.' + m] = eval(m)  # emulate imports

#     # from rosetta.core.pose import Pose
#     # from rosetta.core.import_pose import pose_from_pdb
#     # from rosetta.core.io.pdb import dump_pdb
#     # from rosetta.core.pose import make_pose_from_sequence

#     # from rosetta.core.id import *
#     # from rosetta.core.conformation import *
#     # from rosetta.core.chemical import *
#     # from rosetta.core.scoring import *
#     # from rosetta.core.kinematics import *
#     # from rosetta.core.fragment import *
#     # from rosetta.core.pack.task import *
#     # from rosetta.core.pack.task.operation import *

#     # from rosetta.protocols.moves import *
#     # from rosetta.protocols.simple_moves import *
#     # from rosetta.protocols.abinitio import *
#     # from rosetta.protocols.docking import *
#     # from rosetta.protocols.loops import *
#     # from rosetta.protocols.relax import *

#     # from rosetta.protocols.simple_moves import *

#     # from rosetta.PyMolLink import *

# Double-checked right order...
import rosetta.utility
import rosetta.utility.py
import rosetta.utility.excn
import rosetta.utility.file

if config['monolith']:
    if config['basic']:
        import rosetta.basic.resource_manager
        import rosetta.basic.datacache
    if config['core']:
        import rosetta.core.scoring.func

if config['low_memory_mode']:
    import rosetta.core.pose
    import rosetta.core.graph
    import rosetta.core.scoring.methods
    import rosetta.protocols.init
    from rosetta.core.pose import Pose
    from rosetta.core.scoring import *

else:
    if config['numeric']:
        import rosetta.numeric
        import rosetta.numeric.geometry.hashing

    if config['basic']:
        import rosetta.basic
        import rosetta.basic.datacache
        import rosetta.basic.resource_manager


    if config['core']:
        import rosetta.core
        import rosetta.core.graph
        import rosetta.core.chemical
        import rosetta.core.chemical.orbitals
        import rosetta.core.scoring
        import rosetta.core.scoring.methods
        import rosetta.core.scoring.constraints
        import rosetta.core.scoring.etable
        import rosetta.core.kinematics

        import rosetta.core.io.silent
        import rosetta.core.pose

        import rosetta.core.graph
        import rosetta.core.conformation
        import rosetta.core.id
        import rosetta.core.io
        import rosetta.core.io.pdb
        import rosetta.core.fragment
        import rosetta.core.kinematics
        import rosetta.core.pack
        import rosetta.core.pack.task
        import rosetta.core.scoring.hbonds

        from rosetta.core.id import *
        from rosetta.core.conformation import *
        from rosetta.core.chemical import *
        from rosetta.core.pose import Pose

        from rosetta.core.import_pose import pose_from_pdb
        from rosetta.core.io.pdb import dump_pdb
        from rosetta.core.pose import make_pose_from_sequence

        import rosetta.core.pose
        from rosetta.core.scoring import *
        from rosetta.core.kinematics import *
        from rosetta.core.fragment import *
        from rosetta.core.pack.task import *
        from rosetta.core.pack.task.operation import *


    if config['protocols']:
        import rosetta.protocols
        import rosetta.protocols.moves
        import rosetta.protocols.canonical_sampling
        import rosetta.protocols.simple_moves
        import rosetta.protocols.jumping
        import rosetta.protocols.jd2
        import rosetta.protocols.jd2.archive
        import rosetta.protocols.abinitio

        import rosetta.protocols.filters
        import rosetta.protocols.docking
        import rosetta.protocols.init

        import rosetta.protocols.loops
        import rosetta.protocols.wum
        import rosetta.protocols.relax
        import rosetta.core.pose.signals

        import rosetta.protocols.simple_moves

        from rosetta.protocols.moves import *
        from rosetta.protocols.simple_moves import *
        from rosetta.protocols.abinitio import *
        from rosetta.protocols.docking import *
        from rosetta.protocols.loops import *
        from rosetta.protocols.relax import *

        #from rosetta.protocols.rigid import *
        from rosetta.protocols.simple_moves import *

    # PyMOLMover and associated methods.
    if config['protocols']: from rosetta.PyMolLink import *


import version

if config['basic']: import rosetta.logging_support as logging_support


###############################################################################
# Constants and globals

__version__ =  version.commit_id + ':' + version.commit  # rosetta.core.minirosetta_svn_version()

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
    #def __init__(self): pass

    def __str__(self):
        return 'PyRosettaException'


class PythonPyExitCallback(rosetta.utility.py.PyExitCallback):
    def __init__(self):
        rosetta.utility.py.PyExitCallback.__init__(self)

    def exit_callback(self):
        raise PyRosettaException()



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
        if '__file__' in vars(): candidate_paths.append(os.path.join(os.path.dirname(__file__), "..", database_name))

        #Home directory database
        if 'HOME' in os.environ: candidate_paths.append(os.path.join(os.environ['HOME'], database_name))

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


_ROSETTA_DATABASE_PATH_ = None
def get_rosetta_database_path(): return _ROSETTA_DATABASE_PATH_

# rosetta.init()
def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=True):
    import rosetta

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

    if config['basic']: logging_support.initialize_logging()

    #global _python_py_exit_callback
    _python_py_exit_callback = PythonPyExitCallback()
    rosetta.utility.py_xinc_ref(_python_py_exit_callback)
    rosetta.utility.py.PyExitCallback.set_PyExitCallBack(_python_py_exit_callback)
    #rosetta.utility.set_pyexit_callback()

    if set_logging_handler: logging_support.set_logging_handler()

    args = ['PyRosetta'] + options.split() + extra_options.split()

    # Attempt to resolve database location from environment if not present, else fallback
    # to rosetta's standard resolution
    if not "-database" in args:
        database = rosetta_database_from_env()
        if database is not None: args.extend(["-database", database])

    _ROSETTA_DATABASE_PATH_ = database

    v = rosetta.utility.vector1_string()
    v.extend(args)

    logger.info(version())
    if config['protocols']: rosetta.protocols.init.init(v)
    elif config['core']: import rosetta.core.init; rosetta.core.init.init(v)


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



def version():
    return "PyRosetta 2014 [Rosetta 2014 " + __version__ + \
           "] retrieved from: %s" % rosetta.core.minirosetta_svn_url() + \
           "\n(C) Copyright Rosetta Commons Member Institutions." + \
           "\nCreated in JHU by Sergey Lyskov and PyRosetta Team.\n"


###############################################################################
# Modifications to Rosetta.
# Add iter property to Pose.
def _Pose_residue_iterator(obj):
    def __pose_iter():
        for i in range(obj.total_residue()):
            yield obj.residue(i + 1)
    return __pose_iter()


if config['core']:
    Pose.__iter__ = _Pose_residue_iterator
    '''
    # Add get() method to Pose (for compatibility with Rosetta Mover containers).
    def _get(self):
        """
        Returns the Pose itself.

        The entire purpose of this method is so that custom Movers in PyRosetta do
        not crash when apply(pose) is called from within Mover containers (such as
        TrialMover, RepeatMover, etc.)  When these C++ Movers call the apply()
        method of the custom PyRosetta Mover, they send a PoseAP, which Python can-
        not handle without converting to a raw pointer with get().  The addition of
        this get() method to Pose allows for one to always use get() in the apply()
        method of the custom Mover without concern for whether the apply() was
        called from within Python code or C++ code.
        """
        return self
    '''
    Pose.get = lambda x: x


# Vector compatibility.
def _extend_func(vec, othervec):
    for i in othervec:
        vec.append(i)


def _add_extend(vectype):
    vectype.extend = _extend_func


for k, v in rosetta.utility.__dict__.items():
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
    """Creates a Vector1 object, deducing type from the given list."""
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
        #tf.push_back(IncludeCurrent())
        tf.push_back(NoRepackDisulfides())
        return tf


def standard_packer_task(pose):
        tf = standard_task_factory()
        task = tf.create_task_and_apply_taskoperations(pose)
        return task


def add_extra_options():
    rosetta.protocols.abinitio.AbrelaxApplication.register_options()
    rosetta.protocols.abinitio.IterativeAbrelax.register_options()
    rosetta.protocols.abinitio.register_options_broker()


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
    score_manager = rosetta.core.scoring.ScoringManager.get_instance()
    etable_ptr = score_manager.etable( sfxn.energy_method_options().etable_type() )
    etable = etable_ptr.get()
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
        rmsd = CA_rmsd(native_pose, pose)
        output_line = output_line + " rmsd: " + str(round(rmsd, 2))

    with open(scorefilepath, 'a') as f:
        if additional_decoy_info:
            f.write(output_line + ' ' + score_line + ' '+additional_decoy_info + '\n')
        else:
            f.write(output_line + ' ' + score_line + '\n')

# New classes.
class PyJobDistributor:
    def __init__(self, pdb_name, nstruct, scorefxn):
        self.pdb_name = pdb_name
        self.nstruct = nstruct
        self.current_num = 0		      # Current decoy number
        self.current_name = " "		      # Current decoy name
        self.job_complete = False	      # Job status
        self.scorefxn = scorefxn	      # Used for final score calculation
        self.native_pose = None		      # Used for rmsd calculation
        self.additional_decoy_info = None     # Used for any additional decoy
                                               # information you want stored
        self.start_decoy()		      # Initializes the job distributor

    def start_decoy(self):
        if self.job_complete:
            return
        i = 1
        file_exists = True
        while (file_exists and i <= self.nstruct):
            current_name = self.pdb_name + "_" + str(i) + ".pdb"
            if not os.path.exists(current_name):
                current_name_temp = current_name + ".in_progress"
                if not os.path.exists(current_name_temp):
                    file_exists = False	 # If such a file is not found, i is
                                         # the current decoy #.
                    with open(current_name_temp, 'w') as f:
                        f.write("This decoy is in progress.")
                    self.current_name = current_name
            i += 1
        self.current_num = i - 1
        if file_exists:
            self.job_complete = True

    def output_decoy(self, pose):
        current_name_temp = self.current_name + ".in_progress"
        if not os.path.exists(current_name_temp):
            return

        dump_pdb(pose, self.current_name)  # Outputs pdb file.
        os.remove(current_name_temp)

        score_tag = ".fasc"
        if not pose.is_fullatom():
            score_tag = ".sc"

        scorefile = self.pdb_name + score_tag
        output_scorefile(pose, self.pdb_name, self.current_name, scorefile, self.scorefxn, \
                     self.nstruct, self.native_pose, self.additional_decoy_info)

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


if config['core']:
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
            sts = rosetta.core.vector1_ScoreType()
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
                    self.scoreType = max(s.methods.keys() or [s.first - 1]) + 1
                    if self.scoreType > s.last:
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
