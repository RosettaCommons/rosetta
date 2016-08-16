#!/usr/bin/env python

################################################################################
# A GENERAL EXPLANATION

"""
packer_task.py

This sample script explains the PackerTask object frequently changed in
PyRosetta protocols. It demonstrates how to set up a proper PackerTask for
basic applications and applies this for sidechain repacking and design. Two
additional methods are provided for generating resfiles for a PackerTask.

Sidechain packing (or just packing) refers to the concept of protein sidechain
conformations attaining a physically stable structure dependent upon all
residues involved (such as the hydrophobic core). Van der Waals forces provide
the constraints defining what sets of sidechain conformations are or are not
sufficiently "packed" (hydrophobic separation). Since Rosetta sidechain
conformation selection uses a fixed backbone, at least in any individual step,
the process of selecting low scoring sidechain conformations (those obeying
Van der Waals constraints) is termed "packing" or "repacking".


Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D050_Packer_task.py

    from within python/ipython              [1]: run D050_Packer_task.py

Author: Evan H. Baugh
    revised and motivated by Robert Schleif

Last updated by Boon Uranukul, 6/9/12

References:
    R. L. Dunbrack Jr. and F. Cohen, "Bayesian statistical analysis of protein
        side-chain rotamer preferences," Prot. Science 6 1661-1681 (1997).
"""

################################################################################
# THE BASIC PROTOCOL, packer_task

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing all methods below)

The method packer_task:
1.  sets up a PackerTask for performing sidechain repacking
2.  performs packing
        -display the change in score
        -write the packed pose to a PDB file (packed.pdb)
3.  sets up a PackerTask for design on a set of residues
4.  performs design
        -display the change in score
        -display the change in sequence
        -write the designed pose to a PDB file (designed.pdb)

"""

import optparse    # for option sorting
import rosetta.core.pack.task    # for using resfiles

from rosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


def packer_task(pose, PDB_out = False):
    """
    Demonstrates the syntax necessary for basic usage of the PackerTask object
		performs demonstrative sidechain packing and selected design
		using  <pose>  and writes structures to PDB files if  <PDB_out>
		is True

    """
    # create a copy of the pose
    test_pose = Pose()
    test_pose.assign(pose)

    # this object is contained in PyRosetta v2.0 and above
    pymover = PyMOL_Mover()

    # create a standard ScoreFunction
    scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')

    ############
    # PackerTask
    # a PackerTask encodes preferences and options for sidechain packing, an
    #    effective Rosetta methodology for changing sidechain conformations, and
    #    design (mutation)
    # a PackerTask stores information on a per-residue basis
    # each residue may be packed or designed
    # PackerTasks are handled slightly differently in PyRosetta
    ####pose_packer = PackerTask()    # this line will not work properly
    pose_packer = standard_packer_task(test_pose)
    # the pose argument tells the PackerTask how large it should be

    # sidechain packing "optimizes" a pose's sidechain conformations by cycling
    #    through (Dunbrack) rotamers (sets of chi angles) at a specific residue
    #    and selecting the rotamer which achieves the lowest score,
    #    enumerating all possibilities for all sidechains simultaneously is
    #    impractically expensive so the residues to be packed are individually
    #    optimized in a "random" order
    # packing options include:
    #    -"freezing" the residue, preventing it from changing conformation
    #    -including the original sidechain conformation when determining the
    #        lowest scoring conformation
    pose_packer.restrict_to_repacking()    # turns off design
    pose_packer.or_include_current(True)    # considers original conformation
    print pose_packer

    # packing and design can be performed by a PackRotamersMover, it requires
    #    a ScoreFunction, for optimizing the sidechains and a PackerTask,
    #    setting the packing and design options
    packmover = PackRotamersMover(scorefxn, pose_packer)

    scorefxn(pose)    # to prevent verbose output on the next line
    print '\nPre packing score:', scorefxn(test_pose)
    test_pose.pdb_info().name('original')    # for PyMOL_Mover
    pymover.apply(test_pose)

    packmover.apply(test_pose)
    print 'Post packing score:',  scorefxn(test_pose)
    test_pose.pdb_info().name('packed')    # for PyMOL_Mover
    pymover.apply(test_pose)
    if PDB_out:
        test_pose.dump_pdb('packed.pdb')

    # since the PackerTask specifies how the sidechains change, it has been
    #    extended to include sidechain constitutional changes allowing
    #    protein design, this method of design is very similar to sidechain
    #    packing; all rotamers of the possible mutants at a single residue
    #    are considered and the lowest scoring conformation is selected
    # design options include:
    #    -allow all amino acids
    #    -allow all amino acids except cysteine
    #    -allow specific amino acids
    #    -prevent specific amino acids
    #    -allow polar amino acids only
    #    -prevent polar amino acids
    #    -allow only the native amino acid
    # the myriad of packing and design options can be set manually or, more
    #    commonly, using a specific file format known as a resfile
    #    resfile syntax is explained at:
    #    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_resfiles.html
    # manually setting deign options is tedious, the methods below are handy
    #    for creating resfiles
    # mutate the "middle" residues
    center = test_pose.total_residue() / 2
    specific_design = {}
    for i in range(center - 2, center + 3):
        specific_design[i] = 'ALLAA'
    # write a resfile to perform these mutations
    generate_resfile_from_pose(test_pose, 'sample_resfile', False,
        specific = specific_design)

    # setup the design PackerTask, use the generated resfile
    pose_design = standard_packer_task(test_pose)
    rosetta.core.pack.task.parse_resfile(test_pose, pose_design, 'sample_resfile' )
    print pose_design

    # prepare a new structure
    test_pose.assign(pose)

    # perform design
    designmover = PackRotamersMover(scorefxn, pose_design)
    print '\nDesign with all proteogenic amino acids at (pose numbered)\
        residues', center - 2, 'to', center + 2
    print 'Pre-design score:', scorefxn(test_pose)
    print 'Pre-design sequence: ...' + \
        test_pose.sequence()[center - 5:center + 4] + '...'


    designmover.apply(test_pose)    # perform design
    print '\nPost-design score:', scorefxn(test_pose)
    print 'Post-design sequence: ...' + \
        test_pose.sequence()[center - 5:center + 4] + '...'
    test_pose.pdb_info().name( 'designed' )    # for PyMOL_Mover
    pymover.apply(test_pose)
    if PDB_out:
        test_pose.dump_pdb('designed.pdb')

################################################################################
# PERIPHERAL METHODS

# These methods are not called by the main protocol but address a
#    relevant topic, objects which work intimately with pose.

# writes a specified resfile for the user, defaults to packing w/ input SC
def generate_resfile_from_pose(pose, resfilename,
        pack = True, design = False, input_sc = True,
        freeze = [], specific = {}):
    """
    Writes a resfile for  <pose>  named  <resfilename>
       <pack> = True allows packing by default
       <design> = True allows design using all amino acids by default
       <input_sc> = True allows usage of the original side chain conformation
       <freeze> is an optional list of (pose) residue numbers to exclude
            (preserve the side chain conformations of these residues)
       <specific> is an optional dictionary with (pose) residue numbers as keys
            and resfile keywords as corresponding values
            (for setting individual residue options, it may be easier to add
            these numbers to freeze and edit the resfile manually)

    example:
        generate_resfile_from_pose(pose,'1YY8.resfile')
    See also:
        Pose
        PackRotamersMover
        TaskFactory

    """
    # determine the header, default settings
    header = ''
    if pack:
        if not design:
            header += 'NATAA\n'
        else:
            header += 'ALLAA\n# ALLAA will NOT work on bridged Cysteines\n'
    else:
        header += 'NATRO\n'
    if input_sc:
        header += 'USE_INPUT_SC\n'
    to_write = header + 'start\n'
    # add  <freeze>  list to  <specific>  dict
    for i in freeze:
        specific[i] = 'NATRO'
    #  <specific>  is a dictionary with keys() as pose resi numbers
    #    and values as resfile keywords (PIKAA
    # use PDBInfo object to write the resfile
    info = pose.pdb_info()
    # pose_from_sequence returns empty PDBInfo, Pose() makes NULL
    if info and info.nres():
        for i in specific.keys():
            num = pose.pdb_info().number(i)
            chain = pose.pdb_info().chain(i)
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    else:
        for i in specific.keys():
            num = i
            chain = ' '
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    f = open(resfilename,'w')
    f.write(to_write)
    f.close()

# this is silly, as with cleanCRYS, later implement a Biopy PDBParser method
def generate_resfile_from_pdb(pdbfilename, resfilename,
        pack = True, design = False, input_sc = True,
        freeze = [], specific = {}):
	"""
    Writes a resfile for the PDB file <pdbfilename>  named  <resfilename>
       <pack> = True allows packing by default
       <design> = True allows design using all amino acids by default
       <input_sc> = True allows usage of the original side chain conformation
       <freeze> is an optional list of (pose) residue numbers to exclude
            (preserve the side chain conformations of these residues)
       <specific> is an optional dictionary with (pose) residue numbers as keys
            and resfile keywords as corresponding values
            (for setting individual residue options, it may be easier to add
            these numbers to freeze and edit the resfile manually)

	example:
	    generate_resfile_from_pdb('1YY8.pdb','1YY8.resfile')
	See also:
	    generate_resfile_from_pose
	    Pose
	    PackRotamersMover
	    TaskFactory

	"""
	p = pose_from_file(pdbfilename)
	generate_resfile_from_pose(p, resfilename, pack,design, input_sc, freeze, specific)

################################################################################
# INTERPRETING RESULTS

"""
The PDB files output demonstrating the PackerTask for packing and design
should be compared to the original PDB file (test_in.pdb for the sample case).
Notice that packed.pdb has the original backbone conformation but many of
the sidechain conformations differ. For test_in.pdb, this packing step provides
the quickest visible change in the score value.
Notice that designed.pdb has the original backbone conformation but the "middle"
residues have different amino acids than the original sequence. If you use a
PDB file other than test_in.pdb, these amino acids may preserve their identity,
it depends on their score evaluation in PyRosetta and it there is an accessible
lower score rotamer.
Since a PackerTask can set options to due packing on some residues and design
on others simultaneously, it is IMPERATIVE that you isolate your desired move.
For example, if the design performed in the sample script above also repacked
several other residues, would the score change be due to repacking or design?
Would the amino acids "mutated" have their lowest scoring conformation? This
confounding possibility is why structures must be made Rosetta friendly before
applying a protocol. Many Rosetta "moves" may lower the score due to trivial
structural changes instead of the protocol of interest.
In some cases, design may change the rotamer, but not the amino acid identity,
indicating that the new rotamer has the lowest score of all available rotamers.
"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls packer_task with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb"
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '../test/data/test_in.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel' )
parser.add_option('--PDB_out', dest = 'PDB_out',
    default = '',    # default to False
    help = 'flag for dumping structures to PDB files' )
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
# create a pose from the desired PDB file
# -create an empty Pose object
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# PDB_out flag
PDB_out = bool(options.PDB_out)

packer_task(pose, PDB_out)
