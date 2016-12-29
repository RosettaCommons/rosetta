#!usr/bin/env python

from __future__ import print_function

################################################################################
# A GENERAL EXPLANATION

"""
loop_modeling.py

This script models a loop between residues 77-85 of the PDB file
"test_in.pdb". Loop modeling uses a combination of fragment insertion and
CCD loop closure in a low resolution (centroid) simulated annealing Monte Carlo
protocol to generate a closed loop. This loop is refined using
high resolution (fullatom) small/shear moves to produce low-energy
loop conformations. This script requires the input PDB file,
"test_in.pdb", and the corresponding 3-residue fragment file, "test3_fragments".

Instructions:

1) ensure that your PDB file and fragment file are in you current directory
2) run the script:
    from commandline                        >python D080_Loop_modeling.py

    from within python/ipython              [1]: run D080_Loop_modeling.py

Author: Evan H. Baugh
    based on an original script by Sid Chaudhury
    revised and motivated by Robert Schleif

Updated by Boon Uranukul, 6/9/12
Simplified special constant seed initialization ~ Labonte

References:
    A. A. Cantescu & R. L. Dunbrack, "Cyclic coordinate descent: A robotics
        alforithm for protein loop closure," Protein Sci. 12, 963-972 (2003).
    C. Wang, P. Bradley & D. Baker, "Protein-protein docking with backbone
        flexibility," J. Mol. Biol. 373, 503-519 (2007).

"""

# WARNING
"""
This sample protocol displays a general way to setup a simulated
annealing process in PyRosetta. The process is explicitly written for the
low resolution mode while the high resolution refinement is handled by a
mysterious LoopMover_Refine_CCD. However, the low resolution protocol is based
on the actions (and code) executed by LoopMover_Refine_CCD. As such,
there are two sets of very similar options throughout this sample script but
they are for DIFFERENT PARTS OF THE PROTOCOL and when building your own
protocols be aware that your low resolution and high resolution steps
DO NOT HAVE TO PERFORM SIMILAR ACTIONS. You are free to mix-and-match
any steps you wish.

"""

################################################################################
# THE BASIC PROTOCOL, sample_single_loop_modeling

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing the sample_single_loop_modeling method)

The method sample_single_loop_modeling:
1.  creates a pose from the desired PDB file
2.  creates a copy of the pose (fullatom) for reference
3.  creates a Loop object defining the loop region
4.  modifies the pose FoldTree using the Loop
5.  sets the cut-point residues as cut-point variants
6.  creates a MoveMap object with all chi torsions and
        the loop region backbone torsions free
7.  sets up a ClassicFragmentMover for inserting fragments backbone
        torsions into the loop region
8.  creates a low resolution (centroid) CCD loop closure Mover
        (for closing centroid loops)
9.  creates low and high resolution ScoreFunctions
10. sets up a PackRotamersMover for sidechain packing
11. sets up a high resolution CCD loop closure Mover
        (for fullatom loop optimization)
12. creates Movers for switching between fullatom and centroid and for
        recovering the original sidechain conformations of the fullatom pose
13. creates a copy of the pose (centroid) for reference
14. creates the geometric "temperature" decrement for simulated annealing
15. creates a PyMOLMover for exporting structures to PyMOL
16. creates a (Py)JobDistributor for managing multiple trajectories
17. performs the loop modeling protocol, for each trajectory:
        a.  reset necessary variables for the new trajectory
                -reload the starting pose (centroid)
                -change the pose's PDBInfo.name, for exporting to PyMOL
                -reset the starting "temperature" (to init_temp)
                -create a MonteCarlo object for this trajectory
        b.  "randomize" the structure by:
                -setting loop phi=-180 and psi=180
                -inserting fragments
        c.  perform low resolution (centroid) modeling in
                rounds (outer_cycles) by:
                -recovering the best structure (lowest scoring)
                -performing several simulated annealing steps (inner_cycles) by:
                    >decreasing the "temperature" by one decrement
                    >inserting a fragment into the loop region
                    >closing the loop region using CCD
                    >assessing the pose using the MonteCarlo object
        d.  convert the best structure (lowest scoring) into fullatom
                -recover the best (centroid) structure (lowest scoring)
                -switch the ResidueTypeSet to fullatom (from centroid)
                -recover the original sidechain conformations
                -perform side chain packing
        e.  perform high resolution loop optimization
        f.  output the decoy structure
                -to a PDB file using the PyJobDistributor
                -to PyMOL using the PyMOLMover

"""

import optparse    # for option sorting

from math import pow    # for decrementing kT during simulated annealing
from rosetta.protocols.loops.loop_closure.ccd import *
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta import *
from pyrosetta import *
init(extra_options = "-constant_seed")
# normally, init() works fine
# for this sample script, we want to ease comparison by making sure all random
#    variables generated by Rosetta in this instance of PyRosetta start from a
#    constant seed
# here we provide the additional argument "-constant_seed" which sets all the
#    random variables generated by Rosetta from a constant seed (google random
#    seed for more information)
# some options can be set after initialization, please see PyRosetta.org FAQs
#    for more information
# WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

import os; os.chdir('.test.output')


#########
# Methods

def sample_single_loop_modeling(pdb_filename,
        loop_begin, loop_end, loop_cutpoint,
        frag_filename, frag_length,
        outer_cycles_low = 2, inner_cycles_low = 5,
        init_temp_low = 2.0, final_temp_low = 0.8,
        outer_cycles_high = 5, inner_cycles_high = 10,
        init_temp_high = 2.2, final_temp_high = 0.6,
        jobs = 1, job_output = 'loop_output'):
    """
    Performs simple single loop construction on the input  <pdb_filename>
        with loop from  <loop_begin>  to  <loop_end>  with a
        cutpoint at  <loop_cutpoint>  using fragments of length  <frag_length>
        in the file  <frag_filename>.  <jobs>  trajectories are performed,
        each using a low resolution (centroid) simulated annealing with
        <outer_cycles>  rounds and  <inner_cycles>  steps per round decrementing
        "temperature" from  <init_temp>  to  <final_temp>  geometrically.
        Output structures are named  <job_output>_(job#).pdb.

    """
    # 1. create a pose from the desired PDB file
    p = Pose()
    pose_from_file(p , pdb_filename)

    # 2. create a reference copy of the pose in fullatom
    starting_p = Pose()
    starting_p.assign(p)

    #### if you are constructing multiple loops simultaneously, changes will
    ####    occur in most of the steps below

    # 3. create the Loop object
    #    (note: Loop objects merely specify residues, they contain no
    #         conformation data)
    my_loop = protocols.loops.Loop(loop_begin, loop_end, loop_cutpoint)
    #### if using multiple loops, add additional Loop objects
    # 4. use the Loop to set the pose FoldTree
    protocols.loops.set_single_loop_fold_tree(p, my_loop)
    #### alternate FoldTree setup, if you uncomment the lines below,
    ####    comment-out the set_single_loop_foldtree line above (line 189)
    #### -create an empty FoldTree
    #ft = FoldTree()
    #### -make it a single edge the length of pose
    #ft.simple_tree(p.total_residue())
    #### -insert a jump corresponding to the single loop region
    #ft.add_jump(loop_begin - 2, loop_end + 2, loop_cutpoint)
    #### -give the pose this FoldTree (set it to this object), this will
    ####     erase any previous FoldTree held by the pose
    #p.fold_tree(ft)
    #### there is also a fold_tree_from_loops method in exposed which sets up
    ####    a FoldTree but it is different from set_single_loop_foldtree in
    ####    that is creates jumps +/- 1 residue from their corresponding loop
    ####    endpoints and requires a third argument, the FoldTree to setup

    # 5. sets the cut-point residues as cut-point variants
    protocols.loops.add_single_cutpoint_variant(p, my_loop)

    # 6. create the MoveMap, allow the loop region backbone and
    #    all chi torsions to be free
    movemap = MoveMap()
    movemap.set_bb_true_range(loop_begin, loop_end)
    movemap.set_chi(True)    # sets all chi torsions free

    # 7. setup the fragment Mover
    # this "try--except" is used to catch improper fragment files
    try:
        fragset = core.fragment.ConstantLengthFragSet(frag_length, frag_filename)
        #### the ConstantLengthFragSet is overloaded, this same
        ####    ConstantLengthFragSet can be obtained with different syntax
        # to obtain custom fragments, see Generating Fragment Files below
    except:
        raise IOError('Make sure frag_length matches the fragments in\n\
            frag_file and that frag_file is valid')
    fragment_mover = protocols.simple_moves.ClassicFragmentMover(fragset, movemap)

    # 8. create a Mover for loop modeling using CCD (low resolution)
    ccd_closure = protocols.loops.loop_closure.ccd.CCDLoopClosureMover(my_loop, movemap)

    # 9. create ScoreFunctions
    # for centroid, use the default centroid ScoreFunction with chainbreak on
    scorefxn_low = create_score_function('cen_std')
    # the chainbreak ScoreType exists to penalize broken bonds
    # try creating a broken pose in the interpreter and use a ScoreFunction
    #    with a chainbreak score to investigate its impact, the score is 0.0
    #    except when a bond is broken
    # this penalizes failures caused by CCD failing to close the loop
    scorefxn_low.set_weight(core.scoring.chainbreak, 1)
    # for fullatom, used for packing and scoring final output
    scorefxn_high = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')

    # 10. setup sidechain packing Mover
    task_pack = core.pack.task.TaskFactory.create_packer_task(starting_p)
    task_pack.restrict_to_repacking()    # prevents design, packing only
    task_pack.or_include_current(True)    # considers original sidechains
    pack = protocols.simple_moves.PackRotamersMover(scorefxn_high, task_pack)

    # 11. setup the high resolution refinement
    # by creating a Loops object,
    #    (note: Loops is basically a list of Loop objects),
    sample_loops = protocols.loops.Loops()
    # giving it the loop to remodel,
    sample_loops.add_loop(my_loop)
    # and creating a fullatom CCD Mover (high resolution)
    # this Mover is somewhat abnormal since it handles everything itself, it:
    #    -creates its own MoveMap for the loop regions
    #    -creates its own ScoreFunction (default to get_fa_scorefxn())
    #    -creates its own FoldTree for the pose based on the loops
    #    -creates its own MonteCarlo object for monitoring the pose
    #    -performs "simulated annealing" with 3 outer cycles and 90 inner
    #        cycles, very similar to the protocol outlined ere
    #    -creates its own backbone Movers (SmallMover, ShearMover)
    #    -creates its own PackRotamersMover, it does NOT restrict repacking
    #        to the loop regions and can alter all sidechain conformations
    loop_refine = LoopMover_Refine_CCD(sample_loops)
    # some of these parameters or objects can be set but the protocol
    #    executed by this Mover is effectively untouchable
    #loop_refine.set_score_function(scorefxn_high)    # in beta v2 and above
    loop_refine.temp_initial(init_temp_high)
    loop_refine.temp_final(init_temp_high)
    loop_refine.outer_cycles(outer_cycles_high)
    loop_refine.max_inner_cycles(inner_cycles_high)

    # 12. create centroid <--> fullatom conversion Movers
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    # and a Mover to recover sidechain conformations
    #    when a protocol samples backbone torsion space in centroid,
    #    the sidechain conformations are neglected, when it is transferred
    #    to fullatom, we typically set the sidechain conformations to their
    #    "original" values and perform sidechain packing,
    #    a ReturnSidechainMover saves a pose's sidechains (in this case
    #    staring_pose) and when applied, inserts these conformations
    #    into the input pose
    recover_sidechains = protocols.simple_moves.ReturnSidechainMover(starting_p)

    # 13. create a reference copy of the pose in centroid
    # the first stage of each trajectory is in centroid
    #    so a centroid reference is needed and the pose must start in centroid
    to_centroid.apply(p)
    starting_p_centroid = Pose()
    starting_p_centroid.assign(p)

    # 14. create the geometric "temperature" increment for simulated annealing
    gamma = pow((final_temp_low/init_temp_low),
        (1.0/(outer_cycles_low*inner_cycles_low)))

    # 15. create a PyMOLMover for exporting structures to PyMOL
    pymov = PyMOLMover()
    # uncomment the line below to load structures into successive states
    #pymov.keep_history(True)
    scorefxn_high(starting_p)    # for exporting the scores
    pymov.apply(starting_p)
    pymov.send_energy(starting_p)

    # 16. create a (Py)JobDistributor
    # a PyJobDistributor uses the job_output argument to name all output files
    #    and performs the specified number (int) of jobs
    # a ScoreFunction is required since the PyJobDistributor output .fasc file
    #    contains scoring information about each output PDB
    jd = PyJobDistributor(job_output, jobs, scorefxn_high)
    jd.native_pose = starting_p

    # 17. perform the loop modeling protocol
    counter = 0    # for exporting to PyMOL
    while not jd.job_complete:
        # a. set necessary variables for the new trajectory
        # -reload the starting pose (centroid)
        p.assign(starting_p_centroid)
        # -change the pose's PDBInfo.name, for exporting to PyMOL
        counter += 1
        p.pdb_info().name(job_output + '_' + str(counter) + '_cen')
        # -reset the starting "temperature" (to init_temp)
        kT = init_temp_low
        # -create a MonteCarlo object for this trajectory
        #    a MonteCarlo object assesses pass/fail by the Metropolis Criteria
        #    and also records information on the lowest scoring pose
        mc = MonteCarlo(p, scorefxn_low, kT)

        # b. "randomize" the loop
        #### this section may change if you intend to use multiple loops or
        ####    alter the sampling method to "randomize" the loop
        # -by breaking it open,
        for i in range(loop_begin , loop_end + 1):
            p.set_phi(i , -180)
            p.set_psi(i , 180)
        pymov.apply(p)
        # -and then inserting fragments
        #    the number of insertions performed is somewhat arbitrary
        for i in range(loop_begin, loop_end + 1):
            fragment_mover.apply(p)
        pymov.apply(p)
        ####

        # low resolution loop modeling:
        # c. simulated annealing incrementing kT geometrically
        #    from init_temp to final_temp
        #### this section may change if you intend to use multiple loops or
        ####    alter the sampling method for low resolution modeling
        for i in range( 1, outer_cycles_low + 1):
            # -start with the lowest scoring pose
            mc.recover_low(p)    # loads mc's lowest scoring pose into p
            # -take several steps of in the simulated annealing by
            for j in range(1, inner_cycles_low + 1):
                # >increasing the "temperature"
                kT = kT * gamma
                mc.set_temperature(kT)
                # >inserting a fragment,
                fragment_mover.apply(p)
                pymov.apply(p)
                # >performing CCD,
                ccd_closure.apply(p)
                pymov.apply(p)
                # >and assessing the Metropolis Criteria
                mc.boltzmann(p)
        ####

        # the LoopMover_Refine_CCD makes A LOT of moves, DO NOT expect to
        #    see useful results if you use the PyMOLMover keep_history option, the large
        #    number of intermediates will slow processing to a halt

        # d. convert the best structure (lowest scoring) into fullatom by:
        # -recovering the best (centroid) structure (lowest scoring),
        mc.recover_low(p)    # loads mc's lowest scoring pose into p
        # -switching the ResidueTypeSet to fullatom (from centroid),
        to_fullatom.apply(p)
        # -recovering the original sidechain conformations,
        recover_sidechains.apply(p)
        # -and packing the result (since the backbone conformation has changed)
        pack.apply(p)
        pymov.apply(p)
        p.pdb_info().name(job_output + '_' + str( counter ) + '_fa')

        # high-resolution refinement:
        #### this section may change if you intend to use multiple loops or
        ####    alter the sampling method for high resolution refinement
        # e. apply the LoopMover_Refine_CCD
        loop_refine.apply(p)

        # f. output the decoy (pose result from this trajectory)
        #    include the loop RMSD (Lrsmd)
        # -output a PDB file using the PyJobDistributor
        lrms = protocols.loops.loop_rmsd(p, starting_p, sample_loops, True)
        jd.additional_decoy_info = ' Lrmsd: ' + str(lrms)
        jd.output_decoy(p)
        # -export the structure to PyMOL
        pymov.apply(p)
        pymov.send_energy(p)

    # this step is not absolutely necessary but it is good practice to remove
    #    pose objects attached to the observer before finishing

################################################################################
# INTERPRETING RESULTS

"""
The (Py)JobDistributor will output the lowest scoring pose for each trajectory
(as a PDB file), recording the loop RMSD and score in
<job_output>.fasc. Generally, the decoy generated with the lowest score
contains the best prediction for the loop conformation. If multiple loops
are to be modeled with individual calls to sample_single_loop_modeling,
these loop conformations can be combined to yield the best overall prediction
(assuming each loop does not significantly affect the score of the other loops).
The CCD method CAN produce physically unrealistic results (rarely). As such,
individual inspection of the conformation (such viewing in PyMOL) should
accompany the interpretation of results.

The PyMOLMover sends intermediate and final structures for each trajectory to
PyMOL. The original input structure and the final fullatom structure
(named <job_output>_(job#)_fa) for each trajectory (colored by per-residue
score) are displayed for comparison of the loop conformation. For each
trajectory, the "randomized" centroid structure and outer cycle
(outer_cycles_low) proposed structure(s) are exported into successive states
(named <job_output>_(job#)_cen). Cycle through these states to observe the
simulated annealing protocol output and the loop closure conformation. For a
short number of cycles (such as the default) the loop closure will most likely
succeed at finding a low scoring conformation early causing successive moves
to increase the score and be rejected (this is confounded by the decreasing
acceptance threshold, kT). Thus, for low cycles, the low-resolution steps may
only display the randomized and closed structures (since MonteCarlo.boltzmann
will likely discard the proposed structure if it increases the score).
"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls sample_single_loop_modeling with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb" with reduced
#    cycles/jobs to provide results quickly
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '../test/data/test_in.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel')
# the loop options
parser.add_option('--loop_begin', dest = 'loop_begin',
    default = '15',    # specific to each inquiry, in this case test_in.pdb
    help = 'the starting residue of the loop region to remodel' )
parser.add_option('--loop_end', dest = 'loop_end',
    default = '19',    # specific to each inquiry, in this case test_in.pdb
    help = 'the last residue of the loop region to remodel')
parser.add_option('--loop_cutpoint' , dest = 'loop_cutpoint',
    default = '',    # specific to each inquiry, in this case test_in.pdb
    help = 'the cutpoint residue for the loop region')
# the fragment file options
parser.add_option('--frag_filename', dest = 'frag_filename',
    default = '../test/data/test3_fragments',    # specific to each PDB (test_in.pdb here)
    help = 'the file containing fragments corresponding to the PDB')
parser.add_option('--frag_length', dest = 'frag_length',
    default = '3',    # must match the frag_filename
    help = 'the length of fragments contained in the frag_file')
# low resolution options
parser.add_option('--outer_cycles_low', dest = 'outer_cycles_low',
    default = '2',    # defaults to low value for speed
    help = 'the number of rounds of simulated annealing to perform\
        for each trajectory, low resolution')
parser.add_option('--inner_cycles_low', dest = 'inner_cycles_low',
    default = '5',    # defaults to low value for speed
    help = 'the number of steps in a single simulated annealing,\
        low resolution')
parser.add_option('--init_temp_low',dest ='init_temp_low',
    default = '2.0',    # commonly used higher "temperature"
    help = 'the initial \"temperature\" of the simulated annealing,\
        low resolution')
parser.add_option('--final_temp_low', dest = 'final_temp_low',
    default = '0.8',    # commonly used lower "temperature"
    help = 'the final \"temperature\" of the simulated annealing,\
        low resolution')
# high resolution options
parser.add_option('--outer_cycles_high', dest = 'outer_cycles_high',
    default = '5',    # defaults to low value for speed
    help = 'the number of rounds of simulated annealing to perform\
        for each trajectory, high resolution')
parser.add_option('--inner_cycles_high', dest = 'inner_cycles_high',
    default = '10',    # defaults to low value for speed
    help = 'the number of steps in a single simulated annealing,\
        high resolution')
parser.add_option('--init_temp_high',dest ='init_temp_high',
    default = '2.2',    # commonly used higher "temperature"
    help = 'the initial \"temperature\" of the simulated annealing,\
        high resolution')
parser.add_option('--final_temp_high', dest = 'final_temp_high',
    default = '0.6',    # commonly used lower "temperature"
    help = 'the final \"temperature\" of the simulated annealing,\
        high resolution')
# the JobDistributor options
parser.add_option('--jobs', dest='jobs',
    default = '1',    # default to single trajectory for speed
    help = 'the number of jobs (trajectories) to perform')
parser.add_option('--job_output', dest = 'job_output',
    default = 'loop_output',    # if a specific output name is desired
    help = 'the name preceding all output, output PDB files and .fasc')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename

# loop options
loop_begin = int(options.loop_begin)
loop_end = int(options.loop_end)
# default the loop cutpoint to the average of loop_begin and loop_end
if options.loop_cutpoint:
    loop_cutpoint = int(options.loop_cutpoint)
else:
    loop_cutpoint = (loop_begin + loop_end) // 2
# fragment  options
frag_filename = options.frag_filename
frag_length = int(options.frag_length)
# low resolution modeling options (simulated annealing)
outer_cycles_low = int(options.outer_cycles_low)
inner_cycles_low = int(options.inner_cycles_low)
init_temp_low = float(options.init_temp_low)
final_temp_low = float(options.final_temp_low)
# high resolution modeling options (a different simulated annealing)
outer_cycles_high = int(options.outer_cycles_high)
inner_cycles_high = int(options.inner_cycles_high)
init_temp_high = float(options.init_temp_high)
final_temp_high = float(options.final_temp_high)
# JobDistributor options
jobs = int(options.jobs)
job_output = options.job_output

# perform the primary method of this script
sample_single_loop_modeling(pdb_filename,
    loop_begin, loop_end, loop_cutpoint,
    frag_filename, frag_length ,
    outer_cycles_low, inner_cycles_low,
    init_temp_low, final_temp_low ,
    outer_cycles_high, inner_cycles_high,
    init_temp_high, final_temp_high,
    jobs, job_output)


################################################################################
# ALTERNATE SCENARIOS

############################
# Generating Fragment Files:
"""
You MUST create a fragment file for each new PDB protein sequence that you
wish to model. New fragment files are obtained from the Robetta Server:

1) obtain the protein sequence (or FASTA file) of the protein contained
    within your PDB of interest
2) Go to:
    http://robetta.bakerlab.org/fragmentsubmit.jsp
    and submit the protein sequence (ignore the Optional section unless
    you are familiar with Robetta)

"""

#################
# A Real Example:
"""
All of the default variables and parameters used above are specific to
the example with "test_in.pdb", which is supposed to be simple,
straightforward, and speedy. Here is a more practical example:

A loop region of Triosephosphate Isomerase is theorized to stabilize an
enediol intermediate formed in the enzyme's active site. Suppose you are
interested in the conformations of this loop (pose numbered residues 159-169)
and decide to model the loop in PyRosetta.

1. Download a copy of RCSB PDB file 3S6D (remove waters and any other HETATM)
2. Make a fragment file of 3-mers using the "Generate Fragment Files"
        (instructions above)
3. Make a directory containing:
        -the PDB file for 3S6D (cleaned of HETATMs and waters)
            lets name it "3S6D.clean.pdb" here
        -the 3-mer fragment file for 3S6D
            lets name it "3S6D.frag3" here
        -this sample script (technically not required, but otherwise the
            commands in 4. would change since loop_modeling.py would't be here)
4. Run the script from the commandline with appropriate arguments:

>python loop_modeling.py --pdb_filename 3S6D.clean.pdb --loop_begin 159 --loop_end 169 --frag_filename 3S6D.frag3 --frag_length 3 --jobs 40 --job_output 3S6D_159_169_loop_output --outer_cycles_low 5 --inner_cycles_low 10 --init_temp_low 2. --final_temp_low .8 --outer_cycles_high 3 --inner_cycles_high 90 --init_temp_high 2. --final_temp_high .8

        -The option --loop_cutpoint was left blank causing it to default to the
            "midpoint" of the loop (159+169)/2 = 164
        -The frag_length MUST match frag_filename (3 in this example)
        -40 trajectories is low, sampling loop conformations is difficult,
            typically hundreds (800-1000) trajectories are attempted
        -There are no common values for outer_cycles_low and inner_cycles_low,
            more steps indicates greater sampling with fragments (in this case)
        -The LoopMover_Refine_CCD defaults to outer_cycles_high = 3 and
            inner_cycles_high = 90, once again, increased sampling searches
            a greater space
        -The simulated annealing "temperatures" of kT=2.0 and kT=.8
            provide a useful range for the Metropolis Criteria for both
            low resolution and high resolution

5. Wait for output, this will take a while (performing 40 trajectories
        of loop modeling involving 5*10 total steps of simulated annealing
        per trajectory)
6. Analyze the results (see INTERPRETING RESULTS above)

Note: this is NOT intended to be used for realistic sampling of
loop conformations, it merely provides a "skeleton" for loop modeling
code in PyRosetta. It may be useful for preliminary investigation but
the best protocols are somewhat protein-specific, there is no current
general loop modeling method.

"""

######################################
# Changing Loop Conformation Sampling:
"""
The protocol here uses two methods for sampling loop conformations,
specifically:
    -"randomization" using fragment insertion
    -loop closure using CCD

Fragment insertion significantly reduces meaningless search space by
only attempting backbone conformations resembling those of the PDB (Database).
More conservative moves can be used here however increasing the search space
drastically increases the computation time required.

Cyclic Coordinate Descent is an effective algorithm for loop closure and
almost always succeeds at closing a loop, however it is not directly
constrained by any realistic scoring. As such, some closed loops score poorly
or represent very unlikely conformations. Several alternatives exist for loop
modeling and the usefulness of results for each algorithm may be
protein specific.

The fullatom loop optimization algorithm, LoopMover_Refine_CCD uses its own
ScoreFunction and is useful for loop optimization. This refinement is
analogous to relaxation (high resolution refinement) in other Rosetta
applications and is a combination of small torsion moves and sidechain
packing. Although expensive, this step (or another high resolution step)
improves results significantly. Try removing this refinement and
observe the difference in loop predictions.

Since low resolution modeling and high resolution refinement sample
stochastically, sufficient sampling should expose the lowest scoring pose
(if you are looking for it). As such, both steps are effectively "seeded"
by the "randomization" step. Searching with fragments should thus provide
an adequate sample space.

Please try alternate sampling methods to better understand how these
algorithms perform and find what moves best suite your problem.

"""

#####################################
# Changing Loop Conformation Scoring:
"""
The protocol here uses out-of-the-box scoring functions for several
applications. The chainbreak ScoreType exists to penalize improper
loops and should be used in loop modeling. The weights used here are
in NO way general or best-suited for any single application. Furthermore
this method is searching using the Metropolis Criteria applied by
the MonteCarlo object.

Please try alternate scoring functions or unique selection methods to better
understand which scoring terms contribute to performance and find what
scoring best suites your problem.

"""

##########################
# Modeling Multiple Loops:
"""
If a protein contains multiple loop regions of interest, it is best to remodel
each loop separately and combine the best results.
Modeling the loops simultaneously either sacrifices searching efficiency or
loops like multiple calls to sample_single_loop_modeling.

If you desire to modify this code to account for multiple loops, be aware
if your Movers are sampling multiple loops simultaneously
and if the MonteCarlo object is thus rejecting potentially useful
loop conformations due to low scores of other loops.

The sample script can be easily modified to include multiple loops.
Steps 3-7, 10, and of course 14 (defined above) should be modified where
appropriate for multiple loops.
To properly model each loop separately, you will need instances of most
Movers for each loop. The setup_single_loop_fold_tree will also be improper,
instead create your own FoldTree with jumps corresponding to each loop region.
It is common to define jumps corresponding to loops with jump_points
2 residues below the lowest numbered loop residue and
2 residue above the highest numbered residue
(remember, this is pose numbering).
for example:
        Loop1 = Loop(77, 85, 81)
        Loop2 = Loop(10, 30, 20)
        Loop3 = Loop(50, 60, 55)
        ft = FoldTree()
        ft.simple_tree(pose.total_residue())
        ft.add_jump(75, 87, 81)
        ft.add_jump(8, 32, 20)
        ft.add_jump(48, 62, 55)
        pose.fold_tree(ft)

"""
