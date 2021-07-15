# @file: predict_ddG.py
# @brief: Calculated the ddG of mutation dependent on membrane depth and pH
# @notes: Adapted for benchmark and PyRosetta4 from Alford & Koehler Leman et al. 2015
# @author: Rebcca Alford (ralford3@jhu.edu)
# Modified on Apr 6th, 2021 by Rituparna Samanta (rsamant2@jhu.edu)

from pyrosetta import *
from pyrosetta.teaching import *
from string import Template

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

import sys, os
import random

import pyrosetta.rosetta.protocols.membrane
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import PDBInfo
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pack.task import TaskFactory

# Make a dictionary for classifying residue types (yes I know, global vars are bad practice)
classification = {'A': "nonpolar", 'C' : "special", 'D' : "n-charged", 'E' : "n-charged", \
                         'F' : "aromatic", 'G' : "special", 'H' : "p-charged", 'I' : "nonpolar", \
                         'K' : "p-charged", 'L' : "nonpolar", 'M' : "nonpolar", 'N' : "polar", \
                         'P' : "special", 'Q' : "nonpolar", 'R' : "p-charged", 'S' : "nonpolar", \
                         'T' : "polar", 'V' : "nonpolar", 'W' : "aromatic", 'Y' : "aromatic" }

# @brief Replace the residue at <resid> in <pose> with <new_res> and allows
# repacking within a given <pack_radius>
def mutate_residue( pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn ):

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

    # Create a packer task (standard)
    task = TaskFactory.create_packer_task( test_pose )

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range( 1 , 21 ):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task( mutant_position
        ).restrict_absent_canonical_aas( aa_bool )

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1, pose.total_residue() + 1 ):
        dist = center.distance_squared( test_pose.residue( i ).nbr_atom_xyz() );
        # only pack the mutating residue and any within the pack_radius
        if i != mutant_position and dist > pow( float( pack_radius ), 2 ) :
            task.nonconst_residue_task( i ).prevent_repacking()
        elif i !=mutant_position and dist <= pow( float( pack_radius ), 2 ) :
            task.nonconst_residue_task( i ).restrict_to_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover( pack_scorefxn , task )
    packer.apply( test_pose )

    return test_pose

def main( args ):

    # Read options from the commandline
    parser = OptionParser( usage="usage: %prog --mutation_list ddGs.txt --energy_fxn mpframework_fa_2007" )
    parser.set_description(main.__doc__)

    parser.add_option('--pdb', '-p',
        action="store",
        help="Input list of mutations",)

    parser.add_option('--spanfile', '-s',
        action="store",
        help="Input list of mutations",)

    parser.add_option('--mutation_list', '-m',
        action="store",
        help="Input list of mutations",)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Energy function weights for this benchmark",)

    parser.add_option('--repack_radius', '-r',
        action="store",
        help="Repack radius for making mutations")

    parser.add_option('--outdir', '-o',
        action="store",
        help="Output directory for the data file",)

    parser.add_option('--restore', '-t', 
        action="store", 
        help="Restore talaris behavior to pre-ref2015 for refernce benchmarks", )

    parser.add_option( '--implicit_lipids', '-i', 
        action="store", 
        help="Use implicit lipids and default parameters when running this benchamrk", )

    # parse arguments
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Error checking
    if ( not Options.mutation_list ):
        print("Missing rquired option mutations list! Exiting...")
        sys.exit()

    if ( not Options.energy_fxn ):
        print("Missing required energy function weights! Exiting...")
        sys.exit()

    # Initialize Pyrosetta with const options
    option_string = "-run:constant_seed -in:ignore_unrecognized_res"
    if ( Options.implicit_lipids ): 
        option_string = option_string + " -mp:lipids:temperature 20.0 -mp:lipids:composition DLPC -mp:lipids:has_pore false"
        #changed from 37 to 20 as per the paper. 
    init( extra_options=option_string )

    # Read database file including mutations (space delimited)
    # Expected header format: Nat Pos Mut PDb Spanfile pH exp_ddG double_mut
    with open( Options.mutation_list, 'r' ) as f:
        content = f.readlines()
    content = [ x.strip() for x in content ]
    mutation_list = [ x.split(' ') for x in content ]

    # Create an energy function
    print("going to make energy function")
    sfxn = create_score_function( Options.energy_fxn )
    sfxn.set_weight(fa_water_to_bilayer, 1.5)
    #Results in the paper are with respect to the weight being 1.50
    print("done")
    # Set the repack radius from the option system
    repack_radius = 8.0
    ##RS: a repack_radius is also taken as command line. 

    # Setup output file
    outfile = Options.outdir + "/ddG_" + Options.energy_fxn + ".dat"
    f = open( outfile, 'a' )
    f.write( "Nat Pos Mut experimental_ddG predicted_ddG class depth\n" )

    for entry in mutation_list:

        # Sanity check that the PDB file path and spanfile path exist
        if ( not os.path.isfile(Options.pdb) or not os.path.isfile(Options.spanfile) ):
            print(Options.pdb + " " + Options.spanfile)
            print("Path to PDB file or spanfile is invalid!")
            sys.exit()

        # Read PDB from table - note, must contain an absolute path
        pose = pose_from_pdb( Options.pdb )

        # Add membrane to pose
        add_memb = rosetta.protocols.membrane.AddMembraneMover( Options.spanfile )
        add_memb.apply( pose )

        # Setup in a topology based membrane
        init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
        init_mem_pos.apply( pose )

        # Calculate the native state, based on whether or not there is a shift to alanine first
        if ( entry[7] == "Y" ):
            repacked_native = mutate_residue( pose, int( entry[1] ), 'A', repack_radius, sfxn )
        else:
            native_res = pose.residue( int( entry[1] ) ).name1()
            repacked_native = mutate_residue( pose, int( entry[1] ), native_res, repack_radius, sfxn )

        # Calculate the ddG of mutation for the given position
        ddG_of_mutation = compute_ddG( repacked_native, sfxn, int( entry[1] ), entry[2], repack_radius )
        print(ddG_of_mutation)
        print_score_labels_to_file( repacked_native, sfxn, "dummy" )

        # Calculate some additional classifications for the mutation
        depth = pose.conformation().membrane_info().residue_z_position( pose.conformation(), int( entry[1] ) )
        rsd_class = classification[ entry[2] ]

        # Write ddG data to output file
        outstr = Template( "$native $pos $mutant $exp_val $predicted $rsd_class $depth" )
        output = outstr.substitute( native=entry[0], pos=entry[1], mutant=entry[2], exp_val=entry[6], predicted=ddG_of_mutation, rsd_class=rsd_class, depth=round(depth,3))
        f.write( output + "\n" )

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa, repack_radius ):

    # Score the native pose and grab the native AA
    native_score = sfxn( pose )
    native_aa = pose.residue( resnum ).name1()

    # Perform the mutation at residue <resnum> to amino acid <aa> and score
    mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )
    mutant_score = sfxn( mutated_pose )
    
    print_ddG_breakdown( pose, mutated_pose, sfxn, resnum, aa, "breakdown.sc" )

    # Calculate the ddG in place
    ddG = round( mutant_score - native_score, 3 )
    return ddG

###############################################################################
#@brief Print ddG breakdown from the pose
# Extract weighted energies from the native and mutated pose. Calculate the ddG
# of each and print the component-wise ddG vlaues
def print_ddG_breakdown( native_pose, mutated_pose, sfxn, resnum, aa, fn ): 

    # Extract scores
    tmp_native = native_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    tmp_mutant = mutated_pose.energies().total_energies().weighted_string_of( sfxn.weights() )

    # Parse out scores
    array_native = list(filter( None, tmp_native.split(' ') ))
    array_mutant = list(filter( None, tmp_mutant.split(' ') ))

    # Pull out only the scores from these arrays
    native_scores = []
    for i in range( len(array_native) ): 
        if ( i % 2 != 0 ): 
            native_scores.append( float( array_native[i] ) )

    mutant_scores = []
    for i in range( len(array_mutant) ): 
        if ( i % 2 != 0 ): 
            mutant_scores.append( float( array_mutant[i] ) )

    # Make a label for the mutation
    native_res = native_pose.residue( int( resnum ) ).name1()
    mut_label = native_res + str(resnum) + aa

    # Calculate ddG of individual components
    ddGs = []
    ddGs.append( mut_label )
    for i in range( len( mutant_scores ) ): 
        ddG_component = mutant_scores[i] - native_scores[i]
        ddGs.append( round( ddG_component, 3 ) )

    ddGs_str = convert_array_to_str( ddGs ) 
    with open( fn, 'a' ) as f:
        f.write( ddGs_str + "\n" )
    f.close()

###############################################################################
#@brief Get header for ddG breakdown output
# Save the score labels, to be printed at the top of the output breakdown file
def print_score_labels_to_file( native_pose, sfxn, fn ): 

    tmp_native = native_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    array_native = list(filter( None, tmp_native.split(' ') ))
    labels = []
    labels.append( 'mutation ' ) # Append field for mutation label
    for i in range( len(array_native) ): 
        if ( i % 2 == 0 ): 
            labels.append( array_native[i].translate(':') )

    labels_str = convert_array_to_str( labels )
    print(labels_str)
#    with file( fn, 'a' ) as f:
#        f.write( labels_str + "\n" )
#    f.close()

###############################################################################
#@brief Convert an array to a space deliminted string
# Save the score labels, to be printed at the top of the output breakdown file
def convert_array_to_str( array ): 

    linestr = ""
    for elem in array: 
        if ( linestr == "" ): 
            linestr = linestr + str( elem )
        else: 
            linestr = linestr + " " + str( elem )

    return linestr

if __name__ == "__main__" : main(sys.argv)
