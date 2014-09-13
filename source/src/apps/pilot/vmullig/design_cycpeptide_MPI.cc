/*

	design_cycpeptide_multistate.cc
	Created 14 September 2012 by Vikram K. Mulligan
	Modified 21 September 2012 by VKM:
		--Changed scoring function
		--Prevented mutations to cysteine
		--Added RMS weighting
	Modified 26 September 2012 by VKM:
		--Added v_ignoreresidue option
	Modified 28 September 2012 by VKM:
		--There was a problem in that positive state energy in the numerator and denominator of the probability
		function could have different values (due to the stochastic nature of the minimization).  This meant
		that probabilities greater than 1.0 were possible in some cases.  To eliminate this, the positive state
		is now only processed once per iteration (which makes the computation somewhat more efficient.  Note that
		it should NOT be included in the list of negative states loaded with the -in:file:s option OR the flag
		v_removelowrms should be set to true.
	Modified 29 September 2012 by VKM:
		--Fixed a bug in the output of file names.

	This takes a positive backbone state (using the -in:file:native input option) and many negative backbone states
	(using the -in:file:s input option), then uses a Monte Carlo search to try to find a sequence that favours the
	positive state and disfavours the negative states.  The positive state is assumed to be of the form C(X)_nC.  The
	negative states can be both of the form C(X)_nC (from the PDB) or GC(X)_nCG (from Vageli's datasets).  The score
	used in the search is:

	Score = E_pos - k_B*T*ln( sum_i( exp( -E_i/(k_B*T) ) )

	This is equal to -k_B*T*ln( Prob(positive) )

	Modified 12 October 2012 by VKM:
		--The score function is now:
			Score = -k_B*T*ln( (sum_i (exp(-RMS^2/lambda^2)*exp(-E_i/(k_B*T)) / sum_i( RMS^N*exp(-E_i/(k_B*T)) ) )
		--The above now favours ANY states that come close to target.
		--Implementing a short Monte Carlo-based minimization around each state to ensure best minimization.
		--Cleaned up some duplicated code.

	Modified 6 November 2012 by VKM:
		--Parallelizing the algorithm using MPI.
	Modified 21 November 2012 by VKM:
		--Continuing to  parallelize the algorithm using MPI (I was away for a week).
	Modified 26 November 2012 by VKM:
		--Finished parallelizing the algorithm.
	Modified 27 November 2012 by VKM:
		--Improved efficiency of transmission of information between processes somewhat
	Modified 29 November 2012 by VKM:
		--Fixed crash bug with parallelization when FASTA sequence is specified.
	Modified 5 December 2012 by VKM:
		--There's a very hard-to-reproduce deadlock bug somewhere.  After ~660 rounds, the program deadlocked.  On restarting,
		it went about 3000 rounds without deadlocking.  I've switched the MPI_Send command to MPI_Ssend commands (synchronous
		sends) to avoid possible problems from buffering, but I don't know whether I've found the bug or not.
	Modified 15 April 2013 by VKM:
		--Added a final round of Cartesian minimization as an option.
	Considerably rewritten starting 12 June 2013 by VKM:
	DONE:
		--Make support for stripping terminal glycines into an option
		--Change memory handling (store backbone dihedrals rather than full poses; reconstruct poses as needed).
		--Add support for PCA vectors:
			--Read
			--Use in permutations
		--Add support for true cyclic peptides:
			--Option to connect N- and C-termini
			--Option to try cyclic permutations of backbone dihedral vectors
			--Update the use_in_rmsd function (more options)
	Modified 10 July 2013 by VKM:
		--Added option to preserve chirality (use D-amino acids).
		--Added option to include or exclude glycine from possible mutations.
	Modified 16 July 2013 by VKM:
		--Added option to split replicates over processes.
	Modified 5 Aug 2013 by VKM:
		--Added v_savememory option.  If this is on, then only the master proc stores all states and PCA matrices; these are sent as needed to slave procs.
	Modified 12 Aug 2013 by VKM:
		--Added v_scoreonly option.  If true, it scores the input sequence and exits.
	Modified 14 Aug 2013 by VKM:
		--Taking out v_preserve_cys option.  Replacing it with a list of disulfide-bonded cysteines.
	Modified 10 Sept 2013 by VKM:
		--Adding support for beta-amino acids.
	Modified 12 Nov 2013 by VKM:
		--Adding a flag to allow the user to limit the number of PCA vectors stored (to save memory).
	Modified 22 Feb 2014 by VKM:
		--Caught a stupid mistake in the math: DG = -k_b*T*ln(Keq), NOT -k_b*T*ln(prob(folded)).  This did not invalidate any results, but it did mean that during a Monte Carlo search, I had an unduly high probability of accepting a move that increases the energy if I was already very close to a good folder.
	Modified 23 May 2014 by VKM:
		--Adding support for silent file input.
	Modified 3 June 2014 by VKM:
		--Support for silent files now fully implemented.
		--Added support for a user-defined list of additional atoms to use in RMSD
		calculation.
	Modified 29 June 2014 by VKM:
		--Switching MPI_Ssend statements back to MPI_Send for significant performance boost on Blue Gene/Q.
	Modified 22 August 2014 by VKM:
		--Adding some features for Per to allow use to design for rigidity.  We need the following options:
		--v_norepack_positions (positions that will not be repacked during the FastRelax steps.  These are automatically non-designable positions, too.)
		--v_consider_all_replicates (allow all replicates of a given state to be used in evaluation of the funnel shape, rather than just the lowest-energy replicate).
*/

#include <mpi.h>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <devel/init.hh>
//#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
//#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/xyz.functions.hh>
//#include <numeric/xyz.io.hh>
//#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/kinematics/MoveMap.hh>
//#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <numeric/random/random.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <numeric/xyzVector.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/NamedAtomID.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>

//For importing poses:
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

OPT_KEY( Integer, v_trialcount ) //The number of MC moves when searching sequence space
OPT_KEY( Integer, v_relaxrounds ) //The number of rounds of fastrelax to apply to the structures that are rigourously relaxed
OPT_KEY( Real, v_MCtemperature ) //The temperature for the application of the Metropolis criterion when searching sequence space. 1.0 by default.
OPT_KEY( Real, v_kbt ) //The Boltzmann temperature for estimating fractional population of the positive state (i.e. for scoring each sequence). 1.0 by default.
OPT_KEY( Real, v_kbt_minimize) //The Boltzmann temperature for the local Monte Carlo search of backbone conformational space when minimizing each structure. 1.0 by default.
OPT_KEY( Integer, v_MCminimize_replicates ) //The number of MC trajectories to use to minimize each backbone conformation.  Default 5.
OPT_KEY( Integer, v_MCminimize_steps ) //The number of steps in each MC trajectory when minimizing each backbone conformation.  Deafult 5.
OPT_KEY( Real, v_bb_perturbation ) //The mean magnitude of the backbone perturbation vector (vector of changes to backbone dihedrals) when searching for minimum-energy backbone conformation.
OPT_KEY( Real, v_bb_perturbation_rand_fraction ) //The magnitude of the fully random vector added to the perturbation vector, as a fraction of the final perturbation vector size.  Default 0.
OPT_KEY( Boolean, v_rmspenalty ) //Should I weight the terms in the partition function by the rms to the positive state?  False by default.
OPT_KEY( Real, v_rmspenalty_n ) //The exponent to which the RMS penalty should be raised.  Try setting this to the effective number of degrees of freedom of the system.  Set to 1.0 by default.
OPT_KEY( Boolean, v_randseq ) //Should I start with a random sequence?  False by default.  If true, overrides any provided FASTA file.
OPT_KEY( IntegerVector, v_ignoreresidue ) //Should I ignore certain residues when making mutations?
OPT_KEY( IntegerVector, v_ignoreresidue_in_rms) //Should certain residues be excluded from the RMS calculation?
OPT_KEY( Boolean, v_CB_in_rms ) //Should beta carbons be used in the RMS calculation?  True by default.
OPT_KEY( Boolean, v_removelowrms ) //Should I remove negative states with a backbone RMS less than a certain threshhold from the positive state?  True by default.
OPT_KEY( Real, v_removelowrms_threshhold) //RMS threshhold for removing negative states, 0.05 by default.
OPT_KEY( Real, v_lambda) //The bredth of the Gaussian, in Angstroms, that's used to weight states as "positive" for scoring.
OPT_KEY( Boolean, v_use_cartesian_min) //If true, each structure is subjected to a final round of Cartesian minimization.  Default false.
OPT_KEY( Boolean, v_cyclic) //Is this a cyclic peptide?
OPT_KEY( Boolean, v_use_cyclic_permutations) //Should cyclic permutations of backbone dihedral vectors be tried?  Increases number of states n-fold, where n is the number of amino acids.
OPT_KEY( Integer, v_cyclic_permutation_offset) //The size of the offsets, in residues, tried if the v_use_cyclic_permutations flag is activated.  1 by default.
OPT_KEY( Boolean, v_strip_terminal_gly) //Should I strip off terminal glycine residues, if present?  False by default.
OPT_KEY( FileVector, v_PCAfiles) //List of PCA files, in the same order as pdb files.
OPT_KEY( File, v_PCAfiles_list) //A text file with a list of PCA files.  Alternative to v_PCAfiles.
OPT_KEY( File, v_PCAfile_native) //The PCA file for the design (target) state.
OPT_KEY( Boolean, v_allow_gly_mutation) //Should I allow mutations to glycine?  Default true.
OPT_KEY( Boolean, v_preserve_chirality) //If D-amino acids are in the input native structure, should I make sure that I only mutate to D?  Default true.
OPT_KEY( Boolean, v_split_replicates_over_processes) //If true, different replicates of the backbone perturbation of a given structure can be handed to different processes.  False by default.
OPT_KEY( Boolean, v_savememory) //If true, only the master proc stores all states and PCA vectors; these are transmitted as needed to slave procs.  If false, each slave holds the full list.
OPT_KEY( Boolean, v_scoreonly) //If true, scores input sequence and exits.  False by default.
OPT_KEY( IntegerVector, v_disulfide_positions) //List of disulfide-bonded positions.  Empty list by default.
OPT_KEY( String, v_allowed_aminoacids) //List of allowed amino acids, as a string of one-letter codes, no spaces.
OPT_KEY( String, v_allowed_betas) //List of allowed beta amino acids, as a string of one-letter codes, no spaces.
OPT_KEY( Integer, v_PCAlimit) //Max number of PCA vectors stored per state
OPT_KEY( File, v_equivalent_positions) //Equivalent positions file.
OPT_KEY( File, v_covalent_connections) //Covalent connectoins file.
OPT_KEY( File, v_cst_file) //Constraints file
OPT_KEY( String, v_native_tag) //The tag of the structure, read in from silent files, that will be counted as the positive design state.
OPT_KEY( StringVector, v_extra_rms_atoms) //Extra atoms to use in the RMSD calculation.
OPT_KEY( IntegerVector, v_norepack_positions) //Positions that will not be repacked during FastRelax steps.  These will also not be designed (i.e. no mutations will be permitted here).  Default empty list.
OPT_KEY( Boolean, v_consider_all_replicates) //If true, then all replicates are considered in evaluating the objective function that estimates delta G.  If false, just the lowest-energy pose from each replicate is used.  False by default.

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::Size> empty_vector;
	utility::vector1<std::string> empty_filelist;
	NEW_OPT( v_trialcount ,"" , 2000 );
	NEW_OPT( v_relaxrounds ,"" , 1 );
	NEW_OPT( v_MCtemperature ,"The temperature for the application of the Metropolis criterion when searching sequence space. 1.0 by default." , 1.0 );
	NEW_OPT( v_kbt ,"The Boltzmann temperature for estimating fractional population of the positive state (i.e. for scoring each sequence). 1.0 by default." , 1.0 );
	NEW_OPT( v_kbt_minimize, "The Boltzmann temperature for the inner Monte Carlo search of backbone conformational space when minimizing each structure. 1.0 by default.", 1.0);
	NEW_OPT( v_MCminimize_replicates, "", 5);
	NEW_OPT( v_MCminimize_steps, "", 5);
	NEW_OPT( v_bb_perturbation, "", 3.0);
	NEW_OPT( v_bb_perturbation_rand_fraction, "When using PCA-based backbone perturbation, a small, fully random vector can also be added to the perturbation vector.  This value specifies the magnitude of the random vector as a fraction of the magnitude of the total perturbation vector.  Default 0.0; only used if PCA vectors are specified and are available (otherwise the perturbation is fully random).", 0.0);
	NEW_OPT( v_rmspenalty, "", false);
	NEW_OPT( v_rmspenalty_n, "", 1.0 );
	NEW_OPT( v_randseq, "If true, the sequence is randomized at the start (with the exception of positions protected with the -v_ignoreresidue flag).  Default false.", false);
	NEW_OPT( v_ignoreresidue, "Residues to be protected from mutation.  Default empty list.", empty_vector);
	NEW_OPT( v_ignoreresidue_in_rms, "Residues to be ignored when carrying out the RMSD calculation.  Default empty list.", empty_vector);
	NEW_OPT( v_CB_in_rms, "If true, beta carbons are used in calculating RMS values.  True by default.", true);
	NEW_OPT( v_removelowrms, "", true);
	NEW_OPT( v_removelowrms_threshhold, "", 0.1);
	NEW_OPT( v_lambda, "", 1.5);
	NEW_OPT( v_use_cartesian_min, "If true, each structure is subjected to a final round of Cartesian minimization.  Default false.", false);
	NEW_OPT( v_cyclic, "If true, the N- and C-termini are connected by a peptide bond.  Default false.", false);
	NEW_OPT( v_use_cyclic_permutations, "If true, cyclic permutations of the backbone conformations are tried.  Default false.  Requires -v_cyclic true if used.  This increases the number of conformations considered by a factor of N, where N is the number of amino acids in the sequence.", false);
	NEW_OPT( v_cyclic_permutation_offset, "The size of the shift, in residues, that will be used when trying cyclic permutations if the -v_use_cyclic_permutations flag is specified.  1 by default.", 1);
	NEW_OPT( v_strip_terminal_gly, "If true, terminal glycines (if present) are removed.  False by default.", false);
	NEW_OPT( v_PCAfiles, "The list of optional PCA files.  If used, these must be specified in the order corresponding to states specified using the -in:file:s flag.  Default empty list (i.e. default unused).", empty_filelist);
	NEW_OPT( v_PCAfiles_list, "A text file with a list of optional PCA files, as an alternative to the v_PCAfiles option.  If used, these must be specified in the order corresponding to states specified using the -in:file:s or -in:file:l flag.  Default unused.", "");
	NEW_OPT( v_PCAfile_native, "The PCA file corresponding to the native state.  Must be specified if and only if PCA files are specified for other states.  Defaults to nothing (i.e. defaults to unused).", "");
	NEW_OPT( v_allow_gly_mutation, "If true, mutations to glycine are permitted.  True by default.", true);
	NEW_OPT( v_preserve_chirality, "If true, any D-amino acid positions in the input native structure remain D-amino acids when mutated.  Default true.  (Otherwise, all mutations are to L-amino acids).", true);
	NEW_OPT( v_split_replicates_over_processes, "If true, different replicates of the backbone perturbation of a given structure can be handed to different processes.  This is useful on massively parallel architectures.  True by default.", true);
	NEW_OPT( v_savememory, "If true, only the master proc stores all states and PCA vectors; these are transmitted as needed to slave procs.  If false, each slave holds the full list.  True by default.", true);
	NEW_OPT( v_scoreonly, "If true, this scores the input sequence and then exits.  False by default.", false);
	NEW_OPT( v_disulfide_positions, "A list of positions that are disulfide-bonded, with pairs representing disulfide-bonded residues.  For example, -v_disulfide_positions 5 15 9 23 would mean that C5 and C15 are disulfide-linked, and C9 and C23 are disulfide-linked.  Listed positions must be cysteines, and must be in the -v_ignoreresidue list!  Default empty list.", empty_vector);
	NEW_OPT( v_allowed_aminoacids, "A list of the allowed amino acids (1-letter codes, no spaces).  If this option is used, it overrides v_allow_gly_mutation.  The 19 standard amino acids (minus cysteine) are used if this is not specified.", "");
	NEW_OPT( v_allowed_betas, "A list of the allowed beta-amino acids (1-letter codes, no spaces).  If this option is used, it overrides v_allow_gly_mutation.  The 19 standard L-beta-3-amino acids (minus cysteine) are used if this is not specified.", "");
	NEW_OPT( v_PCAlimit, "If specified, this is the maximum number of PCA vectors stored per state.  This can save a fair bit of memory.  Not used unless a value is specified (i.e. by default, all PCA vectors are stored).", 0);
	NEW_OPT( v_equivalent_positions, "An optional file that lists equivalent positions in the protein, which, when mutated, are all mutated together to the same amino acid residue.  This ASCII file must consist of a series of lines, each of which consists of a series of space-separated integers.  Each line represents a single set of equivalent positions.  Not used if not specified.", "");
	NEW_OPT( v_covalent_connections, "An optional file that lists covalent connections between residues (one per line, in the format \"resnum1 atomname1 resnum2 atomname2\").  Not used if not specified.", "");
	NEW_OPT( v_cst_file, "An optional constraints file applied to all structures.  If used, constraints weights are automatically set to 1.0 unless otherwise specified in the weights file.  Not used if not specified.", "");
	NEW_OPT( v_native_tag, "An alternative to -in:file:native for specifying the native state.  If this option is used with a string corresponding to a silent file tag, the silent file pose with that tag is used as the positive design state.  Default unused.", "");
	NEW_OPT ( v_extra_rms_atoms, "A list of additional atoms to use in the RMSD calculation, in the format \"residue:atomname residue:atomname residue:atomname\".  For example, \"-v_extra_rms_atoms 7:SG 12:CG 12:CD 12:CE 12:NZ 14:OG\".  Default empty list.", "");
	NEW_OPT( v_norepack_positions, "Positions that will not be repacked during FastRelax steps.  These will also not be designed (i.e. no mutations will be permitted here).  Default empty list.", empty_vector);
	NEW_OPT( v_consider_all_replicates, "If true, then all replicates are considered in evaluating the objective function that estimates delta G.  If false, just the lowest-energy pose from each replicate is used.  False by default.", false);

}

//Functions to set the backbone dihedral angles of beta-amino acid peptides:
void betapeptide_setphi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 1 ), angle );
	return;
}

void betapeptide_settheta (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 2 ), angle );
	return;
}

void betapeptide_setpsi (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 3 ), angle );
	return;
}

void betapeptide_setomega (
	core::pose::Pose &pose,
	core::Size resnumber,
	core::Real angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 4 ), angle );
	return;
}

//Function to return beta-amino acid names given a one-letter code
void get_beta_name(
	const char onelettercode,
	std::string &threelettercode
) {
	threelettercode = "B3";
	threelettercode+=onelettercode;
	return;
}

//Function to check whether a value is in a vector
bool is_in_list (
	const core::Size val,
	const utility::vector1 < core::Size > &vect
) {
	core::Size vsize = vect.size();
	if(vsize>0) {
		for(core::Size i=1; i<=vsize; i++) {
			if(vect[i]==val) return true;
		}
	}
	return false;
}

//Function to check whether a value is in a list of vectors, and return which entry in the outer list contains the value if it is.  Returns 0 otherwise.
core::Size is_in_list_of_lists (
	const core::Size val,
	const utility::vector1 < utility::vector1 < core::Size > > &vect
) {
	core::Size vsize = vect.size();
	if(vsize>0) {
		for(core::Size i=1; i<=vsize; i++) {
			if(is_in_list(val, vect[i])) return i;
		}
	}
	return 0;
}

//Function to return D-amino acid names given a one-letter code
void get_D_name(
	const char onelettercode, //Input
	std::string & threelettercode //Output
) {
	switch (onelettercode) {
		case 'A':
			threelettercode="DALA";
			break;
		case 'C':
			threelettercode="DCYS";
			break;
		case 'D':
			threelettercode="DASP";
			break;
		case 'E':
			threelettercode="DGLU";
			break;
		case 'F':
			threelettercode="DPHE";
			break;
		case 'G':
			threelettercode="GLY"; //Keep this regular glycine for now.  At some point, I might add a "D-glycine" for glycine in the context of a bunch of D-amino acids.
			break;
		case 'H':
			threelettercode="DHIS";
			break;
		case 'I':
			threelettercode="DILE";
			break;
		case 'K':
			threelettercode="DLYS";
			break;
		case 'L':
			threelettercode="DLEU";
			break;
		case 'M':
			threelettercode="DMET";
			break;
		case 'N':
			threelettercode="DASN";
			break;
		case 'P':
			threelettercode="DPRO";
			break;
		case 'Q':
			threelettercode="DGLN";
			break;
		case 'R':
			threelettercode="DARG";
			break;
		case 'S':
			threelettercode="DSER";
			break;
		case 'T':
			threelettercode="DTHR";
			break;
		case 'V':
			threelettercode="DVAL";
			break;
		case 'W':
			threelettercode="DTRP";
			break;
		case 'Y':
			threelettercode="DTYR";
			break;
		default:
			threelettercode="UNK";
			break;
	}

	return;
}

void parse_extra_atom_list ( utility::vector1 < core::id::NamedAtomID > &extra_atom_list )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!option[v_extra_rms_atoms].user()) return; //Do nothing if there's no list of extra atoms.
	extra_atom_list.clear();

	utility::vector1 <std::string> tags = option[v_extra_rms_atoms]();

	//Parse each tag:
	//printf("The following additional atoms will be used in the RMSD calculation:\n");
	for(core::Size i=1, imax=tags.size(); i<=imax; ++i) { //Loop through each tag
		core::Size colonposition = tags[i].find( ':' );
		std::string resstring = tags[i].substr( 0, colonposition );
		core::Size res = static_cast<core::Size>( atoi( resstring.c_str() ) ); //The residue number
		std::string atomname = tags[i].substr( colonposition + 1); //The atom name
		//printf("\tResidue %lu, Atom %s\n", res, atomname.c_str());
		extra_atom_list.push_back( core::id::NamedAtomID( atomname, res ) );
	}

	//printf("\n");

	return;
}

//Function to determine whether to use an atom or not in the RMS calculation.
//Currently, this is set to use backbone heavy atoms and CB atoms (if present), with the exception of CB in the terminal cysteines.
bool use_in_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size atomno,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!option[v_cyclic]() && resno==pose1.n_residue() && pose1.residue(resno).has("O") && atomno==pose1.residue(resno).atom_index( "O")) return false;

	if(option[v_ignoreresidue_in_rms]().size()>0) { //Return false for residues that are to be ignored.
		for(core::Size i=1; i<=option[v_ignoreresidue_in_rms]().size(); i++) {
			if((core::Size)option[v_ignoreresidue_in_rms]()[i]==resno) return false;
		}
	}

	if(option[v_norepack_positions]().size()>0) { //If positions that can't repack have been specified, these can't be mutated either.
		for(core::Size i=1, imax=option[v_norepack_positions]().size(); i<=imax; ++i) {
			if((core::Size)option[v_norepack_positions]()[i]==resno) return false;
		}
	}
	
	if(pose1.residue(resno).has( "N") && pose2.residue(resno).has( "N") && pose1.residue(resno).atom_index( "N")==atomno) return true;
	if(pose1.residue(resno).has("CA") && pose2.residue(resno).has("CA") && pose1.residue(resno).atom_index("CA")==atomno) return true;
	if(pose1.residue(resno).has("CM") && pose2.residue(resno).has("CM") && pose1.residue(resno).atom_index("CM")==atomno) return true; //Beta-3-amino acids
	if(pose1.residue(resno).has( "C") && pose2.residue(resno).has( "C") && pose1.residue(resno).atom_index( "C")==atomno) return true;
	if(pose1.residue(resno).has( "O") && pose2.residue(resno).has( "O") && pose1.residue(resno).atom_index( "O")==atomno) return true;
	if(option[v_CB_in_rms]() && pose1.residue(resno).has("CB") && pose2.residue(resno).has("CB") && pose1.residue(resno).atom_index("CB")==atomno) return true;

	//Check extra atoms:
	for(core::Size i=1, imax=extra_atom_list.size(); i<=imax; ++i) {
		if(resno!=extra_atom_list[i].rsd()) continue;
		if(pose1.residue(resno).has( extra_atom_list[i].atom() ) && pose2.residue(resno).has( extra_atom_list[i].atom() ) && pose1.residue(resno).atom_index( extra_atom_list[i].atom() )==atomno) return true;
	}

	return false;
}

//Function that returns the backbone and CB RMS between two poses.
core::Real get_distance_measure(
		const core::pose::Pose & pose1,
		const core::pose::Pose & pose2,
		utility::vector1 <core::id::NamedAtomID> const &extra_rms_atoms
) {
	std::map< core::id::AtomID, core::id::AtomID > amap;
	for(int ir = 1; ir <= (int)pose1.n_residue(); ++ir) {
		for(int ia = 1; ia <= (int)pose1.residue(ir).nheavyatoms(); ++ia) {
			if(use_in_rmsd(pose1,pose2,ir,ia,extra_rms_atoms)) {
				using core::id::AtomID;
				amap.insert(std::pair< core::id::AtomID, core::id::AtomID> (AtomID(ia,ir), AtomID(ia, ir)));
			}
		}
	}
	return core::scoring::rms_at_all_corresponding_atoms(pose1,pose2,amap);
}

//Function to calculate the score of a given sequence, based on the energies of the sequence in all possible backbone conformations:
core::Real calc_sequence_score (utility::vector1 < core::Real > &energies, utility::vector1 < core::Real > &rmsvalues)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real scorenumerator = 0.0;
	core::Real scoredenominator = 0.0;

	const core::Real kbt = option[v_kbt]();

	for(core::Size i=1; i<=energies.size(); i++) {
		scorenumerator += exp(-1.0*pow(rmsvalues[i]/option[v_lambda](),2.0))*exp(-1.0*energies[i]/kbt);
		if (option[v_rmspenalty]()) scoredenominator += pow((1.0+rmsvalues[i]), option[v_rmspenalty_n]())*exp(-1.0*energies[i]/kbt);
		else scoredenominator += exp(-1.0*energies[i]/kbt);
	}
	
	return -kbt*log(scorenumerator/(scoredenominator-scorenumerator));
}

//Function to set the termini of a non-cyclic peptide:
void set_up_termini (core::pose::Pose &pose) {
	//printf("Setting up termini.\n"); fflush(stdout);
	//Loop through and add lower termini and upper termini appropriately:
	for(core::Size ir=1, irmax=pose.n_residue(); ir<=irmax; ++ir) {
		if(pose.residue(ir).is_polymer() && !pose.residue(ir).has_lower_connect() )
			core::pose::add_lower_terminus_type_to_pose_residue(pose, ir);
		if(pose.residue(ir).is_polymer() && !pose.residue(ir).has_upper_connect() )
			core::pose::add_upper_terminus_type_to_pose_residue(pose, ir);
	}
	//printf("Termini setup complete.\n"); fflush(stdout); //DELETE ME
	return;
}

//Function to add user-specified constraints (specified with a CST file) to a pose:
void add_user_constraints (
	core::pose::Pose &mypose
) {
	using namespace protocols::simple_moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(!option[v_cst_file].user()) return; //Do nothing if no constraint file is specified.

	ConstraintSetMoverOP cst_maker = new ConstraintSetMover();
	std::string cstfile = option[v_cst_file]();

	cst_maker->constraint_file(cstfile);
	cst_maker->add_constraints(true); //Add constraints to anything else already there.
	cst_maker->apply(mypose);	

	return;
}

//Function to add cyclic constraints to a pose:
void addcyclicconstraints (core::pose::Pose &mypose) {
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.n_residue());
	
	const std::string hstring = (mypose.residue(1).name1()=='P' ? "CD" : "H");
	const core::Size hindex = mypose.residue(1).atom_index(hstring);
	
	//Rebuild the N-terminal proton.  This has to be done in a slightly irritating way because Rosetta doesn't really like the fact
	//that the last residue is connected to the first:
	{
		core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::pose::Pose dialanine;
		std::string diala_seq = "";
		if(mypose.residue(mypose.n_residue()).has("CM")) diala_seq+="A[B3A]"; else diala_seq+="A";
		if(mypose.residue(1).has("CM")) diala_seq+="A[B3A]"; else diala_seq+="A";

		{//Make the pose:
			core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			core::chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( diala_seq, *standard_residues, false );
			for(core::Size ir=1, numres=requested_types.size(); ir<=numres; ir++){
				core::chemical::ResidueType const & rsd_type = *requested_types[ ir ];
				core::conformation::ResidueOP new_rsd = core::conformation::ResidueFactory::create_residue( rsd_type );
				if(ir==1) dialanine.append_residue_by_jump( *new_rsd, 1, "", "", true );
				else dialanine.append_residue_by_bond(*new_rsd, true);
			}
		}
		//core::pose::make_pose_from_sequence(dialanine, "AA", *standard_residues, true); //The termini are OPEN.

		core::Real omegaval = numeric::dihedral_degrees(
				mypose.residue(mypose.n_residue()).xyz( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ),
				mypose.residue(mypose.n_residue()).xyz("C"),
				mypose.residue(1).xyz("N"),
				mypose.residue(1).xyz("CA")
			);

		if(dialanine.residue(1).has("CM")) betapeptide_setomega(dialanine, 1, omegaval);
		else dialanine.set_omega(1, omegaval);
	
		core::id::AtomID_Map< core::id::AtomID > amap;
		core::pose::initialize_atomid_map(amap,dialanine,core::id::BOGUS_ATOM_ID);
		amap[AtomID(dialanine.residue(1).atom_index( (dialanine.residue(1).has("CM")?"CM":"CA") ),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index( (dialanine.residue(1).has("CM")?"CM":"CA") ),mypose.n_residue());
		amap[AtomID(dialanine.residue(1).atom_index("C"),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index("C"),mypose.n_residue());
		amap[AtomID(dialanine.residue(1).atom_index("O"),1)] = AtomID(mypose.residue(mypose.n_residue()).atom_index("O"),mypose.n_residue());
		amap[AtomID(dialanine.residue(2).atom_index("N"),2)] = AtomID(mypose.residue(1).atom_index("N"),1);
		amap[AtomID(dialanine.residue(2).atom_index("CA"),2)] = AtomID(mypose.residue(1).atom_index("CA"),1);
		core::scoring::superimpose_pose( dialanine, mypose, amap );

		mypose.conformation().set_xyz(AtomID(hindex, 1), dialanine.residue(2).xyz("H")); //This won't be right for the proline side-chain, but it will be rebuilt by a repack in any case.
	}

	mypose.conformation().declare_chemical_bond(1, "N", mypose.n_residue(), "C"); //Declare a chemical bond between the N and C termini.

	{//Peptide bond length constraint:
		FuncOP harmfunc1 = new HarmonicFunc( 1.3288, 0.02);
		ConstraintCOP distconst1 = new AtomPairConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ) ,
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				harmfunc1
			);
		mypose.add_constraint (distconst1);
	}

	{	//Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1 = new CircularHarmonicFunc( PI, 0.02);
		ConstraintCOP dihedconst1 = new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index(hstring) , 1) ,
				circharmfunc1
			);
		ConstraintCOP dihedconst2 = new DihedralConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ) , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				circharmfunc1
			);

		mypose.add_constraint (dihedconst1);
		mypose.add_constraint (dihedconst2);
	}

	{	//Peptide bond angle constraints:
		FuncOP circharmfunc2a = new CircularHarmonicFunc( CNCa_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2b = new CircularHarmonicFunc( CNH_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2c = new CircularHarmonicFunc( CaCN_ANGLE/180.0*PI, 0.02);
		FuncOP circharmfunc2d = new CircularHarmonicFunc( OCN_ANGLE/180.0*PI, 0.02);

		ConstraintCOP angleconst1 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index("CA") , 1) ,
				circharmfunc2a
			);
		ConstraintCOP angleconst2 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				AtomID( mypose.residue(1).atom_index(hstring) , 1) ,
				circharmfunc2b
			);
		ConstraintCOP angleconst3 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ) , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				circharmfunc2c
			);
		ConstraintCOP angleconst4 = new AngleConstraint (
				AtomID( mypose.residue(mypose.n_residue()).atom_index("O") , mypose.n_residue() ),
				AtomID( mypose.residue(mypose.n_residue()).atom_index("C") , mypose.n_residue() ),
				AtomID( mypose.residue(1).atom_index("N") , 1) ,
				circharmfunc2d
			);

		mypose.add_constraint (angleconst1);
		mypose.add_constraint (angleconst2);
		mypose.add_constraint (angleconst3);
		mypose.add_constraint (angleconst4);
	}

	return;
}

/*Function to jitter the backbone and relax it.
This function does v_repeats number of short Monte Carlo trajectories (of length v_MCsteps), then returns the lowest-energy structure found.
Currently, it does naive backbone perturbations, changing each backbone dihedral by a random amount before relaxing.  Expensive FastRelax steps
are done for every MC step (v_repeats*v_MCsteps FastRelax calls per call of the function).*/
//Altered on 18 June 2013 to use PCA vectors.
void perturb_bb_and_relax (
	const core::Size replicates, //Number of replicates (MC trajectories, each starting at the input pose)
	const core::Size MCsteps, //Number of Monte Carlo steps per trajectory
	core::pose::Pose & pose, //The input pose (changed by the function)
	utility::vector1 < std::pair < core::Size, core::Size >  > &disulf_positions,
	numeric::random::RandomGenerator & RandGen, //A random number generator
	const core::Real bb_perturbation, //The amplitude of the backbone perturbation
	protocols::relax::FastRelax & frlx, //A fastrelax mover
	core::scoring::ScoreFunctionOP const & sfxn, //The owning pointer for the scorefunction
	const core::Real MCkbt, //The Boltzmann temperature for the search
	const utility::vector1 <core::Real> &PCAvariances, //The PCA variances, normalized to the largest element
	const utility::vector1 <utility::vector1 <core::Real> > &PCAvectors //The PCA vectors, with each normalized
) {
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	//printf("starting perturb_bb_and_relax\n"); fflush(stdout); //DELETE ME

	const core::Real bb_perturbation_rand_fraction = option[v_bb_perturbation_rand_fraction]();

	//Variables for replicates:
	core::pose::Pose lowestEpose; //Lowest energy pose of any REPLICATE

	//Variables for MC search:
	core::pose::Pose lastpose_MC; //Last pose accepted in the current MC trajectory
	core::pose::Pose lowestEpose_MC; //Lowest energy pose encountered in the current MC trajectory
	core::Size MCcounter = 0; //The counter for Monte Carlo trajectory steps.

	//Variables and objects for Cartesian minimization:
	core::Real const cartweight = sfxn->get_weight(cart_bonded);
	core::Real const procloseweight = sfxn->get_weight(pro_close);
	core::optimization::MinimizerOptions minoptions("dfpmin_armijo_nonmonotone", 0.000001, true, false, false);
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_bb_true_range(1, pose.n_residue());
	mm->set_chi_true_range(1, pose.n_residue());
	mm->set_jump(false);
	core::optimization::CartesianMinimizer cminimizer;

	for(core::Size replicate = 1; replicate <= replicates; replicate++) {
		MCcounter = 0;
		lastpose_MC = pose;

		do { //Looping through steps in the Monte Carlo trajectory.  The following code executes at least once.
			core::pose::Pose temppose;
				
			//Generate a random vector:
			utility::vector1<core::Real> pertvect; //The perturbation vector
			utility::vector1<core::Real> randvect; //An entirely random vector
			core::Real accumulator = 0.0;
			core::Real gausscoeff;

			//printf("starting perturbation\n"); fflush(stdout); //DELETE ME

			while(true) { //Loop until the perturbation does NOT cause atoms to overlap
				randvect.clear();
				pertvect.clear();
				gausscoeff = RandGen.gaussian();
				accumulator = 0.0;
				temppose = lastpose_MC; //Last pose is the input pose if MC is off, or the last pose accepted for MC moves
				
				//Generate a random perturbation vector:
				for(core::Size ir=1, nres=temppose.n_residue(), counter=1; ir<=nres; ir++) { //Loop through all residues.
					if(!temppose.residue(ir).type().is_alpha_aa() && !temppose.residue(ir).type().is_beta_aa()) continue; //Skip non-amino acids for now
					core::Size ntors=3; //Number of backbone torsions for this residue.
					if(temppose.residue(ir).is_lower_terminus() || ir==1 || temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().lower_connect_id())!=ir-1 ) ntors--; //if this is a lower terminus OR a noncanonical connection, remove phi
					if(temppose.residue(ir).is_upper_terminus() || temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().upper_connect_id())!=ir+1) ntors-=2; //if this is an upper terminus OR a noncanonical connection, remove psi and omega
					if(temppose.residue(ir).type().is_beta_aa()) ntors++;

					for(core::Size itors=1; itors<=ntors; itors++) { //Loop through this residue's torsions
						randvect.push_back(RandGen.gaussian());
						if(itors == ntors && ir<nres) randvect[counter] *= 0.025; //Omega perturbations should be SMALL.
						accumulator+=pow(randvect[counter],2); //At the end of the outer for loop, accumulator will be the square of the perturbation vector
						counter++;
					} //End loop through this residue's torsions
				} //End loop through all residues
				if(option[v_cyclic]()) { //Three extra ones for true cyclic peptides
					for(core::Size j=1; j<=3; j++) {
						randvect.push_back(RandGen.gaussian());
						if(j==2) randvect[randvect.size()] *= 0.025; //Omega perturbations should be SMALL.
						accumulator+=pow(randvect[randvect.size()],2); //At the end of the outer for loop, accumulator will be the square of the perturbation vector
					}
				}

				//Normalize the perturbation vector:
				accumulator = sqrt(accumulator);
				for(core::Size i=1; i<=randvect.size(); i++) randvect[i] *= (core::Real)gausscoeff/accumulator;

				//printf("normalized accumulator\n"); fflush(stdout); //DELETE ME

				if(PCAvariances.size()>0) { //If there are PCA vectors
					accumulator=0.0; //Reset this for reuse.
					gausscoeff=RandGen.gaussian(); //Pick a new random coefficient.
					utility::vector1 <core::Real> PCApertvect; //Perturbation from PCA vectors.
					PCApertvect.resize(PCAvectors[1].size(), 0.0); //Initialize to an array of zeroes.

					for(core::Size i=1; i<=PCAvariances.size(); i++) { //Loop through all variances (one for each PCA vector)
						core::Real coeff=RandGen.gaussian();
						for(core::Size j=1; j<=PCApertvect.size(); j++) {
							PCApertvect[j]+=coeff*PCAvariances[i]*PCAvectors[i][j];
							if(i==PCAvariances.size()) accumulator+=pow(PCApertvect[j],2); //If this is the last time through this loop, figure out the magnitude of the perturbation vector.
						}
					}
					//Normalize:
					accumulator=sqrt(accumulator);
					pertvect.resize(randvect.size());
					for(core::Size i=1; i<=PCApertvect.size(); i++) {
						PCApertvect[i] *= (core::Real)gausscoeff/accumulator;
						pertvect[i] = bb_perturbation_rand_fraction*randvect[i] + (1.0-bb_perturbation_rand_fraction)*PCApertvect[i];
					}
				} else pertvect = randvect; //If there are no PCA vectors, just use the random vector.
				//printf("generated PCA perturbation vector\n"); fflush(stdout); //DELETE ME

				//printf("pertvect_size=%lu, randvect_size=%lu\n", pertvect.size(), randvect.size()); fflush(stdout);  //DELETE ME!

				//Scale pertvect by bb_perturbation:
				for(core::Size i=1; i<=pertvect.size(); i++) pertvect[i] *= bb_perturbation;

				//Apply the perturbation
				for(core::Size ir=1, nres=temppose.n_residue(), counter=1; ir<=nres; ir++) {
					if(temppose.residue(ir).type().is_beta_aa()) { //beta-amino acid
						if(!temppose.residue(ir).is_lower_terminus() && ir>1 && temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().lower_connect_id())==ir-1 ) {
							betapeptide_setphi(temppose, ir, pertvect[counter++]+temppose.residue(ir).mainchain_torsion(1)); //phi
						}
						betapeptide_settheta(temppose, ir, pertvect[counter++]+temppose.residue(ir).mainchain_torsion(2)); //theta
						if(!temppose.residue(ir).is_upper_terminus() && temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().upper_connect_id())==ir+1) {
							betapeptide_setpsi(temppose, ir, pertvect[counter++]+temppose.residue(ir).mainchain_torsion(3)); //psi
							betapeptide_setomega(temppose, ir, pertvect[counter++]+temppose.residue(ir).mainchain_torsion(4)); //omega
						}
					} else if(temppose.residue(ir).type().is_alpha_aa()) { //alpha-amino acid
						if(!temppose.residue(ir).is_lower_terminus() && ir>1 && temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().lower_connect_id())==ir-1 ) {
							temppose.set_phi(ir, pertvect[counter++]+temppose.phi(ir));
						}
						if(!temppose.residue(ir).is_upper_terminus() && temppose.residue(ir).connected_residue_at_resconn(temppose.residue(ir).type().upper_connect_id())==ir+1) {
							temppose.set_psi(ir, pertvect[counter++]+temppose.psi(ir));
							temppose.set_omega(ir, pertvect[counter++]+temppose.omega(ir));
						}
					} else {} //default case -- do nothing
				}

				///printf("applied perturbation\n"); fflush(stdout); //DELETE ME

				//Additional perturbations for cyclic peptides:
				if(option[v_cyclic]()) {
					addcyclicconstraints(temppose);
					core::Real ang1 = numeric::dihedral_degrees (
							temppose.residue(temppose.n_residue()).xyz( (temppose.residue(temppose.n_residue()).has("CM")?"CA":"N") ),
							temppose.residue(temppose.n_residue()).xyz( (temppose.residue(temppose.n_residue()).has("CM")?"CM":"CA") ),
							temppose.residue(temppose.n_residue()).xyz("C"),
							temppose.residue(temppose.n_residue()).xyz("O")
						);
					core::Real ang2 = numeric::dihedral_degrees (
							temppose.residue(1).xyz((temppose.residue(1).name1()=='P'?"CD":"H")),
							temppose.residue(1).xyz("N"),
							temppose.residue(1).xyz("CA"),
							temppose.residue(1).xyz( (temppose.residue(1).has("CM")?"CM":"C") )
						);
				temppose.conformation().set_torsion_angle(//Set psi of last residue
					core::id::AtomID(temppose.residue(temppose.n_residue()).atom_index( (temppose.residue(temppose.n_residue()).has("CM")?"CA":"N") ), temppose.n_residue()),
					core::id::AtomID(temppose.residue(temppose.n_residue()).atom_index( (temppose.residue(temppose.n_residue()).has("CM")?"CM":"CA") ), temppose.n_residue()),
					core::id::AtomID(temppose.residue(temppose.n_residue()).atom_index("C"), temppose.n_residue()),
					core::id::AtomID(temppose.residue(temppose.n_residue()).atom_index("O"), temppose.n_residue()),
					(ang1+pertvect[pertvect.size()-2])/180.0*PI);
				temppose.conformation().set_torsion_angle( //Set phi of first residue
					core::id::AtomID(temppose.residue(1).atom_index((temppose.residue(1).name1()=='P'?"CD":"H")), 1),
					core::id::AtomID(temppose.residue(1).atom_index("N"), 1),
					core::id::AtomID(temppose.residue(1).atom_index("CA"), 1),
					core::id::AtomID(temppose.residue(1).atom_index( (temppose.residue(1).has("CM")?"CM":"C") ), 1),
					(ang2+pertvect[pertvect.size()])/180.0*PI);
				}

				temppose.update_residue_neighbors();

				//printf("applied cyclic perturbation\n"); fflush(stdout); //DELETE ME

				bool goodperturbation = true;
				if(option[v_cyclic]()) { //If this is a cyclic peptide, check the length of the peptide bond linking the terminal residues.
					const numeric::xyzVector <core::Real> v1 = temppose.residue(1).atom("N").xyz();
					const numeric::xyzVector <core::Real> v2 = temppose.residue(temppose.n_residue()).atom("C").xyz();
					if( sqrt(pow(v1[0]-v2[0], 2.0) + pow(v1[1]-v2[1], 2.0) + pow(v1[2]-v2[2], 2.0)) < 1.0e-5) goodperturbation = false; //Bad perturbation if peptide bond length is zero.
				}
				for(core::Size ir=1; (ir<=temppose.n_residue()-1 && goodperturbation); ir++) { //Check disulfide lengths
					for(core::Size jr=ir; (jr<=temppose.n_residue() && goodperturbation); jr++) {
						if(core::conformation::is_disulfide_bond(temppose.conformation(), ir, jr) ) {
							const numeric::xyzVector <core::Real> v1 = temppose.residue(ir).atom("SG").xyz();
							const numeric::xyzVector <core::Real> v2 = temppose.residue(jr).atom("SG").xyz();
							if( sqrt(pow(v1[0]-v2[0], 2.0) + pow(v1[1]-v2[1], 2.0) + pow(v1[2]-v2[2], 2.0)) < 1.0e-5) goodperturbation = false; //Bad perturbation if the sulphurs of two disulfide-linked cysteines overlap.
						}
					}
				}
				if(goodperturbation) break;
			}


			if(disulf_positions.size() > 0) temppose.conformation().fix_disulfides(disulf_positions);
			//if(option[v_cyclic]()) addcyclicconstraints(temppose);
			//printf("about to relax\n"); fflush(stdout); //DELETE ME
			//temppose.dump_pdb("temp.pdb"); exit(1); //DELETE ME
			//protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn); //Create the RepackSidechains mover and set the score function -- DELETE ME
			//repack_sc.apply(temppose); --DELETE ME

			frlx.apply(temppose); //Relax the pose copy.

			//printf("relaxed\n"); fflush(stdout); //DELETE ME

			if(option[v_use_cartesian_min]()) { //If we're doing a final round of Cartesian minimization.
				sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.
				sfxn->set_weight(pro_close, 0.0); //Turn off pro_close.
				cminimizer.run(temppose, *mm, *sfxn, minoptions);
				sfxn->set_weight(cart_bonded, cartweight); //Reset cart_bonded.
				sfxn->set_weight(pro_close, procloseweight); //Reset pro_close.
			}
			(*sfxn)(temppose); //Score the relaxed pose copy without cart_bonded.

			if(MCsteps > 1) { //If we're doing more than one Monte Carlo round
				if(MCcounter == 0) { //If this is the first move, initialize lastpose_MC and lowestEpose_MC:
					lastpose_MC = temppose;
					lowestEpose_MC = temppose;
				} else { //Else if we're past the first MC move, do MC tests:
					if(temppose.energies().total_energy() < lastpose_MC.energies().total_energy()) {
						//Accept the move if the energy drops:
						lastpose_MC = temppose;
						//If this is the lowest energy encountered in the MC search, store the pose:
						if(lastpose_MC.energies().total_energy() < lowestEpose_MC.energies().total_energy())
							lowestEpose_MC = lastpose_MC;
					} else {
						//Apply the Metropolis criterion if the energy goes up:
						if(exp(-1.0*( temppose.energies().total_energy() - lastpose_MC.energies().total_energy() )/MCkbt) > RandGen.uniform())
							lastpose_MC = temppose; //Accept the move (no else required to reject it).
					}
				}
			}
			else lowestEpose_MC = temppose;

			MCcounter++; //Increment the number of Monte Carlo moves
		} while (MCcounter < MCsteps); //Loop back to the "do" if we haven't done MCsteps Monte Carlo moves.

		if(replicate == 1 || lowestEpose_MC.energies().total_energy() < lowestEpose.energies().total_energy()) lowestEpose = lowestEpose_MC;
	}

	pose = lowestEpose;

	return;
}

//Function to transmit x, y, z coordinates of all atoms to ONE other proc.
//This assumes that each proc contains a pose with identical SEQUENCE.
void vMPI_transmitpose(
	core::pose::Pose &pose,
	int thisproc,
	int sendingproc,
	int receivingproc
	//core::scoring::ScoreFunctionOP sfxn
) {
	if((thisproc!=receivingproc)&&(thisproc!=sendingproc)) return; //Don't do anything if this is neither the sending nor the receiving proc.

	numeric::xyzVector <double> atomxyz;
	MPI_Status mpistatus;

	//Loop through all residues:
	for(core::Size ir=1; ir<=pose.n_residue(); ir++) {
		//Loop through all atoms in the residue
		for(core::Size ia=1; ia<=pose.residue(ir).natoms(); ia++) {
			core::id::AtomID curatom(ia, ir);
			if(thisproc==sendingproc) {
				atomxyz= pose.xyz(curatom);
				MPI_Send(&atomxyz[0], 3, MPI_DOUBLE, receivingproc, 0, MPI_COMM_WORLD);
			} else if (thisproc==receivingproc) {
				MPI_Recv(&atomxyz[0], 3, MPI_DOUBLE, sendingproc, 0, MPI_COMM_WORLD, &mpistatus);
				pose.set_xyz(curatom, atomxyz);
			}
		}
	}
	//if (thisproc==receivingproc) (*sfxn)(pose); //Score the transmitted pose if appropriate.
	return;
}

//Function to transmit x, y, z coordinates of all atoms to all other procs.
//This assumes that each proc contains a pose with identical SEQUENCE.
void vMPI_Bcastpose(
	core::pose::Pose &pose,
	int rootproc,
	int thisproc
	//core::scoring::ScoreFunctionOP sfxn
) {
	numeric::xyzVector <double> atomxyz;
	
	//Loop through all residues:
	for(core::Size ir=1; ir<=pose.n_residue(); ir++) {
		//Loop through all atoms in the residue
		for(core::Size ia=1; ia<=pose.residue(ir).natoms(); ia++) {
			core::id::AtomID curatom(ia, ir);
			if(thisproc==rootproc) {
				atomxyz= pose.xyz(curatom);
			}
			MPI_Bcast(&atomxyz[0], 3, MPI_DOUBLE, rootproc, MPI_COMM_WORLD);
			if(thisproc!=rootproc) {
				pose.set_xyz(curatom, atomxyz);
			}
		}
	}
	//if (thisproc!=rootproc) (*sfxn)(pose); //Score the transmitted pose if appropriate.
	return;
}

//Function to extract the phi, psi, and omega angles from a pose and store them in a core::Real vector
void store_backbone (
	const core::pose::Pose &mypose,
	utility::vector1 < core::Real > &outvect
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	outvect.clear(); //Empty the output vector

	for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++) {
		if(mypose.residue(ir).type().is_beta_aa()) { //beta-amino acid
			if(!mypose.residue(ir).is_lower_terminus()) outvect.push_back(mypose.residue(ir).mainchain_torsion(1)); //phi
			outvect.push_back(mypose.residue(ir).mainchain_torsion(2)); //theta
			if(!mypose.residue(ir).is_upper_terminus()) {
				outvect.push_back(mypose.residue(ir).mainchain_torsion(3)); //psi
				outvect.push_back(mypose.residue(ir).mainchain_torsion(4)); //omega
			}
		} else if (mypose.residue(ir).type().is_alpha_aa()) { //alpha-amino acid
			if(!mypose.residue(ir).is_lower_terminus()) outvect.push_back(mypose.phi(ir));
			if(!mypose.residue(ir).is_upper_terminus()) {
				outvect.push_back(mypose.psi(ir));
				outvect.push_back(mypose.omega(ir));
			}
		} else {} //Default case -- do nothing
	}

	if (option[v_cyclic]) {
		outvect.push_back(//Psi of last residue
			numeric::dihedral_degrees( 
				mypose.residue(mypose.n_residue()).xyz( (mypose.residue(mypose.n_residue()).has("CM")?"CA":"N") ),
				mypose.residue(mypose.n_residue()).xyz( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ),
				mypose.residue(mypose.n_residue()).xyz("C"),
				mypose.residue(1).xyz("N")
			)
		);
		outvect.push_back(//Omega between last residue and first residue
			numeric::dihedral_degrees(
				mypose.residue(mypose.n_residue()).xyz( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ),
				mypose.residue(mypose.n_residue()).xyz("C"),
				mypose.residue(1).xyz("N"),
				mypose.residue(1).xyz("CA")
			)
		);
		outvect.push_back(//Phi of first residue
			numeric::dihedral_degrees(
				mypose.residue(mypose.n_residue()).xyz("C"),
				mypose.residue(1).xyz("N"),
				mypose.residue(1).xyz("CA"),
				mypose.residue(1).xyz( (mypose.residue(1).has("CM")?"CM":"C") )
			)
		);
	}
	
	return;
}

//Function to store side-chain dihedral values:
/*void store_sidechains(
	const core::pose::Pose &mypose,
	utility::vector1 < utility::vector1 < core::Real > > &outvector
) {
	outvector.clear();
	for(core::Size ir=1; ir<=mypose.n_residue(); ir++) {
		utility::vector1 < core::Real > chivector;
		if(mypose.residue(ir).nchi()>0) {
			for(core::Size ichi=1; ichi<=mypose.residue(ir).nchi(); ichi++) chivector.push_back(mypose.chi(ichi, ir));
		}
		outvector.push_back(chivector); //Store the vector of chi values, even if empty.
	}
	return;
}*/

//A function to strip terminal glycines.  This returns "true" if there were glycines at both termini (that were stripped off), and false otherwise.
bool strip_terminal_glycines (core::pose::Pose &mypose) {
	if(mypose.residue(1).name1()=='G' && mypose.residue(mypose.n_residue()).name1() == 'G')
	{
		core::pose::Pose temp_pose;
		temp_pose.append_residue_by_jump(mypose.residue(2),1); //Put in the first residue.
		for(core::Size ir = 3; ir < mypose.n_residue(); ir++) temp_pose.append_residue_by_bond(mypose.residue(ir)); //Put in the rest.
		core::pose::add_lower_terminus_type_to_pose_residue(temp_pose, 1);
		core::pose::add_upper_terminus_type_to_pose_residue(temp_pose, temp_pose.n_residue());
		temp_pose.conformation().detect_disulfides(); //Detect the disulfides
		mypose = temp_pose;		
		return true;
	}
	return false;
}

//Function to set the backbone dihedral values of a pose to input values.  This also sets the side-chains of residues that aren't to be mutated.
void set_to_state(
	core::pose::Pose &mypose, //The pose to be altered (input/output)
	const utility::vector1 < core::Real > &bb_list //The list of backbone dihedral values, starting with psi1, omega1, phi2, psi2, omega2...
	//const utility::vector1 < utility::vector1 < core::Real > > &chi_list //The list of chi values, by amino acid.  Only used if v_preserve_cys is true.
	//const utility::vector1 < std::pair <core::Size, core::Size> > &disulf_positions //The list of disulfide-bonded positions
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//bool set_sc;

	for(core::Size ir=1, nres=mypose.n_residue(), icount=1; ir<=nres; ir++) { //Loop through all residues
		if(mypose.residue(ir).type().is_beta_aa()) { //beta-amino acids
			if(!mypose.residue(ir).is_lower_terminus()) betapeptide_setphi(mypose, ir, bb_list[icount++]); //phi
			betapeptide_settheta(mypose, ir, bb_list[icount++]); //theta
			if(!mypose.residue(ir).is_upper_terminus()) {
				betapeptide_setpsi(mypose, ir, bb_list[icount++]); //psi
				betapeptide_setomega(mypose, ir, bb_list[icount++]); //omega
			}
		} else if (mypose.residue(ir).type().is_alpha_aa()) { // alpha amino acids
			if(!mypose.residue(ir).is_lower_terminus()) mypose.set_phi(ir, bb_list[icount++]);
			if(!mypose.residue(ir).is_upper_terminus()) {
				mypose.set_psi(ir, bb_list[icount]);
				icount++;
				mypose.set_omega(ir, bb_list[icount]);
				icount++;
			} else {} // default case -- do nothing
		}
	}

	//if(disulf_positions.size()>0) { //If there are disulfide-bonded positions
		//mypose.conformation().fix_disulfides( disulf_positions );
	//}

	if(option[v_cyclic]()) {
		core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
		core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.n_residue());
		const std::string hstring = (mypose.residue(1).name1()=='P' ? "CD" : "H");
		mypose.conformation().set_torsion_angle(//Set psi of last residue
			core::id::AtomID(mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CA":"N") ), mypose.n_residue()),
			core::id::AtomID(mypose.residue(mypose.n_residue()).atom_index( (mypose.residue(mypose.n_residue()).has("CM")?"CM":"CA") ), mypose.n_residue()),
			core::id::AtomID(mypose.residue(mypose.n_residue()).atom_index("C"), mypose.n_residue()),
			core::id::AtomID(mypose.residue(mypose.n_residue()).atom_index("O"), mypose.n_residue()),
			(bb_list[3*mypose.n_residue()-2]-180.0)/180.0*PI);
		mypose.conformation().set_torsion_angle( //Set phi of first residue
			core::id::AtomID(mypose.residue(1).atom_index(hstring), 1),
			core::id::AtomID(mypose.residue(1).atom_index("N"), 1),
			core::id::AtomID(mypose.residue(1).atom_index("CA"), 1),
			core::id::AtomID(mypose.residue(1).atom_index( (mypose.residue(1).has("CM")?"CM":"C") ), 1),
			(bb_list[3*mypose.n_residue()]-180.0)/180.0*PI);
	} //else {
		//core::pose::add_lower_terminus_type_to_pose_residue(mypose, 1);
		//core::pose::add_upper_terminus_type_to_pose_residue(mypose, mypose.n_residue());
	//}

	mypose.update_residue_neighbors();
	//mypose.conformation().detect_disulfides();

	return;
}

//Function to mutate a pose to a given sequence:
void mutate_to_sequence (
	const std::string &sequence,
	core::pose::Pose &mypose,
	const utility::vector1 <bool> &dpositions,
	const utility::vector1 <bool> &betapositions,
	bool const force_mutation
	//const utility::vector1 < std::pair < core::Size, core::Size > > &disulf_positions
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(sequence.length()!=mypose.n_residue()) {
		printf("Error!  Sequence length/structure size mismatch.  Crashing.\n"); fflush(stdout);
		exit(1);
	}

	if(option[v_cyclic]()) mypose.remove_constraints(); //Remove all constraints

	for(core::Size ir=1; ir<=mypose.n_residue(); ir++) {
		bool continue_on = false;
		if(!force_mutation) {
			if(option[v_ignoreresidue]().size()>0) {
				//Don't change residues in the ignore list.
				for(core::Size i=1; i<=option[v_ignoreresidue]().size(); i++) {
					if((core::Size)option[v_ignoreresidue][i]==ir) {
						continue_on=true;
						break; //out of the for loop
					}
				}
			}
			if(!continue_on && option[v_norepack_positions]().size()>0) {
				//Don't mutate residues that can't be repacked.
				for(core::Size i=1, imax=option[v_norepack_positions]().size(); i<=imax; ++i) {
					if((core::Size)option[v_norepack_positions]()[i]==ir) {
						continue_on=true;
						break; //out of the for loop
					}
				}
			}
		} else {
			if(!mypose.residue(ir).type().is_polymer()) continue_on=true; //Don't mutate non-polymer positions.
			if(!mypose.residue(ir).type().is_alpha_aa() && !mypose.residue(ir).type().is_beta_aa()) continue_on=true; //Don't mutate positions that aren't alpha or beta amino acids, for now.
			if(mypose.residue(ir).name3()=="CYX" || mypose.residue(ir).name3()=="ASX" || mypose.residue(ir).name3()=="LYX" || mypose.residue(ir).name3()=="ORX" || mypose.residue(ir).name3()=="GLX") continue_on=true; //Ugly hack.  Fix this with something more general!
		}
		if(continue_on) continue;

		if(option[v_preserve_chirality]() && dpositions[ir]) {
			std::string aaname;
			get_D_name(sequence[ir-1], aaname);
			protocols::simple_moves::MutateResidue mutres(ir, aaname);
			mutres.apply(mypose);
		} else if (betapositions[ir]) { //beta-3-amino acid (note: D-beta-3-amino acids are not currently supported)
			//printf("about to mutate beta\n"); fflush(stdout); //DELETE ME
			std::string aaname;
			get_beta_name(sequence[ir-1], aaname);
			//printf("betaname=%s\n", aaname.c_str()); fflush(stdout); //DELETE ME
			protocols::simple_moves::MutateResidue mutres(ir, aaname);
			mutres.apply(mypose);
			//printf("mutated beta\n"); fflush(stdout); //DELETE ME
		} else {
			protocols::simple_moves::MutateResidue mutres(ir, sequence[ir-1]);
			mutres.apply(mypose);
		}
	}

	if(option[v_cyclic]()) {
		core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
		core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.n_residue());
		addcyclicconstraints(mypose); //Add back the constraints (need to do this in case of !P->P or P->!P mutations at N-terminus.
	} //else {
		//core::pose::add_lower_terminus_type_to_pose_residue(mypose, 1);
		//core::pose::add_upper_terminus_type_to_pose_residue(mypose, mypose.n_residue());
	//}

	//if(disulf_positions.size()>0) { //If there are disulfide-bonded positions
		//mypose.conformation().fix_disulfides( disulf_positions );
	//}

	mypose.update_residue_neighbors();

	return;
}

//Function to transmit additional information
void transmit_additional_info (
	const int fromproc,
	const int toproc,
	const int thisproc,
	const int numelements,
	utility::vector1 < core::Real > &datavector //input or output
) {
	MPI_Status mpistatus;
	double *dataarray = new double [numelements];
	if(thisproc==fromproc) { //If this is the sending proc
		for (core::Size i=1; i<=(core::Size)numelements; i++) dataarray[i-1]=(double)datavector[i]; //Store the data in the data array for transmission
		MPI_Send(dataarray, numelements, MPI_DOUBLE, toproc, 0, MPI_COMM_WORLD); //Send the data
	} else if (thisproc==toproc) { //If this is the receiving proc
		MPI_Recv(dataarray, numelements, MPI_DOUBLE, fromproc, 0, MPI_COMM_WORLD, &mpistatus); //Receive the data
		datavector.resize(numelements); //Resize the vector to store the data
		for (core::Size i=1; i<=(core::Size)numelements; i++) datavector[i]=(core::Real)dataarray[i-1]; //Store the data in the data vector 
	}
	delete [] dataarray;

	//Do nothing if this is not fromproc or toproc.
	return;
}


//Function to receive a silent structure as a string from another proc (sent using transmit_silent_struct()) and to convert it to a pose:
void receive_silent_struct (
	core::pose::Pose &pose, //receiving pose
	const int fromproc,
	const int toproc,
	const int thisproc,
	bool const broadcast //Is this from a broadcast?
) {
	if(!broadcast && (thisproc!=toproc)) return; //Do nothing if this isn't the receiving proc.
	unsigned long charlength=0;
	//unsigned long taglength=0;

	MPI_Status mpistatus;

	utility::vector1 < std::string > tagvector;
	
	if(broadcast) {
		MPI_Bcast( &charlength, 1, MPI_UNSIGNED_LONG, fromproc, MPI_COMM_WORLD); //Receive the length of the char buffer sent by MPI_Bcast
	} else {
		MPI_Recv( &charlength, 1, MPI_UNSIGNED_LONG, fromproc, 0, MPI_COMM_WORLD, &mpistatus); //Receive the length of the char buffer
	}
	char* received_char = new char[charlength];

	if(broadcast) {
		MPI_Bcast( received_char, charlength, MPI_CHAR, fromproc, MPI_COMM_WORLD); //Recieve the char sent by MPI_Bcast
	} else {
		MPI_Recv( received_char, charlength, MPI_CHAR, fromproc, 0, MPI_COMM_WORLD, &mpistatus); //Recieve the char
	}
	std::string received_string( received_char, charlength );
	std::istringstream received_istringstream( received_string );

	//printf("Proc %i received string %s\n.", thisproc, received_string.c_str()); fflush(stdout); //DELETE ME.

	core::io::silent::SilentFileData silentfiledata;
	silentfiledata.read_stream( received_istringstream, tagvector, true, "suppress_bitflip" );
	pose.clear();
	silentfiledata[silentfiledata.tags()[1]]->fill_pose(pose);

	//printf("Proc %i filled pose.\n", thisproc); fflush(stdout); //DELETE ME.

	delete [] received_char;

	add_user_constraints(pose);

	return;
}


//Function to convert a silent structure to a string and transmit it to another proc
void transmit_silent_struct (
	const int fromproc,
	const int toproc, //Not used if broadcast=true
	const int thisproc,
	core::io::silent::SilentStructOP &silentstruct,
	bool const broadcast //Should this be transmitted to a single other proc, or broadcast to all?
) {
	if(thisproc!=fromproc) return; //Do nothing if this isn't the transmitting proc.

	core::io::silent::SilentFileData silentfiledata;
	std::stringbuf sb;
	std::ostream outstream(&sb);

	//Add the SEQUENCE line, which for some reason isn't written to the stream by default except on file output:
	std::stringstream header;
	silentstruct->print_header( header );
	outstream << header.str();

	silentfiledata.write_silent_struct( (*silentstruct), outstream );
	std::string outstring = sb.str();
	//printf("transmit_silent_struct from proc %i will send %s\n", thisproc, outstring.c_str()); fflush(stdout); //DELETE ME
	unsigned long strlength = outstring.length()+1;
	char *outchar = new char [strlength];
	sprintf(outchar, "%s", outstring.c_str());

	if(broadcast) {
		MPI_Bcast(&strlength, 1, MPI_UNSIGNED_LONG, fromproc, MPI_COMM_WORLD); //Send the length of the string that I'm about to send.
		MPI_Bcast(outchar, strlength, MPI_CHAR, fromproc, MPI_COMM_WORLD); //Send the data
	} else { //if not broadcast (specific receiver)
		MPI_Send(&strlength, 1, MPI_UNSIGNED_LONG, toproc, 0, MPI_COMM_WORLD); //Send the length of the string that I'm about to send.
		MPI_Send(outchar, strlength, MPI_CHAR, toproc, 0, MPI_COMM_WORLD); //Send the data
	}

	delete [] outchar;
	return;
}


//Overloaded function to convert a pose to a silent structure, then to convert a silent structure to a string and transmit it to another proc
void transmit_silent_struct (
	const int fromproc,
	const int toproc,
	const int thisproc,
	std::string const &tag,
	core::pose::Pose const &pose,
	bool const broadcast
) {
	if(thisproc!=fromproc) return; //Do nothing if this isn't the transmitting proc.
	core::io::silent::SilentStructOP silentstruct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
	silentstruct->fill_struct(pose, tag);
	transmit_silent_struct( fromproc, toproc, thisproc, silentstruct, broadcast);
	return;
}


//Function to divide the calculation of states' energies over many processors, and to calculate energies of each state.  This function also transmits the current sequence from the master proc to all slaves.
void scoreall (
	const int procnum, //The process number
	const int proccount, //The total number of procs
	const core::Size totalstatecount, //The total number of states
	const core::Size PCAvectorcomponentcount, //The number of entries in each PCA vector in the PCA matrices
	std::string &sequence, //The current sequence
	utility::vector1<bool> &dpositions, //The list of D-amino acid positions
	utility::vector1<bool> &betapositions, //The list of beta-amino acid positions
	utility::vector1 < core::io::silent::SilentStructOP  > &allstates_master, //Vector of silent structures -- MASTER
	utility::vector1 < std::pair < core::Size, core::Size > > &disulf_positions, //The list of pairs of disulfide-bonded cysteines.
	//utility::vector1 < utility::vector1 < utility::vector1 < core::Real > > > &allstates_master_sidechains, //Vector of chi values (for each amino acid, for each state) -- MASTER, only used if v_preserve_cys is true
	utility::vector1 < core::Real > &energyvals_current, //Vector of energy values -- output
	utility::vector1 < core::Real > &rmsvals_current, //Vector of rms values -- output
	const core::pose::Pose &positive_master, //The positive master state -- MASTER (input)
	core::pose::Pose &positive_current, //For proc 0, this will be the energy-minimized positive state -- output.  For other procs, it's a temporary buffer.
	numeric::random::RandomGenerator &RG, //A random number generator
	protocols::relax::FastRelax &frlx, //A fastrelax mover
	core::scoring::ScoreFunctionOP const &sfxn, //The owning pointer for the scorefunction
	utility::vector1 < utility::vector1 <core::Real> > &PCAvariances, //The list of variance vectors (for each state)
	utility::vector1 < utility::vector1 < utility::vector1 < core::Real > > > &PCAmatrices, //The list of PCA matrices (for each state)
	utility::vector1 <core::id::NamedAtomID> const &extra_rms_atoms //List of additional atoms to use in RMSD calculation
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1 < core::Size > procassignments; //Which process does which state?  The nth entry is for the nth state, and stores the associated process number.
	int mpibuffer;
	MPI_Status mpistat;

	//Transmit/receive the current sequence (note that the length of "sequence" passed to all procs must be the same!):
	char charbuffer [256];
	if(procnum==0) strcpy (charbuffer, sequence.c_str());
	MPI_Bcast(&charbuffer, sequence.length()+1, MPI_CHAR, 0, MPI_COMM_WORLD);
	sequence=charbuffer;
	//if(procnum!=0) { printf("\tProc %i received sequence %s.\n", procnum, sequence.c_str()); fflush(stdout); } //DELETE ME.

	//Variables needed to split replicates over multiple processess:
	const core::Size numstates = option[v_split_replicates_over_processes]() ? option[v_MCminimize_replicates]() * totalstatecount : totalstatecount; //Number of states is multiplied by number of replicates iff the split over procs option is used.
	//utility::vector1 < core::pose::Pose > positivepose_list;
	core::Size lowestE_positivepose = 1;
	core::Real lowestE_positivepose_energy = 0.0;

	/********** PROC 0 SENDS AND RECEIVES JOBS **********/
	if(procnum==0) {
		for (core::Size statecount=1; statecount<=numstates+(core::Size)proccount-1; statecount++) {

			//Receive requests from each process for states to compute, and keep track of which process does which state
			MPI_Recv(&mpibuffer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistat); //A slave process will send an int as a query.
			if(statecount>numstates) { //If we've exhausted the available states
				//printf("Sending 0 to proc %i.\n", mpistat.MPI_SOURCE); fflush(stdout); //DELETE ME
				mpibuffer=0; //Send "0" as a message to indicate that all states are exhausted.
				MPI_Send(&mpibuffer, 1, MPI_INT, mpistat.MPI_SOURCE, 0, MPI_COMM_WORLD); //Send "0" to the slave that sent the query.
			} else { //If there are still states to send out to slave processes
				mpibuffer=statecount;
				MPI_Send(&mpibuffer, 1, MPI_INT, mpistat.MPI_SOURCE, 0, MPI_COMM_WORLD); //Send the number of the next uncalculated state to the slave that sent the query.

				if(option[v_savememory]()) {	//If procs aren't storing the full state list:
					const core::Size thisstate = option[v_split_replicates_over_processes]() ? (core::Size)(div( (int)statecount-1, (int) option[v_MCminimize_replicates]() ).quot) + 1 : statecount;
					//printf("Assigning state %i to proc %i.\n", (int)thisstate, (int)mpistat.MPI_SOURCE); fflush(stdout); //DELETE ME
					//mpibuffer=allstates_master[thisstate].size();

					//MPI_Send(&mpibuffer, 1, MPI_INT, mpistat.MPI_SOURCE, 0, MPI_COMM_WORLD); //Send the size of the allstates_master vector.					

					//transmit_additional_info(0, mpistat.MPI_SOURCE, 0, mpibuffer, allstates_master[thisstate]); //Transmit the allstates_master vector (mpibuffer is still the size of the PCAvariances vector)

					transmit_silent_struct( 0, mpistat.MPI_SOURCE, 0, allstates_master[thisstate], false );

					if(option[v_PCAfiles].user() || option[v_PCAfiles_list].user()) { //If there are PCAfiles:
						mpibuffer=PCAvariances[thisstate].size();

						MPI_Send(&mpibuffer, 1, MPI_INT, mpistat.MPI_SOURCE, 0, MPI_COMM_WORLD); //Send the number of PCA vectors
						if(mpibuffer>0) { //(mpibuffer is still the size of the PCAvariances vector).
							transmit_additional_info(0, mpistat.MPI_SOURCE, 0, mpibuffer, PCAvariances[thisstate]); //Transmit the PCA variances (mpibuffer is still the size of the PCAvariances vector).
							const core::Size variancescount = (core::Size)mpibuffer;
							for(core::Size ii=1; ii<=variancescount; ii++) { //variancescount is the number of PCA vectors (or size of the PCAvariances vector)
								transmit_additional_info(0, mpistat.MPI_SOURCE, 0, PCAvectorcomponentcount, PCAmatrices[thisstate][ii]); //Transmit each PCA vector (mpibuffer is still the size of the PCAvariances vector).
							}
						}

					}
				}
				procassignments.push_back(mpistat.MPI_SOURCE); //Store the number of the slave.
			} //end if
		}//end for

		//printf("Proc %i has reached the barrier.\n", procnum); fflush(stdout); //DELETE ME
		MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side

		printf("Collecting computed energies from slave processes...\n"); fflush(stdout);

		if(option[v_consider_all_replicates]()) {
			//If we're scoring ALL replicates (not just the lowest-energy one for each state), then
			//these two vectors need to store values for all replicates, in the order state1-rep1,
			//state1-rep2, state1-rep3... state1-repN, state2-rep1, state2-rep2... etc.
			energyvals_current.resize(numstates);
			rmsvals_current.resize(numstates);
		} else {
			//If we're scoring just the lowest-energy replicate for each state, then these two vectors
			//just store one value per state, in order of state (state1, state2, state3... etc.)
			energyvals_current.resize(totalstatecount);
			rmsvals_current.resize(totalstatecount);
		}

		//********************************************************************************
		//********************************************************************************
		//TODO: REWRITE THIS SECTION TO ALLOW USE OF THE v_consider_all_replicates OPTION!
		//********************************************************************************
		//********************************************************************************

		//These variables are necessary if minimizations are split over multiple processes:
		core::Real curstate_curenergy, curstate_currms;
		core::Size curstate = 0;

		if(option[v_split_replicates_over_processes]) { //If we ARE dividing replicates over different processes:
			for(core::Size statecount=1; statecount<=numstates; statecount++) {
				mpibuffer=1;
				MPI_Send(&mpibuffer, 1, MPI_INT, procassignments[statecount], 0, MPI_COMM_WORLD); //Send an inquiry to the proc
				MPI_Recv(&curstate_curenergy, 1, MPI_DOUBLE, procassignments[statecount], 0, MPI_COMM_WORLD, &mpistat); //Recieve the energy from the proc
				MPI_Recv(&curstate_currms, 1, MPI_DOUBLE, procassignments[statecount], 0, MPI_COMM_WORLD, &mpistat); //Recieve the rms from the proc
				if(option[v_consider_all_replicates]()) {
					if((statecount-1) % option[v_MCminimize_replicates]() == 0) curstate++; //If this is the first replicate for this state, increment the current state
					energyvals_current[statecount] = curstate_curenergy; //Store the energy value as the state's energy value
					rmsvals_current[statecount] = curstate_currms; //Store the rms as the state's rms
					if(curstate==1 && (statecount==1 || curstate_curenergy < lowestE_positivepose_energy)) { //Check whether this is the lowest-energy POSITIVE state encountered so far, and if so, store its index.
						lowestE_positivepose_energy=curstate_curenergy;
						lowestE_positivepose=statecount;
					}
				} else {
					if((statecount - 1) % option[v_MCminimize_replicates]() == 0) { //If this is the first replicate for this state
						curstate++; //Increment the number of the current state
						energyvals_current[curstate] = curstate_curenergy; //Store the energy value as the state's energy value
						rmsvals_current[curstate] = curstate_currms; //Store the rms as the state's rms
						//printf("ping%i\n", curstate); fflush(stdout); //DELETE ME
					} else { //If this is not the first replicate for this state
						if(curstate_curenergy < energyvals_current[curstate]) { //Check whether this replicate returned a lower energy.  If so, overwrite the energy and rms for this state.
							if(curstate==1) lowestE_positivepose = statecount; //If this is the first state, store the index of the lowest energy replicate
							energyvals_current[curstate] = curstate_curenergy;
							rmsvals_current[curstate] = curstate_currms;
						}
					}
				}
			}
			//positive_current = positivepose_list[lowestE_positivepose]; //Store the lowest-energy pose for the positive state (of all the replicates).
		} else { //if we're NOT splitting replicates over different processes, preserve old behaviour: 
			for(core::Size statecount=1; statecount<=procassignments.size(); statecount++) {
				mpibuffer=1;
				MPI_Send(&mpibuffer, 1, MPI_INT, procassignments[statecount], 0, MPI_COMM_WORLD); //Send an inquiry to the proc
				MPI_Recv(&energyvals_current[statecount], 1, MPI_DOUBLE, procassignments[statecount], 0, MPI_COMM_WORLD, &mpistat); //Recieve the energy from the proc
				MPI_Recv(&rmsvals_current[statecount], 1, MPI_DOUBLE, procassignments[statecount], 0, MPI_COMM_WORLD, &mpistat); //Recieve the rms from the proc
			}
		}

		//Figure out which positive pose was lowest-energy.  Send that proc number out, and have that proc return its energy-minimized positive pose structure.
		if(option[v_split_replicates_over_processes]()) mpibuffer = (int)procassignments[lowestE_positivepose]; //If we're splitting replicates over procs, I want the proc that did the lowest-energy replicate of the positive state.
		else mpibuffer = (int)procassignments[1]; //If we're not splitting replicates over procs, I want the proc that did the positive state.
		MPI_Bcast(&mpibuffer, 1, MPI_INT, 0, MPI_COMM_WORLD); //mpibuffer is the number of the proc from which the master wants a structure.
		positive_current.clear(); //Clear the pose.
		receive_silent_struct( positive_current, mpibuffer, 0, 0, false );
	} //if procnum==0
	/********** ALL OTHER PROCS DO THE FOLLOWING **********/
	else {
		mpibuffer=1;
		energyvals_current.clear(); //Used for temporary storage.
		rmsvals_current.clear(); //Used for temporary storage.
		core::Size this_proc_did_positive=false; //Did this proc do the positive state?
		core::Real positive_lowestEvalue=0.0; //The energy of the lowest-energy replicate of the positive state done by this proc.

		//If we're conserving memory, we need places to put the backbone vector and PCA vectors for this state when they're transmitted from the master proc:
		utility::vector1 <core::Real> thisstate_phipsivector; //A container for this state's phi and psi information, if it needs to be transmitted to this proc.
		utility::vector1 <core::Real> thisstate_PCAvariances; //A container for this state's PCA variances, if they need to be transmitted to this proc.
		utility::vector1 < utility::vector1 < core::Real > > thisstate_PCAmatrix; //A container for this state's PCA matrix, if it needs to be transmitted to this proc.

		//Pointers to the current data (might give const-access problems):
		//utility::vector1 < core::Real > * curphipsivector; 
		utility::vector1 < core::Real > * curPCAvariances;
		utility::vector1 < utility::vector1 < core::Real > > * curPCAmatrix;

		do { //will keep doing until mpibuffer == 0
			MPI_Send(&mpibuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); //Send an int as a query (the value doesn't matter).
			MPI_Recv(&mpibuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistat); //Recieve the next state to calculate from the master process.
			if(mpibuffer==0) break; //Exit the do loop if the master process sends "0".

			//Determine which state we're actually minimizing:
			const core::Size thisstate = option[v_split_replicates_over_processes]() ? (core::Size)(div( (int)mpibuffer-1, (int) option[v_MCminimize_replicates]() ).quot) + 1 : mpibuffer;

			//A temporary pose to receive the state:
			core::pose::Pose temppose;
			temppose.clear();

			//If we don't have the data for this state in this proc, yet, we need to receive it from the master proc:
			if(option[v_savememory]()) {
				receive_silent_struct(temppose, 0, procnum, procnum, false);
				if(option[v_PCAfiles].user() || option[v_PCAfiles_list].user()) { //If there are PCAfiles:
					MPI_Recv(&mpibuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistat); //Receive the size of the PCA variances vector
					if(mpibuffer>0) { //(mpibuffer is still the size of the PCAvariances vector).
						transmit_additional_info(0, procnum, procnum, mpibuffer, thisstate_PCAvariances); //Recive the PCA variances (mpibuffer is the size of the PCA variances vector).
						for(core::Size ii=1; ii<=(core::Size)mpibuffer; ii++) { //mpibuffer is still the number of PCA vectors (the size of the variances vector)
							utility::vector1 < core::Real > PCAvector; //A container for the PCA vector
							transmit_additional_info(0, procnum, procnum, PCAvectorcomponentcount, PCAvector); //Transmit each PCA vector (mpibuffer is still the size of the PCAvariances vector).
							thisstate_PCAmatrix.push_back(PCAvector); //Store the vector in the PCA matrix
						}
					}
				}
			}

			//Set the pointers:
			if(!option[v_savememory]()) { 
				//curphipsivector = &(allstates_master[thisstate]); 
				curPCAvariances = &(PCAvariances[thisstate]);
				curPCAmatrix = &(PCAmatrices[thisstate]);
			} else {
				//curphipsivector = &thisstate_phipsivector;
				curPCAvariances = &thisstate_PCAvariances;
				curPCAmatrix = &thisstate_PCAmatrix;
			}

			mutate_to_sequence(sequence, temppose, dpositions, betapositions, false); //Mutate the copy of the state in this proc to the current sequence

			//Set the conformation:
			//set_to_state(temppose, (*curphipsivector) /*, disulf_positions*/);

			//Energy-minimize and score the negative state:
			if(option[v_PCAfiles].user() || option[v_PCAfiles_list].user()) {
				//printf("%i Perturb_and_relax\n", procnum); fflush(stdout); //DELETE ME!
				//printf("%i pcavariances[%i][1]=%.4f\n", procnum, thisstate, PCAvariances[thisstate][1]); fflush(stdout); //DELETE ME!
				perturb_bb_and_relax(option[v_split_replicates_over_processes]() ? 1 : option[v_MCminimize_replicates](), option[v_MCminimize_steps](), temppose, disulf_positions, RG, option[v_bb_perturbation], frlx, sfxn, option[v_kbt_minimize](), (*curPCAvariances), (*curPCAmatrix) );
			} else {
				//If there are no PCA vectors provided by the user, I need to pass some dummy variables.
				utility::vector1 <core::Real> PCAvariances_dummy;
				utility::vector1 <utility::vector1 <core::Real> > PCAmatrices_dummy;
				perturb_bb_and_relax(option[v_split_replicates_over_processes]() ? 1 : option[v_MCminimize_replicates](), option[v_MCminimize_steps](), temppose, disulf_positions, RG, option[v_bb_perturbation], frlx, sfxn, option[v_kbt_minimize](), PCAvariances_dummy, PCAmatrices_dummy );
			}

			//Debug output -- DELETE ME!:
			//{
			//	char tempfilename[128];
			//	sprintf(tempfilename, "PROC_%02i_struct_%03lu.pdb", procnum, thisstate);
			//	temppose.dump_scored_pdb( tempfilename, *sfxn );
			//}
	
			//Determine whether the proc did the positive pose:
			if(thisstate==1) {
				if(!this_proc_did_positive) { //Case 1: this is the first time the proc did a replicate of the positive state
					this_proc_did_positive = true; //Set this to true.
					positive_lowestEvalue = temppose.energies().total_energy(); //Store the energy
					positive_current = temppose; //Store the pose.
				} else { //Case 2: the proc has done other replicates of the positive state
					if(temppose.energies().total_energy() < positive_lowestEvalue) { //If the energy is the lowest of any of the positive pose replicates...
						positive_lowestEvalue = temppose.energies().total_energy(); //Store the energy
						positive_current = temppose; //Store the pose.
					} //...otherwise do nothing.
				}
			}

			energyvals_current.push_back(temppose.energies().total_energy());
			rmsvals_current.push_back( get_distance_measure(temppose, positive_master, extra_rms_atoms)  ); //Calculate the rms to target, and store it.

			//printf("Proc %i is finished with state %i.\n", (int)procnum, (int)thisstate); fflush(stdout); //DELETE ME!
		} while(true); //Keep looping; only break out when the master proc transmits "0".

		//printf("Proc %i has reached the barrier.\n", procnum); fflush(stdout); //DELETE ME
		MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side

		if(energyvals_current.size() > 0) {
			for(core::Size ii=1; ii<=energyvals_current.size(); ii++) {
				MPI_Recv(&mpibuffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistat); //Receive an inquiry from the master proc
				MPI_Send(&energyvals_current[ii], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //Send the energy value to the master proc
				MPI_Send(&rmsvals_current[ii], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //Send the energy value to the master proc
			} //Do this once for each energy value that this proc calculated.
		}

		//Receive a proc num from the master (broadcast to all); if this is that proc, transmit the positive pose
		MPI_Bcast(&mpibuffer, 1, MPI_INT, 0, MPI_COMM_WORLD); //mpibuffer is the number of the proc from which the master wants a structure.
		if(procnum==mpibuffer) { //If this is the proc named by the master, transmit the structure.
			//Since the master wants the very lowest-energy replicate of the positive state, it figures out which proc did the lowest-energy replicate, and asks for its lowest-energy replicate.
			transmit_silent_struct( mpibuffer, 0, procnum, "lowest_energy_positive", positive_current, false );
		}

	} //if procnum!=0

	return;
} //scoreall

/********************************************************************************
	Function to set up covalent connections in a pose based on a list of
	pairs of NamedAtomIDs.
********************************************************************************/
void set_up_covalent_connections (
	core::pose::Pose &pose,
	utility::vector1 < std::pair < core::id::NamedAtomID, core::id::NamedAtomID > > const &conlist
) {
	core::Size const conlist_size = conlist.size();
	if(conlist_size==0) return;

	for(core::Size i=1; i<=conlist_size; ++i) {
		runtime_assert_string_msg(
			pose.residue(conlist[i].first.rsd()).has(conlist[i].first.atom()) &&
			pose.residue(conlist[i].second.rsd()).has(conlist[i].second.atom()),
			"In set_up_covalent_connections: pose is missing one or more atoms specified in the covalent connections file."
		);
		pose.conformation().declare_chemical_bond(
			conlist[i].first.rsd(),
			conlist[i].first.atom(),
			conlist[i].second.rsd(),
			conlist[i].second.atom()
		);
		/*printf("Added covalent connection between residue %lu, atom %s and residue %lu, atom %s.\n",
			conlist[i].first.rsd(),
			conlist[i].first.atom().c_str(),
			conlist[i].second.rsd(),
			conlist[i].second.atom().c_str()
		);*/
	}

	return;
}

/********************************************************************************
	Function to read in a file that lists covalent connections and stores
	these in a list of pairs of NamedAtomIDs.
********************************************************************************/
void read_covalent_connections_file (
	utility::vector1 < std::pair < core::id::NamedAtomID, core::id::NamedAtomID > > &conlist,
	core::Size const procnum
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	conlist.clear();

	//Open the list file for reading
	std::ifstream lfile;
	std::string connectionfile = option[v_covalent_connections]();
	printf("Reading connection file %s\n", connectionfile.c_str()); fflush(stdout);
	lfile.open(connectionfile.c_str(), std::ifstream::in);
	std::string line;
	while (std::getline(lfile, line)) //Read each line
	{
		utility::vector1 < core::Size > curline;
		char tempchar [512]; sprintf(tempchar, "%s", line.c_str()); //Convert string->char
		char * curcharpointer = &tempchar[0];

		char * token1 = strtok( curcharpointer, " " );
		char * token2 = strtok( NULL, " " );
		char * token3 = strtok( NULL, " " );
		char * token4 = strtok( NULL, " " );

		core::Size res1 = static_cast<core::Size>( strtol(token1, &token1, 10) );
		core::Size res2 = static_cast<core::Size>( strtol(token3, &token3, 10) );	
		std::string atom1(token2);
		std::string atom2(token4);
		core::id::NamedAtomID atid1(atom1, res1);
		core::id::NamedAtomID atid2(atom2, res2);
		conlist.push_back( std::pair <core::id::NamedAtomID,core::id::NamedAtomID>(atid1, atid2) );
	}

	//Close the file:
	lfile.close();	

	//Output a summary of the file that we just read:
	if(procnum==0) {
		printf("The following covalent connections will be added:\n");
		printf("Res1\tAtom1\tRes2\tAtom2\n");
		for (core::Size i=1, imax=conlist.size(); i<=imax; ++i) {
			printf("%lu\t%s\t%lu\t%s\n", conlist[i].first.rsd(), conlist[i].first.atom().c_str(), conlist[i].second.rsd(), conlist[i].second.atom().c_str());
		}
		printf("\n"); fflush(stdout);
	}
	
	return;
}


/********************************************************************************
	Function to read in a file that lists positions that are equivalent (i.e.
	positions that must be mutated to the same residue when mutations are
	introduced).  This file consists of several lines, each of which has a
	list of numbers.  All numbered positions in a line are equivalent.
********************************************************************************/
void read_equivalent_position_file (
	const char* filename,
	utility::vector1 < utility::vector1 < core::Size > > &equivalent_positions
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	equivalent_positions.clear();

	//Open the list file for reading
	std::ifstream lfile;
	lfile.open(filename, std::ifstream::in);

	std::string line;
	while (std::getline(lfile, line)) //Read each line
	{
		utility::vector1 < core::Size > curline;
		char tempchar [512]; sprintf(tempchar, "%s", line.c_str()); //Convert string->char
		char * curcharpointer = &tempchar[0];

		//Pull integers out of the string:
		while(true) {
			core::Size curint = (core::Size)strtol(curcharpointer, &curcharpointer, 10);
			if(curint == 0) break;
			curline.push_back(curint);
			if((*curcharpointer)=='\0') break;
		}

		equivalent_positions.push_back(curline); //Store the current set of integers.
	}

	//Close the file:
	lfile.close();

	return;

}

/********************************************************************************
	Function to read in a PCA file and output a variances vector and a
	vector of principal component vectors.
********************************************************************************/
core::Size read_PCAfile (
	const char* filename,
	utility::vector1 <core::Real> &variances, //Output -- variances data stored here
	utility::vector1 < utility::vector1 < core::Real > > &PCA_vectors, //Output -- PCA vector data stored here
	bool normalize_variances //If true, the largest variance is set to 1, and all others are rescaled in proportion
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Clear the storage boxes:
	variances.clear();
	PCA_vectors.clear();

	//Open the PCA file for reading
	FILE* pcafile = fopen(filename, "r");

	//Get the size of each vector, and the number of vectors
	core::Size realvectorcount, vectorcount, vectorcomponentcount;
	int readresult = fscanf(pcafile, "%lu %lu", &vectorcomponentcount, &realvectorcount);
	runtime_assert_string_msg(readresult, "Error in read_PCAfile(): Read failed on PCA file initial line.");
	vectorcount=realvectorcount;

	//If there are no principal component vectors for this state, return without doing anything.
	//The sizes of the variances and PCA_vectors vectors will be 0, in this case.
	if(vectorcount==0) return 0;

	//If we're limiting the number of PCA vectors, do so here:
	if(option[v_PCAlimit].user() && (core::Size)option[v_PCAlimit]() > vectorcount) vectorcount = (core::Size)option[v_PCAlimit]();

	//Resize the variances and PCA_vectors vectors:
	variances.resize(vectorcount);
	PCA_vectors.resize(vectorcount);

	//Loop through and read the variances vector:
	core::Real realbuffer;
	for(core::Size i=1; i<=realvectorcount; i++) {
		readresult=fscanf(pcafile, "%lf", &realbuffer);
		runtime_assert_string_msg(readresult, "Error in read_PCAfile(): Read failed on PCA file variances vector line.");		
		if(i<=vectorcount) variances[i]=realbuffer;
	}

	//Loop through and populate the PCA matrix:
	for(core::Size i=1; i<=vectorcount; i++) {
		PCA_vectors[i].resize(vectorcomponentcount);
		for(core::Size j=1; j<=vectorcomponentcount; j++) {
			readresult=fscanf(pcafile, "%lf", &realbuffer);
			runtime_assert_string_msg(readresult, "Error in read_PCAfile(): PCA file read failed while reading the PCA matrix.");		
			PCA_vectors[i][j]=realbuffer;
		}
	}

	//Close the file:
	fclose(pcafile);

	//Normalize the variances vector to the size of the first entry:
	if(normalize_variances) {for(core::Size i=variances.size(); i>=1; i--) variances[i]=variances[i]/variances[1];}

	return vectorcomponentcount; //Return the number of entries in each PCA vector.
}

/********************************************************************************
	Function to parse a text file containing a list of files, and to put that
	list of files into an array of strings.
********************************************************************************/
void parsefilelist (
	const std::string &listfile, //Input -- name of file that contains list of files
	utility::vector1 < std::string > &filevector //Output -- list of files
) {
	//Clear the file vector
	filevector.clear();

	//Open the list file for reading
	std::ifstream lfile;
	lfile.open(listfile.c_str(), std::ifstream::in);

	std::string line;
	while (std::getline(lfile, line)) //Read each line
	{
		filevector.push_back(line); //Store each line
	}

	//Close the file:
	lfile.close();
}

//MAIN
int main(int argc, char *argv[]) {
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace core::id;

	numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

	//Initialize MPI:
	int procnum, totalprocs; //The current process number and the total number of processes.
	MPI_Init(&argc, &argv); //Initialize
	MPI_Comm_rank(MPI_COMM_WORLD, &procnum); //Get the current process number.
	MPI_Comm_size(MPI_COMM_WORLD, &totalprocs); //Get the total number of processes.
	char currentseq [256];

	//The following few lines are only for debugging with gdb.  Comment these out otherwise.
	/*if(procnum==1) {
		 int i = 0;
		 char hostname[256];
		 gethostname(hostname, sizeof(hostname));
		 printf("PID %d on %s ready for attach\n", getpid(), hostname);
		 fflush(stdout);
		 while (0 == i)
		     sleep(5);
	}*/

	
	if(procnum == 0){
		printf ("Starting design_cycpeptide_MPI.cc\n");
		printf ("Pilot app created 14 September 2012 by Vikram K. Mulligan, Baker Laboratory.\n");
		printf ("Starting %i parallel processes.\n", totalprocs);
		fflush(stdout);
	}
	register_options();
	devel::init(argc, argv); //Initialize Rosetta
	core::scoring::ScoreFunctionOP sfxn;
	if(procnum!=0) {
		sfxn = core::scoring::get_score_function(); //NOTE -- proc 0 never creates the scorefunction.  It SHOUOLDN'T ever use it, as a memory-saving thing.
		if(option[v_cst_file].user()) { //If a constraints file has been specified by the user, turn on the atom_pair, angle, and dihedral constraint weights unless otherwise on.
			if(sfxn->get_weight(atom_pair_constraint) < 1.0e-6) sfxn->set_weight(atom_pair_constraint, 1.0);
			if(sfxn->get_weight(angle_constraint) < 1.0e-6) sfxn->set_weight(angle_constraint, 1.0);
			if(sfxn->get_weight(dihedral_constraint) < 1.0e-6) sfxn->set_weight(dihedral_constraint, 1.0);
		}
	}

	//Message about splitting processes:
	if(procnum==0 && option[v_split_replicates_over_processes]()) {
		printf("The -v_split_replicates_over_processes option was used.  Replicates will be spit over processes.\n"); fflush(stdout);
	} else {
		if(procnum==0 && option[v_consider_all_replicates]()) {
			printf("Error!  If the -v_consider_all_replicates option is used, then the -v_split_replicates_over_processes option must also be used."); fflush(stdout); exit(1);
		}
	}

	//Message about saving memory:
	if(procnum==0 && option[v_savememory]()) {
		printf("The -v_savememory option was used.  Only the master proc will store state information.  This will be tranmitted on-the-fly to slave procs.\n");  fflush(stdout);
	} else if (procnum==0) {
		printf("The -v_savememory option was not used.  All procs will store state information.  Warning -- memory could fill up!\n");  fflush(stdout);
	}

	//Message about PCA limit:
	if(procnum==0 && option[v_PCAlimit].user()) {
		if(option[v_PCAlimit]() < 1) {
			printf("Error!  The -v_PCAlimit flag value cannot be less than 1.\n");  fflush(stdout); exit(1);
		} else {
			printf("The -v_PCAlimit flag was used. %lu PCA vector(s) will be stored per state.\n", (core::Size)option[v_PCAlimit]());  fflush(stdout);
		}
	}

	//Check that amino acids have been specified, if the -v_allowed_aminoacids flag is used.
	if(procnum==0 && option[v_allowed_aminoacids].user() && option[v_allowed_aminoacids]().length()==0 ) {
		printf("Error!  If the -v_allowed_aminoacids flag is used, a list of allowed amino acid one-letter codes must be provided (with no spaces -- e.g. \"ADEFHKLNPQRSTVW\").\nCrashing.\n");
		fflush(stdout); exit(1);
	}

	//Check that beta amino acids have been specified, if the -v_allowed_betas flag is used.
	if(procnum==0 && option[v_allowed_betas].user() && option[v_allowed_betas]().length()==0 ) {
		printf("Error!  If the -v_allowed_betas flag is used, a list of allowed beta-3-amino acid one-letter codes must be provided (with no spaces -- e.g. \"ADEFKLNPQRSTVW\").\nCrashing.\n");
		fflush(stdout); exit(1);
	}

	//Message about disulfides:
	if(procnum==0 && option[v_disulfide_positions]().size() > 0) {
		if (option[v_disulfide_positions]().size() % 2 != 0) {
			printf("Error!  PAIRS of disulfide-bonded cysteines must be specified if the -v_disufide_positions flag is used.  Crashing.\n");
			fflush(stdout); exit(1);
		} else {
			printf("Disulfide-bonded residues were specified with the -v_disulfide_positions flag:\n");
		}
	}
	utility::vector1 < std::pair < core::Size, core::Size > > disulf_positions; //List of pairs of disulfide-bonded cysteines.
	if(option[v_disulfide_positions].user()) {
		for(core::Size i=1; i<=option[v_disulfide_positions]().size(); i+=2) {
			disulf_positions.push_back( std::make_pair( (core::Size)option[v_disulfide_positions]()[i], (core::Size)option[v_disulfide_positions]()[i+1] ) ); //Store the list of disulfide-bonded cysteines.
			if(procnum==0) printf("\tCYS%lu--CYS%lu\n", (core::Size)option[v_disulfide_positions]()[i], (core::Size)option[v_disulfide_positions]()[i+1] );
		}
	}
	if(procnum==0) fflush(stdout);

	//Alter the scoring function for v_cyclic:
	if(option[v_cyclic]()) {
		if(procnum==0) printf("Setting constraint weights for a peptide bond between the N- and C-termini (-v_cyclic flag).\n");
		else { //The scorefunction doesn't exist in proc 0; only the other procs can set it.
			sfxn->set_weight(atom_pair_constraint, 1.0);
			sfxn->set_weight(dihedral_constraint, 1.0);
			sfxn->set_weight(angle_constraint, 1.0);
		}
		if(procnum==0 && option[v_use_cyclic_permutations]()) {
			printf("Cyclic permutations of backbone conformations will be tried (-v_use_cyclic_permutations flag).  An offset of %lu residues will be used.\n", (core::Size)option[v_cyclic_permutation_offset]());
		}
	} else {
		if(procnum==0 && option[v_use_cyclic_permutations]()) {
			printf("Error!  The -v_use_cyclic_permutations flag may only be used with the -v_cyclic flag.  Crashing.\n");
			fflush(stdout); exit(1);
		}
	}

	//Check for extra atoms to use in RMSD calculation:
	utility::vector1 <core::id::NamedAtomID> extra_rms_atoms;
	if(option[v_extra_rms_atoms].user()) {
		parse_extra_atom_list(extra_rms_atoms);
		if(procnum==0) {
			printf("The following additional atoms will be used in the RMSD calculation:\nRESIDUE\tATOM\n");
			for(core::Size i=1, imax=extra_rms_atoms.size(); i<=imax; ++i) {
				printf("%lu\t%s\n", extra_rms_atoms[i].rsd(), extra_rms_atoms[i].atom().c_str());
			}
			fflush(stdout);
		}
	}

	//Checks for PDB files:
	//utility::vector1 < std::string> negfiles; //All of the input PDBs are NEGATIVE states.
	if(!option[in::file::s].user() && !option[in::file::l].user() && !option[in::file::silent].user()) {
			printf("ERROR!  The user must specify input PDB files for the negative design states with the -in:file:s or -in:file:l flags, or silent files with the -in:file:silent flag.\nCrashing.\n");
			fflush(stdout);
			exit(1);
	}

	//if(option[in::file::s].user()) { //If PDBs have been specified with -in:file:s
	//	negfiles = option[in::file::s]();
	//} else if(option[in::file::l].user()) { //If PDBs have been specified with -in:file:l
	//	std::string pdbfileslistfile = option[in::file::l]()[1];
	//	if(procnum==0 || !option[v_savememory]()) parsefilelist(pdbfileslistfile, negfiles); //Put the list of PCA files into pcafilenames.
	//} else {}

	//Checks about PCA files:
	utility::vector1 <std::string> pcafilenames; //Vector to store all the PCA file names.
	if(option[v_PCAfiles]().size() > 0 || option[v_PCAfiles_list].user()) { //If we're using PCA files
		if(option[v_PCAfiles]().size() > 0) { //If PCAfiles were specified with the v_PCAfiles flag
			if(procnum==0 || !option[v_savememory]()) pcafilenames=option[v_PCAfiles](); //Copy the file list
		} else { //If PCA files were specified with the v_PCAfiles_list option
			if(procnum==0 || !option[v_savememory]() ) {
				std::string pcafileslistfile = option[v_PCAfiles_list]();
				parsefilelist(pcafileslistfile, pcafilenames); //Put the list of PCA files into pcafilenames.
			}
		}

		if(procnum==0) { //Checks done by proc 0:
			//This check is now done at file load.
			//if (pcafilenames.size() != negfiles.size()) {
			//	printf("ERROR!  If PCA files are specified, one such file must be specified for each negative state (and in the same order).\nCrashing.\n");
			//	fflush(stdout);
			//	exit(1);
			//}
			if (!option[v_PCAfile_native].user() ) {
				printf("ERROR!  If PCA files are specified, one must be specifed for the native state, too.\nCrashing.\n");
				fflush(stdout);
				exit(1);
			}
		}
	} else if(procnum==0 && option[v_PCAfile_native].user() && option[v_PCAfiles]().size()==0 && !option[v_PCAfiles_list].user()) {
		printf("ERROR!  If a PCA file is specified for the native state, PCA files must be specifed for all of the other states, too.\nCrashing.\n");
		fflush(stdout);
		exit(1);
	}

	//Variables for making random mutations:
	core::Size target_pos = 0;
	string all_aa = "ADEFHIKLMNPQRSTVWY"; //Possible mutations -- no cysteines or glycines by default!
	string all_betas = "ADEFIKLMNPQRSTVWY"; //Possible mutations for beta-amino acids -- no cysteines, glycines, or histidines by default!
	if(option[v_allowed_aminoacids].user()) all_aa=option[v_allowed_aminoacids](); //If the user has specified allowed amino acid mutations, use this list.
	else if(option[v_allow_gly_mutation]()) all_aa+="G"; //Add glycine to the possible mutations if the option to do so is set to "true" and the user hasn't specified a list.
	if(option[v_allowed_betas].user()) all_betas=option[v_allowed_betas](); //If the user has specified allowed beta-3-amino acid mutations, use this list.
	else if(option[v_allow_gly_mutation]()) all_betas+="G"; //Add glycine to the possible beta-mutations if the option to do so is set to "true" and the user hasn't specified a list.
	char target_aa = 'X';

	//The FastRelax mover that I'll use:
	protocols::relax::FastRelax frlx(sfxn, option[v_relaxrounds]()); //Create the FastRelax mover and set it to do a number of repeats given by option[v_relaxrounds]() (and give it the score function).

	{	// Add a taskoperation to frlx indicating which residues can be moved if the
		// -v_norepack_positions option has been used.
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		TaskFactoryOP frlx_tasks = new TaskFactory();
		frlx_tasks->push_back(new RestrictToRepacking());
		if(option[v_norepack_positions]().size()>0) {
			//If the user has specified that certain positions should not repack, set up a task operation here.
			PreventRepackingOP norepack = new PreventRepacking();
			for(core::Size i=1, imax=option[v_norepack_positions]().size(); i<=imax; ++i) {
				norepack->include_residue( option[v_norepack_positions]()[i] );
			}
			frlx_tasks->push_back( norepack );
		}
		frlx.set_task_factory( frlx_tasks );
	}

	//Out file prefix:
	string outprefix;
	if(option[out::file::o].user()) outprefix = option[out::file::o]();
	else outprefix = "out";

	std::string posfile = ""; if(option[in::file::native].user()) posfile=option[in::file::native](); //The "native" input PDB is the POSITIVE state.
	std::string posstruct_tag = option[v_native_tag](); //Alternatively, the index of the structure to be used as the positive state can be specified with this flag.
	if(procnum==0) { //Check that one of these has been specified:
		runtime_assert_string_msg( posfile!="" || posstruct_tag!="", "The positive state must be supplied as a PDB file with the -in:file:native flag, or its tag in the input silent file must be specifed with the -v_native_tag flag." );
	}

	utility::vector1 < core::io::silent::SilentStructOP > allstates_master; //Minimial pose data for all conformations (with the positive in position 1), unrelaxed.  This is a list of binary silent structure owning pointers.  Binary silent structures will be converted to strings for transmission to/from slave procs.
	//utility::vector1 < utility::vector1 < core::Real > > allstates_master; //All conformations (with the positive in position 1), unrelaxed.  This is a list of vectors of backbone (phi, psi, omega) values.
	//utility::vector1 < utility::vector1 < utility::vector1 < core::Real > > > allstates_master_sidechains; //The chi values of side chain residues.  This is a list of states, containing lists of amino acid positions, containing lists of chi values.
	utility::vector1 < string > allfiles; //All of the file names.
	utility::vector1 < utility::vector1 < utility::vector1 < core::Real > > > PCAmatrices; //A list of PCA matrices, one for each state.  This list remains empty unless PCA files are specified.
	utility::vector1 < utility::vector1 < core::Real > > PCAvariances; //A list of variance vectors, one for each state.  This list remains empty unless PCA files are specified.
	core::Size PCAvectorcomponentcount=0; //The number of entries in each PCA vector in the PCA matrices.

	core::pose::Pose positive_master; //A pose to store the positive backbone configuration.
	core::io::silent::SilentStructOP positive_master_silentstruct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");

	char curstructtag_positive[256];
	if(procnum == 0 ) {
		if(posfile!="") {
			printf("Importing POSITIVE backbone configuration from %s.\n", posfile.c_str()); fflush(stdout);
			sprintf(curstructtag_positive, "POSITIVE_%s", posfile.c_str());
			//Only the master proc imports the positive state, then transmits it to all others.
			core::import_pose::pose_from_pdb(positive_master, posfile);
			//Store the file name:
			allfiles.push_back(posfile);
		} else if (posstruct_tag!="") { //If the index of the positive state in the silent files has been specified.
			//core::Size curstruct_index = 1;
			bool foundit = false;
			sprintf(curstructtag_positive, "POSITIVE_%s", posstruct_tag.c_str());
			core::import_pose::pose_stream::MetaPoseInputStream input_stream = core::import_pose::pose_stream::streams_from_cmd_line( false ); //Get all inputs, NOT renumbering decoys.
			utility::vector1 < core::import_pose::pose_stream::PoseInputStreamOP > posestreams = input_stream.get_input_streams(); //Split the list of input streams.
			for(core::Size i=1, imax=posestreams.size(); i<=imax; ++i) {
				core::import_pose::pose_stream::SilentFilePoseInputStreamOP silentstream ( dynamic_cast<core::import_pose::pose_stream::SilentFilePoseInputStream*>(posestreams[i].get()) );
				core::import_pose::pose_stream::PDBPoseInputStreamOP pdbstream ( dynamic_cast<core::import_pose::pose_stream::PDBPoseInputStream*>(posestreams[i].get()) );
				if(silentstream) {
					for(core::Size ii=1, iimax=silentstream->silent_file_data()->structure_list().size(); ii<=iimax; ++ii) {
						core::io::silent::SilentStructOP curstruct = silentstream->silent_file_data()->structure_list()[ii];
						//printf("curstruct->decoy_tag()=%s\n", curstruct->decoy_tag().c_str()); fflush(stdout); //DELETE ME
						if(curstruct->decoy_tag() == posstruct_tag) {
							printf("Importing POSITIVE backbone configuration with tag \"%s\".\n", curstruct->decoy_tag().c_str()); fflush(stdout);
							curstruct->fill_pose( positive_master );
							foundit=true;
							break;
						}
					}
					silentstream->reset();
				} else if (pdbstream) {
					utility::vector1<utility::file::FileName> filenames = pdbstream->get_filenames();
					for(core::Size i=1, imax=filenames.size(); i<=imax; ++i) {
						if(filenames[i].base() == posstruct_tag) {
							printf("Importing POSITIVE backbone configuration from %s.", filenames[i].name().c_str());
							core::import_pose::pose_from_pdb(positive_master, filenames[i]);
							foundit=true;
							break;
						}
					}
				}
			}
			runtime_assert_string_msg(foundit, "Positive pose tag was not found!");
			//Store the tag name:
			allfiles.push_back(posstruct_tag);
		}
		positive_master_silentstruct->fill_struct(positive_master, curstructtag_positive);
		//printf("Master proc transmitting positive state.\n"); fflush(stdout); //DELETE ME
		transmit_silent_struct (0, 0, 0, positive_master_silentstruct, true); //BROADCAST the silent struct
	} else { //All other procs receive the positive state from the master proc
		//printf("Proc %i waiting to receive...\n", procnum); fflush(stdout); //DELETE ME
		receive_silent_struct (positive_master, 0, procnum, procnum, true); //Receive the silent struct via MPI_Bcast
		//printf("Proc %i has received master state.\n", procnum); fflush(stdout); //DELETE ME
	}

	//Temporary debugging output: dump pdbs for all procs
	//char outfilename_temp[256];                                   //DELETE ME
	//sprintf(outfilename_temp, "POSITIVE_%04i.pdb", procnum);      //DELETE ME
	//positive_master.dump_pdb(outfilename_temp);                   //DELETE ME

	//Barrier for positive state load, just to be on the safe side:
	MPI_Barrier(MPI_COMM_WORLD);

	if(option[v_strip_terminal_gly]()) {//If we're stripping off terminal glycines:
		bool glystripped = strip_terminal_glycines(positive_master);
		if(glystripped && procnum==0) { printf("Removed terminal glycine residues.\n"); fflush(stdout); }
	}

	//Set the termini depending on whether this is a cyclic peptide:
	if(option[v_cyclic]()) {
		core::pose::remove_lower_terminus_type_from_pose_residue(positive_master, 1);
		core::pose::remove_upper_terminus_type_from_pose_residue(positive_master, positive_master.n_residue());
	} else {
		set_up_termini(positive_master);
	}

	//Add appropriate covalent connections if a covalent connections file is specified:
	utility::vector1 < std::pair < core::id::NamedAtomID, core::id::NamedAtomID > > covalentconnections;
	if(option[v_covalent_connections].user()) {
		read_covalent_connections_file(covalentconnections, static_cast<core::Size>(procnum)); //Function reads covalent connections file and populates covalent connections list (done ONCE).
		set_up_covalent_connections(positive_master, covalentconnections); //Function reads covalent connections list and sets up connections.
	}

	//If the user has specified a constraints file, here's the place to add it:
	//add_user_constraints(positive_master); //Checks internally for user options.

	if (procnum==0 || !option[v_savememory]()) {
		//Store backbone conformation:
		//allstates_master.push_back(positive_master_silentstruct);

		//Read PCA matrix for native state, if specified:
		if(option[v_PCAfile_native].user()) {
			std::string PCAfile_native = option[v_PCAfile_native]();
			utility::vector1 <core::Real> variance_array;
			utility::vector1 < utility::vector1 < core::Real > > PCA_matrix;
			core::Size vectcounttemp = read_PCAfile ( PCAfile_native.c_str(), variance_array, PCA_matrix, true);
			if(vectcounttemp>0) PCAvectorcomponentcount=vectcounttemp; //Store the number of entries in each PCA vector
			PCAvariances.push_back(variance_array);
			PCAmatrices.push_back(PCA_matrix);
		}
	}

	if (procnum == 0) { printf("Import complete.  Starting sequence is %s.\n", positive_master.sequence().c_str()); fflush(stdout); }
	strcpy (currentseq,positive_master.sequence().c_str()); //Copy the sequence into currentseq

	//Check whether the user has provided a list of equivalent positions that must be mutated together:
	utility::vector1 < utility::vector1 < core::Size > > equivpositions; //List of lists of equivalent positions
	if(procnum==0 && option[v_equivalent_positions].user()) {
		std::string equivfile = option[v_equivalent_positions]();
		read_equivalent_position_file( equivfile.c_str(), equivpositions);
		printf("The following positions are equivalent, and will always be mutated in tandem to the same amino acid residue:\n");
		if(equivpositions.size()>0) {
			for(core::Size i=1, imax=equivpositions.size(); i<=imax; i++) {
				if(equivpositions[i].size()>0) {
					for(core::Size j=1, jmax=equivpositions[i].size(); j<=jmax; j++) {
						if(j>1 && jmax > 2) printf(", ");
						else if (j>1) printf (" ");
						if(j>1 && j==jmax) printf("and ");
						printf("%lu", equivpositions[i][j]);
					}
					printf(" are equivalent\n");
				}
			}
		} else {
			printf("No equivalent positions found.\n");
		}
		fflush(stdout);
	}

	//Check whether the -v_preserve_chirality flag is used.  If it is, store a list of all of the D-amino acid positions.
	if (procnum == 0 && option[v_preserve_chirality]()) printf("-v_preserve_chirality flag used.  The following positions are D-amino acid positions: ");
	utility::vector1 <bool> dpositions;
	for(core::Size ir=1, nres=positive_master.n_residue(); ir<=nres; ir++) dpositions.push_back( core::chemical::is_canonical_D_aa( positive_master.residue(ir).aa() ) );
	if(procnum==0) {
		core::Size dpositioncount = 0;
		for(core::Size i=1; i<=dpositions.size(); i++) {
			if(dpositions[i]) {
				printf("%lu ", i);
				dpositioncount++;
			}
		}
		if(dpositioncount>0) {
			printf("\n");
		} else {
			printf("No D-amino acid positions found.\n");
		}
		fflush(stdout);
	}

	//Store a list of beta-amino acid positions:
	utility::vector1 <bool> betapositions; //The list of beta-amino acid positions (one entry per residue, true for beta, false otherwise).
	bool hasbetaresidues=false;
	for(core::Size ir=1, nres=positive_master.n_residue(); ir<=nres; ir++) {
		betapositions.push_back( positive_master.residue(ir).type().is_beta_aa() );
		if( positive_master.residue(ir).type().is_beta_aa() ) hasbetaresidues=true;
	}
	if(procnum==0 && hasbetaresidues) {
		printf("The following positions are beta-3-amino acids: ");
		for(core::Size ir=1, imax=betapositions.size(); ir<=imax; ir++) {
			if(betapositions[ir]) printf("%lu ", ir);
		}
		printf("\n"); fflush(stdout);
	}

	if(option[v_randseq]() && procnum == 0) //If the user has turned on the random sequence flag, randomize the sequence.
	{
		printf("-vrandseq flag used.  Randomizing sequence.  "); fflush(stdout);
		for(core::Size ir = 1; ir<=positive_master.n_residue(); ir++)
		{
			//if(option[v_preserve_cys]() && positive_master.residue(ir).name1()=='C') continue; //Don't mutate cysteines if we're supposed to preserve them.
			if(option[v_ignoreresidue]().size()>0) { //Don't change residues in the ignore list.
				bool continue_on = false;
				for(core::Size i=1; i<=option[v_ignoreresidue]().size(); i++) {
					if((core::Size)option[v_ignoreresidue]()[i]==ir) continue_on=true;
				}
				if(continue_on) continue;
			}
			if(option[v_norepack_positions]().size()>0) { //Don't change residues in the list of residues that can't repack.
				bool continue_on = false;
				for(core::Size i=1, imax=option[v_norepack_positions]().size(); i<=imax; i++) {
					if((core::Size)option[v_norepack_positions]()[i]==ir) continue_on=true;
				}
				if(continue_on) continue;
			}
			if(betapositions[ir]) { //beta-amino acids:
				target_aa=all_betas[RG.random_range(0,all_betas.length()-1)]; //Choosing from the allowed mutations.
			} else { //default -- alpha-amino acids (D and L)
				target_aa=all_aa[RG.random_range(0,all_aa.length()-1)]; //Choosing from the allowed mutations.
			}
			currentseq[ir-1] = target_aa;
		}
	}
	else //If the user hasn't specified that a random sequence be used, check for a FASTA file.
	{
		if (option[in::file::fasta].user() && procnum==0)
		{
			string fastaseq = core::sequence::read_fasta_file (option[in::file::fasta]()[1])[1]->sequence();
			printf("Fasta file specified.  Read sequence %s.\n", fastaseq.c_str()); fflush(stdout);
			for(core::Size ir = 1; ir<=positive_master.n_residue(); ir++){
				currentseq[ir-1]=fastaseq[ir-1];
			}
			mutate_to_sequence(currentseq, positive_master, dpositions, betapositions, true); //Mutate the positive master pose in proc 0 to the starting sequence
		}		
	}
	if (procnum==0) {printf("Starting sequence is %s (%lu residues).\n", currentseq, (core::Size)strlen(currentseq)); fflush(stdout);}

	if (procnum==0 || !option[v_savememory]()) {
		//Store backbone conformation:
		//We need to re-generate the silent structure because the sequence might have changed due to a FASTA file.
		core::io::silent::SilentStructOP temp_master_silentstruct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
		temp_master_silentstruct->fill_struct(positive_master, curstructtag_positive);
		allstates_master.push_back(temp_master_silentstruct);
	}


	MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side

	//Broadcast the starting sequence to all other group members:
	MPI_Bcast( &currentseq[0], strlen(currentseq)+1, MPI_CHAR, 0, MPI_COMM_WORLD); //strlen(currenseq) is the same on all procs at this point
	const std::string startingsequence = currentseq;
	std::string currentsequence = currentseq;
	std::string lowestEsequence = currentseq;
	std::string lastacceptsequence=currentseq;
	
	//if (procnum != 0) printf ("Proc %i received %s.\n", procnum, startingsequence.c_str()); //Uncomment to check that each proc is receiving the proper starting sequence.

	MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side

	//More variables for random mutations:
	utility::vector1 < core::Size > all_positions; //A vector of positions that can be mutated.
	for (core::Size i = 1; i<=startingsequence.length(); i++) {
		bool addme=true;
		if(option[v_ignoreresidue]().size() > 0) {
			for(core::Size j=1; j<=option[v_ignoreresidue].size(); j++) {
				if((core::Size)option[v_ignoreresidue][j] == i) {addme = false; break;}
			}
		}
		if(addme && option[v_norepack_positions]().size() > 0) {
			for(core::Size j=1, jmax=option[v_norepack_positions]().size(); j<=jmax; ++j) {
				if((core::Size)option[v_norepack_positions]()[j]==i) { addme=false; break; }
			}
		}
		if(addme) all_positions.push_back(i);
	}
	if (procnum==0) {
		printf("Mutations are possible at these positions:");
		for (core::Size i = 1; i <= all_positions.size(); i++) printf(" %lu", all_positions[i]);
		printf("\n");
		printf("Mutations to these alpha-amino acids are possible: %s\n", all_aa.c_str());
		if(hasbetaresidues) printf("Mutations to these beta-3-amino acids are possible: %s\n", all_betas.c_str());
	}

	MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side

	if (procnum==0) printf("\nImporting states for NEGATIVE design:\n");
	if (procnum==0 || !option[v_savememory]()) {
		core::import_pose::pose_stream::MetaPoseInputStream input_stream = core::import_pose::pose_stream::streams_from_cmd_line( false ); //Get all inputs, NOT resorting the decoys
		utility::vector1 < core::import_pose::pose_stream::PoseInputStreamOP > posestreams = input_stream.get_input_streams(); //Split the list of input streams.

		core::Size ifile=0;

		for(core::Size istream=1, istreammax=posestreams.size(); istream<=istreammax; ++istream) { //Loop through all input streams

			//For each input stream, store a list of tags (from silent files) or a list of filenames (pdbs)
			utility::vector1 <std::string> negstate_names;
			core::import_pose::pose_stream::SilentFilePoseInputStreamOP silentstream ( dynamic_cast<core::import_pose::pose_stream::SilentFilePoseInputStream*>(posestreams[istream].get()) );
			core::import_pose::pose_stream::PDBPoseInputStreamOP pdbstream ( dynamic_cast<core::import_pose::pose_stream::PDBPoseInputStream*>(posestreams[istream].get()) );
			if(silentstream) {
				//Do nothing for a silentstream.  The tags() function doesn't seem to work properly, so we'll get tags further down.
				//negstate_names=silentstream->tags();
			} else if (pdbstream) {
				utility::vector1<utility::file::FileName> filenames = pdbstream->get_filenames();
				for(core::Size i=1, imax=filenames.size(); i<=imax; ++i) {
					negstate_names.push_back(filenames[i].base());
				}
			}


			core::Size posecount = 0;
			while( posestreams[istream]->has_another_pose() ) { //Loop through all input structures in this input stream
				++ifile; //Count states
				++posecount; //Count poses within this stream

				//The name of this negative state.
				std::string negstate_name = "";
				if(negstate_names.size() >= posecount) {
					negstate_name=negstate_names[posecount];
				} else {
					char negstatechar [256];
					sprintf(negstatechar, "NEG_%lu", ifile);
					negstate_name = negstatechar;
				}

				if (procnum==0) printf("\n--BEGIN negative state %lu (%s) --\n", ifile, negstate_name.c_str());			

				core::pose::Pose importpose; //A pose to store the pose as it is imported.
				if (procnum==0) printf("\tImporting %s.\n", negstate_name.c_str());
				//core::import_pose::pose_from_pdb(importpose, negfiles[ifile]);
				if(silentstream) {
					silentstream->silent_file_data()->structure_list()[posecount]->fill_pose( importpose );
					silentstream->next_struct();
				} else posestreams[istream]->fill_pose(importpose);
				if(silentstream) {
					std::string tagstring = silentstream->silent_file_data()->structure_list()[posecount]->decoy_tag();
					char negstatechar[256];
					sprintf(negstatechar, "%s_%s", negstate_name.c_str(), tagstring.c_str());
					negstate_name=negstatechar;
					printf( "\tRead tag %s.  State name is now %s.\n", tagstring.c_str(), negstate_name.c_str() ); fflush(stdout);
				}

				if(option[v_strip_terminal_gly]()) {//If we're stripping off terminal glycines:
					bool glystripped = strip_terminal_glycines(importpose);
					if(glystripped && procnum==0) printf("\tRemoved terminal glycine residues.\n");
				}

				//Set the termini depending on whether this is a cyclic peptide:
				if(option[v_cyclic]()) {
					core::pose::remove_lower_terminus_type_from_pose_residue(importpose, 1);
					core::pose::remove_upper_terminus_type_from_pose_residue(importpose, importpose.n_residue());
				} else {
					set_up_termini(importpose);
				}

				//Add appropriate covalent connections if a covalent connections file is specified:
				if(option[v_covalent_connections].user()) {
					set_up_covalent_connections(importpose, covalentconnections); //Function reads covalent connections list and sets up connections.
				}

				//If the user has specified a constraints file, here's the place to add it:
				//add_user_constraints(importpose); //Checks internally for user options.

				core::Real rmsval=get_distance_measure(importpose, positive_master, extra_rms_atoms);
				//printf("Proc %i: rms=%.2f\n", procnum, rmsval); fflush(stdout); //DELETE ME -- for debugging
				if(option[v_removelowrms]() && (rmsval < option[v_removelowrms_threshhold]()) ) {
					if (procnum==0) {
						printf("\tNegative state %lu has a backbone RMSD of %.2f from the positive design state, which is less than the removal threshhold of %.2f.  Ignoring and moving on.\n", ifile, rmsval, option[v_removelowrms_threshhold]());
					} else {
						printf("\tProc %i also discarded negative state %lu.\n", procnum, ifile); fflush(stdout);
					}
					if(!option[v_savememory]()) MPI_Barrier(MPI_COMM_WORLD); //For testing only -- probably not needed.
					if (procnum==0) {printf("--END negative state %lu--\n", ifile); fflush(stdout);}
					continue; //Don't add this state to the set of negative states.
				}

				if(procnum==0 && option[in::file::fasta].user()) {
					mutate_to_sequence(currentseq, importpose, dpositions, betapositions, true); //Mutate the negative state pose in proc 0 to the starting sequence
				}

				//Store the state descriptor:
				//allfiles.push_back(negfiles[ifile]);
				allfiles.push_back(negstate_name);

				//Store backbone conformation:
				//utility::vector1 < core::Real > bb_vect;
				//store_backbone(importpose, bb_vect);
				core::io::silent::SilentStructOP negstate_silentstruct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
				negstate_silentstruct->fill_struct(importpose, negstate_name);
				allstates_master.push_back(negstate_silentstruct);

				//Read PCA matrix for native state, if specified:
				if(pcafilenames.size()>0) {
					//printf("ifile=%lu pcafilenames.size()=%lu\n", ifile, pcafilenames.size()); fflush(stdout); //DELETE ME
					runtime_assert_string_msg(pcafilenames.size() >= ifile, "If PCA files are specified, one must be specified for each negative state!  Crashing!");
					std::string PCAfile = pcafilenames[ifile];
					if (procnum==0) printf("\tImporting PCA vectors from %s.\n", PCAfile.c_str());
					utility::vector1 <core::Real> variance_array;
					utility::vector1 < utility::vector1 < core::Real > > PCA_matrix;
					core::Size vectcounttemp = read_PCAfile ( PCAfile.c_str(), variance_array, PCA_matrix, true);
					if(vectcounttemp>0) PCAvectorcomponentcount=vectcounttemp; //Store the number of entries in each PCA vector
					PCAvariances.push_back(variance_array);
					PCAmatrices.push_back(PCA_matrix);
				}

				if (procnum==0) printf("--END negative state %lu (%s) --\n\n", ifile, negstate_name.c_str());

			} //End looping through states in this stream
		} //End looping through all streams
	} //End if (procnum==0 || !option[v_savememory]())

	//If v_use_cyclic_permutations is used, define additional states:
	if(option[v_use_cyclic_permutations]() && (procnum==0 || !option[v_savememory]())) {
		printf("The cyclic permutations option has been temporarily disabled, since we are now storing and transmitting silent structures rather than backbone dihedral vectors.\n");
		printf("This will be enabled again at some point in the future, but written more efficiently.\n");  fflush(stdout);
		exit(1);
		/*if(procnum==0) printf("Setting up additional states for cyclic permutations (-v_use_cyclic_permutations flag).\n"); fflush(stdout);
		core::Size basestatecount = allstates_master.size();
		char filename[1024];
		std::string fname;
		for(core::Size state=1; state<=basestatecount; state++) {
			for(core::Size offset=1; offset<startingsequence.length(); offset+=option[v_cyclic_permutation_offset]()) {
				//Store additional file names:
				sprintf(filename, "%s_%lu", allfiles[state].c_str(), offset);
				fname = filename;
				allfiles.push_back(fname);

				//Store additional backbone vectors -- MEMORY-INEFFICIENT!:
				utility::vector1 < core::Real > shifted_bb;
				shifted_bb.resize(allstates_master[state].size());
				for(core::Size i=1; i<=shifted_bb.size(); i++) { //Loop through all elements in i
					signed int i2 = (signed int)i-3*(signed int)offset;
					if(i2<=0) i2+=shifted_bb.size();
					shifted_bb[i]=allstates_master[state][i2];
				}
				allstates_master.push_back(shifted_bb);

				//DO NOT store additional side-chain vectors.  For now, v_preserve_cys and v_use_cyclic_permutations are mutually incompatible.

				//Store additional PCA vectors and variances vectors:
				//MEMORY-INEFFICIENT!  Revisit later!
				utility::vector1 < utility::vector1 < core::Real > > shifted_PCAmatrix;
				PCAvariances.push_back(PCAvariances[state]); //Duplicate the PCA variances.  They don't need to be offset.
				for(core::Size i=1; i<=PCAmatrices[state].size(); i++) { //Loop through all PCA vectors for the state
					utility::vector1 <core::Real> shifted_PCAvector;
					shifted_PCAvector.resize(PCAmatrices[state][i].size());
					for(core::Size j=1; j<=PCAmatrices[state][i].size(); j++) { //Loop through all the elements of the ith PCA vector.
						signed int j2=(signed int)j-3*(signed int)offset;
						if(j2<=0) j2+=PCAmatrices[state][i].size();
						shifted_PCAvector[j]=PCAmatrices[state][i][j];
					}
					shifted_PCAmatrix.push_back(shifted_PCAvector);
				}
				PCAmatrices.push_back(shifted_PCAmatrix);
			}
		}*/
	}

	//Count the states and send the count to the other procs:
	core::Size totalstatecount;
	if(procnum==0) {
		totalstatecount = allstates_master.size();
		int tempstatecount=(int)totalstatecount;
		MPI_Bcast(&tempstatecount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		int tempstatecount;
		MPI_Bcast(&tempstatecount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		totalstatecount = (core::Size)tempstatecount;		
	}
	//printf("(Proc %i) Total state count = %i\n", (int)procnum, (int)totalstatecount); fflush(stdout); //DELETE ME

	//Send the number of entries in each PCA matrix to all the other procs
	if(procnum==0) {
		int compcount=(int)PCAvectorcomponentcount;
		MPI_Bcast(&compcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		int compcount;
		MPI_Bcast(&compcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		PCAvectorcomponentcount = (core::Size)compcount;		
	}

	//Additional variables:
	utility::vector1 < core::Real > energyvals_initial; //Energies for each state (initial sequence).  Only the master proc will hold all energies; all slaves will dynamically resize this array for the subset of states that it calculates.
	utility::vector1 < core::Real > energyvals_current; //Energies for each state (current sequence)
	utility::vector1 < core::Real > energyvals_lowestE; //Energies for each state (lowest energy sequence)
	utility::vector1 < core::Real > energyvals_lastaccept; //Energies for each state (last accepted sequence)
	utility::vector1 < core::Real > rmsvals_initial; //RMS values for each state (initial sequence), measured from target.  Only the master proc will hold this list.
	utility::vector1 < core::Real > rmsvals_current; //RMS values for each state (current sequence), measured from target.  Only the master proc will hold this list.
	utility::vector1 < core::Real > rmsvals_lowestE; //RMS values for each state (lowest energy sequence), measured from target.  Only the master proc will hold this list.
	utility::vector1 < core::Real > rmsvals_lastaccept; //RMS values for each state (last accepted sequence), measured from target.  Only the master proc will hold this list.
	core::pose::Pose positive_current; //The current energy-minimized positive state (with the current sequence)
	core::pose::Pose positive_lowestE; //The lowest-energy energy-minimized positive state (with the lowest-energy sequence)
	core::pose::Pose positive_lastaccept; //The last accepted energy-minimized positive state (with the last accepted sequence)
	if (procnum==0 || !option[v_savememory]()) {
		positive_current = positive_master;
		positive_lowestE = positive_master;
		positive_lastaccept = positive_master;
	}

	MPI_Barrier(MPI_COMM_WORLD); //To be on the safe side, wait until all procs get to this point.

	//At this point, we're set up and ready to go.  First, we'll score the native state, and all of the negative states.
	if (procnum==0) printf("Scoring starting sequence.\n");
	scoreall(procnum, totalprocs, totalstatecount, PCAvectorcomponentcount,
		currentsequence, dpositions, betapositions, allstates_master, disulf_positions,
		energyvals_initial, rmsvals_initial, positive_master, positive_current, RG,
		frlx, sfxn, PCAvariances, PCAmatrices, extra_rms_atoms);

	//Vars for sequence scores:
	core::Real sequence_score_initial=0.0;
	core::Real prob_initial=0.0;
	core::Real sequence_score_current=0.0;
	core::Real sequence_score_lowestE=0.0;
	core::Real sequence_score_lastaccept=0.0;

	//Calculate the score for this sequence, based on the calculated energies:
	if (procnum==0) {
		sequence_score_initial = calc_sequence_score(energyvals_initial, rmsvals_initial);
		prob_initial = exp(-1.0*sequence_score_initial/option[v_kbt]());
		prob_initial = prob_initial / (1.0+prob_initial);
		printf("Initial score = %.6f\n", sequence_score_initial);
		printf("P(init) = %.6f\n", prob_initial);
		printf("STATE\tRMSD\tENERGY\n");
		if(option[v_consider_all_replicates]()) {
			core::Size j=1, k=1;
			for(core::Size i=1, imax=energyvals_initial.size(); i<=imax; ++i) {
				printf("%s_rep%lu\t%.6f\t%.6f\n", utility::file_basename(allfiles[j]).c_str(), k, rmsvals_initial[i], energyvals_initial[i]);
				++k;
				if((i-1) % (core::Size)option[v_MCminimize_replicates]() == 0) { ++j; k=1; }
			}
		} else {
			for(core::Size i=1; i<=allfiles.size(); i++) printf("%s\t%.6f\t%.6f\n", utility::file_basename(allfiles[i]).c_str(), rmsvals_initial[i], energyvals_initial[i]);
		}
		printf("\n\n"); fflush(stdout);

		char outfile[1024];
		sprintf(outfile, "%s_0000.pdb", outprefix.c_str());
		positive_current.dump_pdb(outfile);
		//(*sfxn)(positive_current); //DELETE ME
		//positive_current.dump_scored_pdb(outfile, (*sfxn)); //Ordinarily, don't want to do this, since it necessitates setting up the scorefunction for the master proc, which is a waste of memory.
		printf("Writing out initial positive structure as %s\n", outfile);
	}

	//Store lowest-energy values found so far:
	energyvals_lowestE=energyvals_initial;
	energyvals_lastaccept=energyvals_initial;
	rmsvals_lowestE=rmsvals_initial;
	rmsvals_lastaccept=rmsvals_initial;
	sequence_score_lowestE=sequence_score_initial;
	sequence_score_lastaccept=sequence_score_initial;

	//Monte Carlo search for the best sequence.
	//I will:
	//	--Make a random mutation to positive_trial, and the same mutation to all of the negative_trial states.
	//	--Relax all of these states.
	//	--Score all of these states.
	//	--Calculate my scoring function.
	//	--Accept or reject the mutation based on the Metropolis criterion.	

	bool acceptmutation;
	core::Size acceptcount=0;
	if(procnum==0) {
		if(!option[v_scoreonly]()) printf("\nBeginning Monte Carlo search for the best sequence.\n");
		else printf("\nThe -v_scoreonly option was used.  Initial sequence has been scored.  Exiting.\n");
		fflush(stdout);
	}
	if(!option[v_scoreonly]()) {
		for(core::Size itrial = 1; itrial <= (core::Size)option[v_trialcount](); itrial++) {
			//printf("Proc %i reached barrier.\n", procnum); fflush(stdout);  //DELETE ME
			MPI_Barrier(MPI_COMM_WORLD); //Better safe than sorry.

			if(procnum==0) fflush(stdout); //Flush the printf buffer.
	
			positive_current=positive_master; //Make a copy of the master positive state:

			if(procnum==0) {
				//Picking a random position and mutation:
				//Pick a position of the allowed positions:
				target_pos = all_positions[RG.random_range(1, all_positions.size())];
				while(true) {
					//Pick an amino acid randomly:
					if(betapositions[target_pos]) { //beta-amino acids
						target_aa=all_betas[RG.random_range(0,all_betas.length()-1)]; //Choosing from allowed beta-mutations.
					} else { //default -- alpha-amino acids (D or L)
						target_aa=all_aa[RG.random_range(0,all_aa.length()-1)]; //Choosing from allowed mutations.
					}
					if(currentsequence[target_pos-1]!=target_aa) break; //Continue only if we're not mutating to the same residue.
				}

				currentsequence[target_pos-1]=target_aa; //Make the mutation in the sequence string.

				//If there are equivalent positions, mutate them too:
				core::Size listpos=0;
				if(option[v_equivalent_positions].user()) {
					listpos = is_in_list_of_lists(target_pos, equivpositions);
					if(listpos>0) {
						for(core::Size i=1,imax=equivpositions[listpos].size(); i<=imax; i++) currentsequence[ (equivpositions[listpos][i])-1 ]=target_aa;
					}
				}

				printf ("\nMove %lu: trying effect of %c%lu%c mutation (%s->%s).\n", itrial, lastacceptsequence[target_pos-1], target_pos, target_aa, lastacceptsequence.c_str(), currentsequence.c_str());
				if(listpos>0 && equivpositions[listpos].size()>1) {
					printf("(Also mutating position(s): ");
					for(core::Size i=1,imax=equivpositions[listpos].size(); i<=imax; i++) {
						if(equivpositions[listpos][i]!=target_pos) printf("%lu ", equivpositions[listpos][i]);
					}
					printf(")\n");
				}
				fflush(stdout);
			} //The new sequence is ONLY in the master proc, now.  I need to send it out to the slaves.

			//Send the sequence out to all procs and do the energy-minimization for all states
			scoreall (procnum, totalprocs, totalstatecount, PCAvectorcomponentcount,
				currentsequence, dpositions, betapositions, allstates_master, disulf_positions,
				energyvals_current, rmsvals_current, positive_master, positive_current, RG,
				frlx, sfxn, PCAvariances, PCAmatrices, extra_rms_atoms);

			if(procnum==0) {
				//Calculate the current score for the current sequence:
				sequence_score_current = calc_sequence_score (energyvals_current, rmsvals_current);
				core::Real prob = exp(-1.0*sequence_score_current/option[v_kbt]());
				prob = prob / (1.0+prob);
				printf("Move %lu: Score function = %.6f\nP(current) = %.6f\n", itrial, sequence_score_current, prob);
		
				//Apply the Metropolis criterion:
				acceptmutation = false;
				if(sequence_score_current < sequence_score_lastaccept) {
					acceptmutation = true;
				} else { //if the score is greater than the last one accepted
					//Decide whether or not to accept here based on the Metropolis criterion.
					if( exp(-(sequence_score_current-sequence_score_lastaccept)/option[v_MCtemperature]()) > RG.uniform() ) acceptmutation = true;
				}

				//Accept the mutation.
				if(acceptmutation) {
					printf("Accepting move %lu (%c%lu%c).\n", itrial, lastacceptsequence[target_pos-1], target_pos, target_aa);
					sequence_score_lastaccept=sequence_score_current; //Store the score as the last accepted.
					positive_lastaccept=positive_current; //Store the current pose of the positive state as the last accepted.
					energyvals_lastaccept=energyvals_current; //Store the current list of energy values as the last accepted.
					rmsvals_lastaccept=rmsvals_current; //Store the current list of rms values as the last accepted.
					lastacceptsequence=currentsequence;

					if(sequence_score_lastaccept<sequence_score_lowestE) { //If this is the best-scoring sequence found so far...
						printf("New lowest scorefunction value found with sequence %s\n", currentsequence.c_str());
						sequence_score_lowestE=sequence_score_current; //Store the score as the lowest energy.
						positive_lowestE=positive_current; //Store the current pose of the positive state as the lowest energy.
						energyvals_lowestE=energyvals_current; //Store the current list of energy values as the lowest energy.
						rmsvals_lowestE=rmsvals_current; //Store the current list of rms values as the lowest energy.
						lowestEsequence=currentsequence;
					
						//Print a summary:
						printf("\tSTART\tCURRENT\n");
						printf("Sequence:\t%s\t%s\n", startingsequence.c_str(), currentsequence.c_str());
						printf("FITNESS:\t%f\t%f\n", sequence_score_initial, sequence_score_current);
						printf("P(DesignState):\t%f\t%f\n", prob_initial, prob);
						printf("STATE\tRMSD\tENERGY\n");
						if(option[v_consider_all_replicates]()) {
							core::Size j=1, k=1;
							for(core::Size i=1, imax=energyvals_current.size(); i<=imax; ++i) {
								printf("%s_rep%lu\t%.6f\t%.6f\n", utility::file_basename(allfiles[j]).c_str(), k, rmsvals_current[i], energyvals_current[i]);
								++k;
								if((i-1) % (core::Size)option[v_MCminimize_replicates]() == 0) { ++j; k=1; }
							}
							for(core::Size i=1; i<=allfiles.size(); i++) printf("%s\t%.6f\t%.6f\n", utility::file_basename(allfiles[i]).c_str(), rmsvals_current[i], energyvals_current[i]);
						}
						printf("\n"); fflush(stdout);
		
						//Dump out a PDF:
						acceptcount++;
						char outfile [1024];
						sprintf(outfile, "%s_%04lu.pdb", outprefix.c_str(), acceptcount);
						positive_lowestE.dump_pdb(outfile);
					} //if score_lastaccept<score_lowest
				printf("Sequence is now %s.\n", currentsequence.c_str());
				} //if acceptmutation
				else {
					printf ("Rejecting move %lu (%c%lu%c).\n", itrial, lastacceptsequence[target_pos-1], target_pos, target_aa);
					currentsequence=lastacceptsequence;
				}
			} //All of the acceptance or rejection of the mutation is done in the master proc, and is not transmitted to the slaves.
		} //Looping through trials
	}

	MPI_Barrier(MPI_COMM_WORLD); //Better safe than sorry.
	//TODO: add some output info here, perhaps.
	if(procnum==0) {
		printf("\n*****\nJOB COMPLETED\n*****\n");
		fflush(stdout);
	}

	//Terminate all procs:
	MPI_Finalize();

	return 0;
} //End main()

