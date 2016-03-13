/****************************************************************************************************
	design_truecycpeptide_fromstub.cc

	A pilot app to let me design a cyclic peptide inhibitor in the active site of an enzyme, based
	on starting residues already present in an existing structure.  This is intended for the design
	of a cyclic peptide inhibitor of the antibiotic resistance enzyme New Delhi metallo-beta-lactamase
	1, based on the structure of this enzyme bound to L-captopril (a structural analogue of a
	DCys-LPro dipeptide).
	This program was written by Vikram K. Mulligan, Baker Laboratory.

	USER INPUTS:
	-in::file::s (the input PDB, formatted so that any bound metals are separate chains, as is the
		stub that will be used to grow the cyclic peptide).
	-relaxfirst (Boolean -- do I relax the input PDB first?  Default true.)
	-relaxrounds (Integer -- number of fastrelax rounds.  Default 3.)
	-stubchain (Integer -- the chain that represents the peptide stub that will be extended to
		generate a cyclic peptide.  Default chain 2.)

	HISTORY:
	--File created 21 October 2013 by Vikram K. Mulligan.
	--Tweaked 9 January 2014 to allow enumeration of backbone conformations for a chain of
	L-amino acids only.  TODO: add support for D-amino acids at specific positions, or prolines
	at specific positions.
****************************************************************************************************/

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
//#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.hh>
#include <numeric/model_quality/rms.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/util.hh>
#include <numeric/random/random.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/chemical/AA.hh>

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

using ObjexxFCL::format::F;
using namespace ObjexxFCL;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, acceptfirstmove )
OPT_KEY( StringVector, allowedLaas )
OPT_KEY( StringVector, allowedDaas )
OPT_KEY( Boolean, cyclic )
OPT_KEY( Real, dab_sev_strength )
OPT_KEY( Real, design_farep_mult )
OPT_KEY( Real, design_faatr_mult )
OPT_KEY( Boolean, enumerate_only_mode )
OPT_KEY( Real, jittersize )
OPT_KEY( Integer, looplength )
OPT_KEY( Real, MCtemp )
OPT_KEY( Integer, maxattempts )
OPT_KEY( Integer, maxdesignrounds )
OPT_KEY( Integer, maxiterations )
OPT_KEY( Integer, newloopend )
OPT_KEY( Integer, newloopstart )
OPT_KEY( IntegerVector, norepackpositions )
OPT_KEY( Boolean, relaxbb )
OPT_KEY( Boolean, relaxjumps )
OPT_KEY( Integer, relaxrounds )
OPT_KEY( Real, repackradius )
OPT_KEY( Integer, stubchain )

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

void register_options() {
	utility::vector1 < core::Size > empty_vector;
	utility::vector1 < std::string > defaultaavector;
	defaultaavector.resize(18, "");
	defaultaavector[1]=("ALA");
	defaultaavector[2]=("ASP");
	defaultaavector[3]=("GLU");
	defaultaavector[4]=("PHE");
	defaultaavector[5]=("HIS");
	defaultaavector[6]=("ILE");
	defaultaavector[7]=("LYS");
	defaultaavector[8]=("LEU");
	defaultaavector[9]=("MET");
	defaultaavector[10]=("ASN");
	defaultaavector[11]=("PRO");
	defaultaavector[12]=("GLN");
	defaultaavector[13]=("ARG");
	defaultaavector[14]=("SER");
	defaultaavector[15]=("THR");
	defaultaavector[16]=("VAL");
	defaultaavector[17]=("TRP");
	defaultaavector[18]=("TYR");
	utility::vector1 < std::string > emptystringvector;

	NEW_OPT (acceptfirstmove, "Should I accept the first move by default?  True by default.  Use false if the loop or peptide is already in a good conformation, and you just want to refine it further.", true);
	NEW_OPT (allowedLaas, "List of allowed L-amino acids.  This should be entered as case-sensitive 3-letter codes separated by spaces.  By default, it includes all of the standard amino acids except CYS and GLY.  Use this flag with NONE (\"-allowedLaas NONE\") to disable L-amino acids.", defaultaavector);
	NEW_OPT (allowedDaas, "List of allowed D-amino acids.  This should be entered as case-sensitive 3-letter codes of the equivalent L-amino acids (e.g. ALA, SER, etc.) separated by spaces.  By default, this is an empty list (which disables D-amino acids).", emptystringvector);
	NEW_OPT (cyclic, "Are we designing a cyclic peptide?  False by default.", false);
	NEW_OPT (dab_sev_strength, "If specified, the dab_sev (DAlphaBall_SurfaceVolumeEnergy) score term is used during design and fast relaxation.  A final round of minimization without dab_sev is also applied.  Default unused.", 0.0);
	NEW_OPT (design_farep_mult, "Factor by which to multiply fa_rep during design steps.  Turning this down can help to improve packing.  If not specified, not used (full fa_rep strength used).", 1.0);
	NEW_OPT (design_faatr_mult, "Factor by which to multiply fa_atr during design steps.  Turning this up can help to improve packing.  If not specified, not used (normal fa_atr strength used).", 1.0);
	NEW_OPT (enumerate_only_mode, "If true, then backbone conformations for a chain of L-alanines will be generated.  No design or Monte Carlo search will be attempted; subsequent rounds will discard the previous solution.  Default false.", false);
	NEW_OPT (jittersize, "The average magnitude, in degrees, of the perturbation vector used to perturb the new loop backbone prior to relaxation.  If not specified, no perturbation is applied.", 0.0);
	NEW_OPT (looplength, "The length of the new loop, in amino acid residues.  Must be more than 2.  Default 5.", 5);
	NEW_OPT (maxattempts, "Maximum number of times the program will attempt to perturb start and end residues of a loop closure segment in order to find a clash-free conformation.  Default 100.", 100);
	NEW_OPT (maxiterations, "The number of iterations of kinematic closure and design that will be attempted.  Default 100.", 100);
	NEW_OPT (maxdesignrounds, "The number of rounds of design, perturbation, and relaxation that will be attempted for each kinematic closure.  Default 3.", 3);
	NEW_OPT (MCtemp, "The temperature used for application of the Metropolis criterion to accept or reject relaxations or moves.  Default 1.0.", 1.0);
	NEW_OPT (newloopend, "The residue that will serve as the C-terminal end of the new loop.  This must be in the chain specified with the -stubchain flag.  Mandatory input.", 0);
	NEW_OPT (newloopstart, "The residue that will serve as the N-terminal end of the new loop.  This must be in the chain specified with the -stubchain flag.  Mandatory input.", 0);
	NEW_OPT (norepackpositions, "List of positions that should not be repacked.  Default empty list.", empty_vector);
	NEW_OPT (relaxbb, "Should target backbone be relaxed?  Default true.  (Note that this only applies to the target.  The kinematically-closed chain is always backbone-relaxed.)", true);
	NEW_OPT (relaxjumps, "Should jumps be relaxed?  Default true.", true);
	NEW_OPT (relaxrounds, "Number of rounds of fast relax to apply.  Default 3.  Alternatively, a script file can be specified with the -relax:script flag.", 3);
	NEW_OPT (repackradius, "Radius, in Angstroms, around the loop that target residues are repacked.  Residues specified by the -norepackpositions flag are always ignored.  All residues are repacked (except for those listed with the -norepackpositions flag) if this is flag is not used.", 0.0);
	NEW_OPT (stubchain, "The chain that represents the peptide stub that will be extended to generate a cyclic peptide.  Default chain 2.", 2);
}

//Function to check the input flags for ridiculous user inputs.
bool do_initial_checks ()
{
	//TODO -- add check of allowedLaas and allowedDaas lists to ensure that everything provided is either a standard L-amino acid or a standard D-amino acid
	bool checksfailed = false;

	if(option[enumerate_only_mode]()) {
		printf("The -enumerate_only_mode flag was used.  Backbone conformations for a chain of L-alanines will be generated.  Note that no Monte Carlo search or design will be attempted.  Each round will discard the previous round's solution.\n");
		if(option[MCtemp].user()) {
			printf("Warning: the -MCtemp flag was used with the -enumerate_only_mode flag.  The value of the -MCtemp flag will be ignored, since no Monte Carlo search will be attempted.\n");
		}
	} else {
		if(option[MCtemp]() < 1.0e-12) {
			printf("Error!  The value specified with the -MCtemp flag cannot be negative!\n");
			checksfailed = true;
		} else {
			printf("Using a temperature value of %.4f for the Metropolis criterion (accepting or rejecting moves).\n", option[MCtemp]());
		}
	}

	if(option[OptionKeys::relax::script].user()) {
		printf("Using FastRelax script file %s.\n", option[OptionKeys::relax::script]().name().c_str());
	}

	if(option[jittersize].user()) {
		if(option[jittersize]() < 1.0e-12) {
			printf("Error!  The value specified with the -jittersize flag cannot be negative!\n");
			checksfailed = true;
		} else {
			printf("A small backbone perturbation (of average perturbation vector magnitude %.4f degrees) will be applied prior to relaxation (-jittersize flag).\n", option[jittersize]());
		}
	}

	if(option[dab_sev_strength].user()) {
		printf("Design and FastRelax steps will use a dab_sev (DAlphaBall_SurfaceVolumeEnergy) value of %.4f (-dab_sev_strength flag).  A final round of minimization without dab_sev will also be applied.", option[dab_sev_strength]());
	}

	if(option[design_farep_mult].user()) {
		if(option[enumerate_only_mode]()) {
			printf("Warning: the -design_farep_mult flag was used with the -enumerate_only_mode flag.  The value of -design_farep_mult will be ignored, since no design will be attempted.\n");
		} else {
			if(option[design_farep_mult]() < 1.0e-12) {
				printf("Error!  The value specified with the -design_farep_mult flag cannot be negative!\n");
				checksfailed = true;
			} else {
				printf("The fa_rep (repulsive) score term will be multiplied by a factor of %.2f when designing (-design_farep_mult flag).\n", option[design_farep_mult]());
			}
		}
	}

	if(option[design_faatr_mult].user()) {
		if(option[enumerate_only_mode]()) {
			printf("Warning: the -design_faatr_mult flag was used with the -enumerate_only_mode flag.  The value of -design_faatr_mult will be ignored, since no design will be attempted.\n");
		} else {
			if(option[design_faatr_mult]() < 1.0e-12) {
				printf("Error!  The value specified with the -design_faatr_mult flag cannot be negative!\n");
				checksfailed = true;
			} else {
				printf("The fa_atr (attractive) score term will be multiplied by a factor of %.2f when designing (-design_faatr_mult flag).\n", option[design_faatr_mult]());
			}
		}
	}


	if(option[maxattempts]()<1) {
			printf("Error!  The number of attempts at perturbing start and end residues (-maxattempts flag) must be at least 1.\n");
			checksfailed = true;
	}

	if(option[maxiterations]() < 1) {
			printf("Error!  The number of iterations (-maxiterations flag) must be at least 1.\n");
			checksfailed = true;
	}

	if(option[maxdesignrounds]() < 1) {
			printf("Error!  The number of design rounds (-maxdesignrounds flag) must be at least 1.\n");
			checksfailed = true;
	}

	if(!option[newloopstart].user()) {
			printf("Error!  The user must specify the loop start residue with the -newloopstart flag.\n");
			checksfailed = true;
	} else {
		if(option[newloopstart]() < 1) {
			printf("Error!  The user must provide a value greater than zero with the -newloopstart flag.\n");
			checksfailed = true;
		}
	}

	if(!option[newloopend].user()) {
			printf("Error!  The user must specify the loop end residue with the -newloopend flag.\n");
			checksfailed = true;
	} else {
		if(option[newloopend]() < 1) {
			printf("Error!  The user must provide a value greater than zero with the -newloopend flag.\n");
			checksfailed = true;
		}
	}

	if(option[cyclic]()) {
		if( option[newloopend]() > option[newloopstart]() ) {
			printf("Error!  For a cyclic peptide, the end of the loop must be a lower-numbered residue than the start.\n");
			checksfailed = true;
		}
	} else {
		if( option[newloopend]() < option[newloopstart]() ) {
			printf("Error!  The end of the loop must be a higher-numbered residue than the start.\n");
			checksfailed = true;
		}
	}

	if(option[stubchain]() < 1) {
		printf("Error!  The -stubchain flag cannot take a value less than 1.\n");
		checksfailed = true;
	} else printf("Chain %i is the peptide stub that will be extended to generate a cyclic peptide.\n", option[stubchain]());

	if(!option[in::file::s].user() || option[in::file::s]().size()!=1) {
		printf("Error!  The user must specify a single PDB file with the -in:file:s flag.\n");
		checksfailed = true;
	}

	if(option[relaxrounds]()<1) {
		printf("Error!  The number of rounds of relaxation (-relaxrounds flag) cannot be less than 1.\n");
		checksfailed = true;
	}

	if(option[norepackpositions].user() && option[norepackpositions]().size()<1) {
		printf("Error!  If the -norepackpositions flag is used, the user must provide a list of positions that cannot repack.\n");
		checksfailed = true;
	}

	if(option[repackradius].user()) {
		if(option[repackradius]() < 1e-12) {
			printf("Error!  If the -repackradius flag is used, the user must provide a repack radius greater than zero.\n");
			checksfailed = true;
		} else {
			printf("A repack radius of %.1f Angstroms around the redesigned loop or peptide will be used (-repackradius flag).\n", option[repackradius]());
		}
	}

	return checksfailed;
}

//Function to do some checks after PDB load to make sure that the user inputs are consistent with the input structure:
bool do_postload_checks( const core::pose::Pose &mypose )
{
	bool checksfailed = false;
	core::Size loopstart = (core::Size)option[newloopstart]();
	core::Size loopend = (core::Size)option[newloopend]();
	core::Size stub_chain = (core::Size)option[stubchain]();

	if((core::Size)mypose.chain(loopstart) != stub_chain) {
		printf("Error!  The starting residue for the loop (-newloopstart flag) is not in the stub chain (-stubchain flag).\n");
		checksfailed = true;
	}

	if((core::Size)mypose.chain(loopend) != stub_chain) {
		printf("Error!  The ending residue for the loop (-newloopend flag) is not in the stub chain (-stubchain flag).\n");
		checksfailed = true;
	}

	return checksfailed;
}

//Function to check whether a value is in a list:
bool is_in_list (const core::Size val, const utility::vector1 < core::Size > &vallist)
{
	if (vallist.size()==0) return false;
	for(core::Size i=1, listsize=vallist.size(); i<=listsize; i++) {
		if(vallist[i]==val) return true;
	}
	return false;
}

//Function to check whether a string is in a list of strings:
bool is_in_list (const std::string str, const utility::vector1 < std::string > &stringlist)
{
	if (stringlist.size()==0) return false;
	for(core::Size i=1, listsize=stringlist.size(); i<=listsize; i++) {
		if(stringlist[i]==str) return true;
	}
	return false;
}


void reset_sfxn (
	core::scoring::ScoreFunctionOP sfxn_constrained
) {
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	sfxn_constrained->set_weight(atom_pair_constraint, 0.0);
	sfxn_constrained->set_weight(dihedral_constraint, 0.0);
	sfxn_constrained->set_weight(angle_constraint, 0.0);
}

void constrain_loop_ends (
	core::pose::Pose &mypose,
	const core::Size resC, //The last residue of the loop (provides C)
	const core::Size resN, //The anchor residue (provides N)
	core::scoring::ScoreFunctionOP sfxn_constrained
) {
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::id;
	mypose.remove_constraints();
	sfxn_constrained->set_weight(atom_pair_constraint, 1.0);
	sfxn_constrained->set_weight(dihedral_constraint, 1.0);
	sfxn_constrained->set_weight(angle_constraint, 1.0);

	const std::string Hstring = (mypose.residue(resN).has("H")?"H":"CD");

	//Make the linkage between the N- and C-termini:
	mypose.conformation().declare_chemical_bond(resC, "C", resN, "N"); //Declare a chemical bond between the N and C termini.
	{
		//Peptide bond length constraint:
		mypose.add_constraint (
			new AtomPairConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ) ,
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new HarmonicFunc( 1.3288, 0.01)
			)
		);

		//Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(resC).atom_index("O") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index(Hstring) , resN) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);
		mypose.add_constraint (
			new DihedralConstraint (
				AtomID( mypose.residue(resC).atom_index("CA") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index("CA") , resN) ,
				new CircularHarmonicFunc( PI, 0.02)
			)
		);

		//Peptide bond angle constraints:
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index("CA") , resN) ,
				new CircularHarmonicFunc( CNCa_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				AtomID( mypose.residue(resN).atom_index(Hstring) , resN) ,
				new CircularHarmonicFunc( CNH_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("CA") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new CircularHarmonicFunc( CaCN_ANGLE/180.0*PI, 0.02)
			)
		);
		mypose.add_constraint (
			new AngleConstraint (
				AtomID( mypose.residue(resC).atom_index("O") , resC ),
				AtomID( mypose.residue(resC).atom_index("C") , resC ),
				AtomID( mypose.residue(resN).atom_index("N") , resN) ,
				new CircularHarmonicFunc( OCN_ANGLE/180.0*PI, 0.02)
			)
		);

	}

}

//Function to append a pose to another by a jump.
void append_pose_by_jump (
	core::pose::Pose &mypose,
	const core::pose::Pose &pose_to_append
) {
	if(pose_to_append.empty()) return;
	for(core::Size ir=1,nres=pose_to_append.n_residue(); ir<=nres; ir++) {
		if(ir==1) mypose.append_residue_by_jump( pose_to_append.residue(ir), 1, "", "", true );
		else mypose.append_residue_by_bond( pose_to_append.residue(ir), false);
	}
	return;
}

void append_alanines (
	core::pose::Pose &mypose, //This gets modified.
	const core::Size alacount, //Number of alanines to append.
	const bool append_D
) {
	core::pose::remove_lower_terminus_type_from_pose_residue(mypose, mypose.n_residue());

	std::string seq = "A";
	if (append_D) seq+="[DALA]";
	if(alacount>1) {
		for(core::Size i=2; i<=alacount; i++) {seq+="A"; if (append_D) seq+="[DALA]";}
	}

	core::chemical::ResidueTypeSetCAP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( seq, *standard_residues, false );

	for(core::Size i=1; i<=alacount; i++) {
		core::conformation::ResidueOP new_rsd = core::conformation::ResidueFactory::create_residue ( *requested_types[ i ] );
		if(!mypose.empty()) { //If the input pose is not empty, add the new residue as a bond.
			mypose.append_residue_by_bond( *new_rsd, true);
			mypose.set_psi(mypose.n_residue()-1, 135.0);
			mypose.set_omega(mypose.n_residue()-1, 180.0);
			mypose.set_phi(mypose.n_residue(), -135.0);
		} else { //If the input pose is empty, add the new residue as a jump.
			mypose.append_residue_by_jump(*new_rsd, 1, "", "", true);
		}
	}

	return;
}

//Function to add a series of alanines to the stub chain, replacing a loop if there is one there:
void create_loop (
	core::pose::Pose &stubpose, //This gets modified.
	const core::pose::Pose &masterpose, //For reference only.
	const core::Size stub_chain, //The number of the stub chain in the master pose.
	const bool d_only //Is the new loop only D-residues?
) {
	core::pose::Pose outpose;
	core::Size chaincount = 0;
	const core::Size startres = option[newloopstart]();
	const core::Size endres = option[newloopend]();
	const bool is_cyclic = option[cyclic]();
	bool first_residue_appended = false;

	for(core::Size ir=1, nres=masterpose.n_residue(); ir<=nres; ir++) {
		if(ir == 1 || (ir > 1 && masterpose.chain(ir)!=masterpose.chain(ir-1) ) ) chaincount++; //If this is the first residue of a new chain
		if(chaincount == stub_chain) { //If this is the stub chain
			if(is_cyclic) { //If we are designing a cyclic peptide, do not append anything until endres, then append everything up to startres, then append alanines.
				if(ir>=endres && ir<=startres) {
					if(!first_residue_appended) {
						outpose.append_residue_by_jump(masterpose.residue(ir), 1, "", "", true);
						core::pose::remove_lower_terminus_type_from_pose_residue(outpose, 1);
						first_residue_appended=true;
					} else {
						outpose.append_residue_by_bond(masterpose.residue(ir), false);
						core::pose::remove_upper_terminus_type_from_pose_residue(outpose, outpose.n_residue());
					}
				}
			} else { //if we're NOT designing a cyclic peptide, append everything up to startres, then append alanines, then do nothing up to endres, then append everything from endres to the end of the chain.
				if(ir<=startres) {
					if(!first_residue_appended) {
						outpose.append_residue_by_jump(masterpose.residue(ir), 1, "", "", true);
						first_residue_appended=true;
					} else {
						outpose.append_residue_by_bond(masterpose.residue(ir), false);
					}
					if(ir==startres) {
						core::pose::remove_upper_terminus_type_from_pose_residue(outpose, outpose.n_residue());
						append_alanines(outpose, option[looplength](), d_only);
					}
				} else if(ir>=endres) {
					outpose.append_residue_by_bond(masterpose.residue(ir), false);
				}
			}
		}

	}

	if(is_cyclic) append_alanines(outpose, option[looplength](), d_only);

	//Copy the pose:
	stubpose=outpose;

	return;
}

//Function to initialize the master pose.
void set_up_pose (
	core::pose::Pose &masterpose,
	utility::vector1 < core::Size > &no_repack_positions,
	core::Size &new_startpos, //The new starting position of the loop, after chain rearrangement (OUTPUT)
	core::Size &new_endpos, //The new ending position of the loop, after chain rearrangement (OUTPUT)
	const bool d_only //Is the new loop going to be only D-amino acids?
)
{
	const std::string infile = option[in::file::s]()[1];
	const core::Size stub_chain = (core::Size)option[stubchain]();
	const bool is_cyclic = option[cyclic]();
	const core::Size endpos = option[newloopend]();
	const core::Size startpos = option[newloopstart]();

	printf("Importing pose from %s\n", infile.c_str()); fflush(stdout);
	core::import_pose::pose_from_file(masterpose, infile, core::import_pose::PDB_file);

	utility::vector1 < core::pose::PoseOP > posechains = masterpose.split_by_chain();
	core::pose::Pose newpose;

	//New vector of positions that don't repack (since the order of the chain will change when the stub is added at the end).
	utility::vector1 < core::Size > new_no_repack_positions;
	//Temporarily add startpos and endpos to the no_repack_positions vector to conveniently update these values as chains are moved around.  (These will be stripped off again).
	no_repack_positions.push_back(startpos);
	no_repack_positions.push_back(endpos);
	const core::Size aas_in_stubchain = posechains[stub_chain]->n_residue();
	core::Size first_stubchain_aa=1;
	if(stub_chain > 1) {
		for(core::Size i=1; i<stub_chain; i++) first_stubchain_aa+=posechains[i]->n_residue();
	}
	core::Size aas_after_stubchain = 0;
	if(stub_chain < posechains.size()) { //If the stub chain isn't the last pose chain
		for(core::Size i=stub_chain, imax=posechains.size(); i<=imax; i++) aas_after_stubchain+=posechains[i]->n_residue();
	}
	if(no_repack_positions.size()>0) {
		for(core::Size i=1, imax=no_repack_positions.size(); i<=imax; i++) {
			if((core::Size)masterpose.chain(no_repack_positions[i]) < stub_chain) new_no_repack_positions.push_back(no_repack_positions[i]);
			else if((core::Size)masterpose.chain(no_repack_positions[i]) == stub_chain) {
				if(is_cyclic){//Cyclic case: chain is from endpos to startpos, with anything before endpos in the chain cut off.  It is shifted to the end of the pose.
					core::Size trimmed_residues = endpos-first_stubchain_aa;
					new_no_repack_positions.push_back( no_repack_positions[i] - trimmed_residues + aas_after_stubchain );
				} else { //Not cyclic case: chain is from first residue to startpos, then polyALA is appended, then from endpos to end residue.  Anything between startpos and endpos is cut out.  It is shifted to the end of the pose.
					if(no_repack_positions[i]<=startpos) new_no_repack_positions.push_back(no_repack_positions[i]+aas_after_stubchain);
					else if (no_repack_positions[i]>=endpos) {
						new_no_repack_positions.push_back( no_repack_positions[i] + aas_after_stubchain - (endpos-startpos-1) + (core::Size)option[looplength]() );
					}
				}
			}
			else if((core::Size)masterpose.chain(no_repack_positions[i]) >= stub_chain) new_no_repack_positions.push_back(no_repack_positions[i] - aas_in_stubchain);
		}
	}
	//for(core::Size i=1; i<=no_repack_positions.size(); i++) printf("%lu ", no_repack_positions[i]); printf("\n"); //DELETE ME
	//for(core::Size i=1; i<=new_no_repack_positions.size(); i++) printf("%lu ", new_no_repack_positions[i]); printf("\n"); fflush(stdout); //DELETE ME

	//Copy the updated start and end positions:
	new_startpos=new_no_repack_positions[new_no_repack_positions.size()-1];
	new_endpos=new_no_repack_positions[new_no_repack_positions.size()];
	new_no_repack_positions.resize(new_no_repack_positions.size()-2); //Delete the last two entries, which were the updated startpos and endpos
	no_repack_positions = new_no_repack_positions; //Copy the new vector

	create_loop(*(posechains[stub_chain]), masterpose, stub_chain, d_only); //Add a series of alanines to the stub chain, replacing a loop if there is one there.

	//Ensure that all the chains have proper termini, EXCEPT the stub chain:
	printf("Setting up termini.\n"); fflush(stdout);
	for(core::Size ichain=1, nchain=posechains.size(); ichain<=nchain; ichain++) {
		if((posechains[ichain]->residue(1).is_protein()) ) { //Skip non-protein chains
			if(ichain!=stub_chain) {
				core::pose::add_lower_terminus_type_to_pose_residue((*(posechains[ichain])), 1);
				core::pose::add_upper_terminus_type_to_pose_residue((*(posechains[ichain])), posechains[ichain]->n_residue());
			}
		}
		if(ichain!=stub_chain) append_pose_by_jump(newpose, *(posechains[ichain])); //Append the chain, UNLESS it is the stub chain.
	}

	append_pose_by_jump(newpose, *(posechains[stub_chain])); //Append the stub chain LAST.

	masterpose=newpose;
	masterpose.update_residue_neighbors();
	masterpose.conformation().detect_disulfides();

	printf("Import complete.\n"); fflush(stdout);

	return;
}

//Function to set which sidechains can repack, based on a user-input option:
void set_sc_repacking (core::kinematics::MoveMapOP mm, utility::vector1 < core::Size > &no_repack_positions)
{
	for(core::Size i=1, imax=no_repack_positions.size(); i<=imax; i++) {
		if(no_repack_positions[i] > 0) mm->set_chi(no_repack_positions[i], false);
	}
	return;
}

//Function to copy the position of a residue.:
void copy_residue_position (
	const core::pose::Pose &refpose,
	const core::Size refposition,
	core::pose::Pose &modifiedpose,
	const core::Size modifiedposition
) {
	/*const core::conformation::Residue &refres(refpose.residue(refposition)); //The reference residue.
	core::conformation::Residue modifiedres(modifiedpose.conformation().residue(modifiedposition)); //A copy of the residue to modify.
	core::Size atcount = refres.natoms();
	for(core::Size ia=1; ia<=atcount; ia++) {
		std::string atname = refres.atom_name(ia);
		if(modifiedres.has(atname)) modifiedres.set_xyz(atname, refres.xyz(ia)); //Copy xyz coordinates of each atom.
	}
	modifiedpose.replace_residue(modifiedposition, modifiedres, false);*/
   modifiedpose.replace_residue(modifiedposition, refpose.residue(refposition), false);
	modifiedpose.update_residue_neighbors();
	return;
}

//Function to align two poses.
//Yanked from a private member function in SuperimposeMover.cc.
core::Real superimposebb(
        core::pose::Pose & mod_pose,
        core::pose::Pose const & ref_pose,
        core::Size ref_start,
        core::Size ref_end,
        core::Size target_start,
        core::Size /*target_end*/ ,
	bool use_O
)
{
        core::id::AtomID_Map< core::id::AtomID > atom_map;
        std::map< core::id::AtomID, core::id::AtomID> atom_id_map;
        core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::BOGUS_ATOM_ID );
        for ( core::Size i_target = target_start, i_ref = ref_start; i_ref <= ref_end; ++i_ref, ++i_target ) {

                //if ( ! mod_pose.residue(i_target).has("N") ) continue;
                //if ( ! ref_pose.residue(i_ref).has("N") ) continue;
                core::id::AtomID const id1( mod_pose.residue(i_target).atom_index("N"), i_target );
                core::id::AtomID const id2( ref_pose.residue(i_ref).atom_index("N"), i_ref );
                atom_map.set( id1, id2 );
                atom_id_map.insert( std::make_pair(id1, id2) );

                //if ( ! mod_pose.residue(i_target).has("CA") ) continue;
                //if ( ! ref_pose.residue(i_ref).has("CA") ) continue;
                core::id::AtomID const id3( mod_pose.residue(i_target).atom_index("CA"), i_target );
                core::id::AtomID const id4( ref_pose.residue(i_ref).atom_index("CA"), i_ref );
                atom_map.set( id3, id4 );
                atom_id_map.insert( std::make_pair(id3, id4) );

                //if ( ! mod_pose.residue(i_target).has("C") ) continue;
                //if ( ! ref_pose.residue(i_ref).has("C") ) continue;
                core::id::AtomID const id5( mod_pose.residue(i_target).atom_index("C"), i_target );
                core::id::AtomID const id6( ref_pose.residue(i_ref).atom_index("C"), i_ref );
                atom_map.set( id5, id6 );
                atom_id_map.insert( std::make_pair(id5, id6) );

                //if ( ! mod_pose.residue(i_target).has("HA") ) continue;
                //if ( ! ref_pose.residue(i_ref).has("HA") ) continue;
		if(use_O && mod_pose.residue(i_target).has("O") && ref_pose.residue(i_ref).has("O")) {
			core::id::AtomID const id7( mod_pose.residue(i_target).atom_index("O"), i_target );
			core::id::AtomID const id8( ref_pose.residue(i_ref).atom_index("O"), i_ref );
			atom_map.set( id7, id8 );
			atom_id_map.insert( std::make_pair(id7, id8) );
		} else {
			if(mod_pose.residue(i_target).has("CB") && ref_pose.residue(i_ref).has("CB")) {
				core::id::AtomID const id7( mod_pose.residue(i_target).atom_index("CB"), i_target );
				core::id::AtomID const id8( ref_pose.residue(i_ref).atom_index("CB"), i_ref );
				atom_map.set( id7, id8 );
			     	atom_id_map.insert( std::make_pair(id7, id8) );
			} else { //gly-to-gly alignment:
				if(mod_pose.residue(i_target).has("1HA") && ref_pose.residue(i_ref).has("1HA")) {
					core::id::AtomID const id7( mod_pose.residue(i_target).atom_index("1HA"), i_target );
					core::id::AtomID const id8( ref_pose.residue(i_ref).atom_index("1HA"), i_ref );
					atom_map.set( id7, id8 );
					atom_id_map.insert( std::make_pair(id7, id8) );
				} else {
					printf("Error!  Could not align.\n"); fflush(stdout); exit(1);
				}
			}
		}
        }
        return core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}

//Function to do quick backbone clash checks for a residue against all residues that are not in the same chain.
bool noclash (
	const core::Size position,
	const core::pose::Pose &mypose
) {
	core::Size curchain = (core::Size)mypose.chain(position);
	utility::vector1 < std::string > atlist;
	atlist.push_back("N");
	atlist.push_back("CA");
	atlist.push_back("C");
	atlist.push_back("O");
	atlist.push_back("CB");

	for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++) {
		if((core::Size)mypose.chain(ir)==curchain) continue; //No clash check within chain.

		for(core::Size ia1=1, ia1max=atlist.size(); ia1<=ia1max; ia1++) { //Loop through atoms in this residue
			if(!mypose.residue(ir).has(atlist[ia1])) continue; //Skip if this residue doesn't have this atom (e.g. CB).

			for(core::Size ia2=1, ia2max=atlist.size(); ia2<=ia2max; ia2++) { //Loop through atoms in the residue that we want to compare.
				if(!mypose.residue(position).has(atlist[ia2])) continue; //Skip if the comparison residue doesn't have this atom (e.g. CB).

				core::Real dist_cutoff_squared = mypose.residue(position).atom_type( mypose.residue(position).atom_index( atlist[ia2] ) ).lj_radius();
				dist_cutoff_squared += mypose.residue(ir).atom_type( mypose.residue(ir).atom_index( atlist[ia1] ) ).lj_radius();
				dist_cutoff_squared = dist_cutoff_squared*dist_cutoff_squared;

				core::Real dist_squared = (mypose.residue(position).xyz( atlist[ia2] ) - mypose.residue(ir).xyz( atlist[ia1] )).length_squared();
				if(dist_squared<dist_cutoff_squared) return false;

			}
		}
	}

	printf("Mini clash check for residue %lu passed.\n", position); fflush(stdout);

	return true;
}

//Function to check that all peptide bonds are sensible in a stretch:
bool good_geometry (
	const core::pose::Pose &mypose,
	const core::Size start,
	const core::Size end,
	const core::Size additional
) {
	numeric::xyzVector < core::Real > Nxyz;
	numeric::xyzVector < core::Real > Cxyz;
	core::Real dist;

	for(core::Size ir=start; ir<=end; ir++) {
		Cxyz = mypose.residue(ir).xyz("C");
		if(ir<end) Nxyz = mypose.residue(ir+1).xyz("N");
		else Nxyz = mypose.residue(additional).xyz("N");
		dist = (Cxyz-Nxyz).length();
		if(dist > 1.4 || dist < 1.2) {
			printf("Peptide bond length between residues %lu and %lu is %.2f.\n", ir, ir+1, dist);  fflush(stdout);
			return false;
		}
	}
	return true;
}

//Function to generate new loop geometry by kinematic closure
void kinclose(
	core::pose::Pose &mypose,
	const core::Size startres, //The stub residue that starts the new loop
	const core::Size endres, //The stub residue than ends the new loop
	const core::Size loop_length,
	const bool is_cyclic,
	const bool close_random_subset,
	numeric::random::RandomGenerator &RG,
	core::scoring::ScoreFunctionOP sfxn,
	const bool allowL, const bool allowD, const bool allowGLY
) {
	printf("Attempting kinematic closure (startres=%lu, endres=%lu).\n", startres, endres); fflush(stdout); //DELETE ME

	bool sample_gly = allowGLY;
	if(allowL && allowD) sample_gly = false; //We don't need to sample gly explicitly if we're sampling both L- and D- amino acids.

	while(true) {

		bool perturbationfailed = false;

		core::Size start=startres+1; //First pivot.
		core::Size end=startres+loop_length; //Third pivot.  Note that this is measured from startres, since endres might be before startres (cyclic case).
		core::Size mid=RG.random_range(startres+2, startres+loop_length-1);//Second pivot.  Note that this is measured from startres, since endres might be before startres (cyclic case).
		core::Size beforestart=startres; //The residue before the first pivot.
		core::Size afterend = endres; //The residue after the last pivot.
		core::pose::Pose mypose_copy = mypose; //Make a copy of mypose.

	 	//Mutate mypose_copy to sample L-amino acids, D-amino acids, or glycine:
		for (core::Size ir=start; ir<=end; ir++) {
			std::string aaname = "ALA";
			core::Size dieroll = RG.random_range(1, 100); //Pick a number from 1 to 100.
			if(allowD && !allowL) aaname = "DALA"; //If we're only sampling D-amino acids, change this to D-alanine
			else if (allowD && allowL && dieroll <=50) aaname = "DALA"; //If we're sampling both, flip a coin and decide whether to change this to D-alanine

			if(sample_gly) {
				dieroll = RG.random_range(1,100); //We'll use a one in ten chance of sampling a gly, if we're sampling gly.
				if(dieroll <=10) aaname = "GLY";
			}

			//Mutate this position to L-ala, D-ala, or glycine, appropriately:
			protocols::simple_moves::MutateResidue mutres1(ir, "GLY"); //Mutate to glycine first, since otherwise I fear that the side-chain will be in the wrong place.
			mutres1.apply(mypose_copy);
			protocols::simple_moves::MutateResidue mutres2(ir, aaname);
			mutres2.apply(mypose_copy);
		}
		mypose_copy.update_residue_neighbors();
		//mypose_copy.conformation().detect_disulfides();

		printf("Trial sequence for kinematic closure: ");
		for(core::Size ir=start; ir<=end; ir++) printf("%s ", mypose_copy.residue(ir).name3().c_str());
		printf("\n"); fflush(stdout);

		if(close_random_subset){ //If we're closing a random subset of the new loop, pick random start, mid, and end points.
			mid=RG.random_range(startres+2, startres+loop_length-1);
			start=RG.random_range(startres+1, mid-1);
			end=RG.random_range(mid+1, startres+loop_length);
			beforestart=start-1;

			//The residue after the last residue OF THE LOOP is endres, which might come earlier in sequence (cyclic case).
			//In all other cases, the residue after the last residue of the segment that we're closing is just the i+1 residue:
			afterend=((end==startres+loop_length) ? endres : end+1);
		}

		printf("\tstart=%lu, mid=%lu, end=%lu, afterend=%lu, beforestart=%lu\n", start, mid, end, afterend, beforestart); fflush(stdout);

		{ //Scope 1: Perturb psi, omega, phi of the beforestart-start linkage.
			core::Size numattempts = 1;
			while(true) {
				mypose_copy.set_psi(beforestart, RG.uniform()*360.0-180.0);
				mypose_copy.set_omega(beforestart, 180.0);
				mypose_copy.set_phi(start, RG.uniform()*360.0-180.0);
				mypose_copy.update_residue_neighbors();
				if(noclash(start, mypose_copy) ) break; //Keep perturbing until no clash
				numattempts++;
				if(numattempts > (core::Size)option[maxattempts]()) {
					perturbationfailed = true;
					break;
				}
			}
		}

		if(perturbationfailed) {
			printf("Could not find clash-free perturbation for start residue (-maxattempts flag exceeded).  Trying again with new start, middle, and end positions.\n");
			continue;
		}

		{	//Scope 2: make a dipeptide for the last loop residue plus the end residue, randomly perturb, align to the end residue, and copy the last loop residue's position into mypose_copy.
			//If this is NOT a cyclic peptide, the last part of the last chain needs to be repositioned, too.
			core::pose::Pose dipeptide1;
			dipeptide1.append_residue_by_jump(mypose_copy.residue(end), 1, "", "", true);
			dipeptide1.append_residue_by_bond(mypose_copy.residue(afterend), true); //Append the end residue of the loop, building ideal geometry.
			dipeptide1.append_residue_by_bond(mypose_copy.residue(afterend), true); //Append the end residue of the loop again, building building ideal geometry, just to have another residue.

			//Set psi, omega, phi:
			core::Size numattempts = 1;
			while(true) {
				dipeptide1.set_psi(1, RG.uniform()*360.0-180.0);
				dipeptide1.set_omega(1, 180.0);
				dipeptide1.set_phi(2, RG.uniform()*360.0-180.0);
				dipeptide1.set_psi(2, mypose_copy.psi(afterend));
				dipeptide1.set_omega(2, mypose_copy.omega(afterend));

				superimposebb(dipeptide1, mypose, afterend, afterend, 2, 2, true); //Superimpose the afterend residue

				copy_residue_position(dipeptide1, 2, mypose_copy, afterend); //Copy the position of the afterend residue
				copy_residue_position(dipeptide1, 1, mypose_copy, end); //Copy the position of the end residue

				if(!is_cyclic && mypose_copy.n_residue()>afterend) {
					for(core::Size ir=afterend+1, nres=mypose_copy.n_residue(); ir<=nres; ir++) copy_residue_position(mypose, ir, mypose_copy, ir); //Copy the positions of the rest of the residues in the chain.
				}

				mypose_copy.update_residue_neighbors();

				if(noclash(end, mypose_copy)) break; //Keep perturbing until no clash
				numattempts++;
				if(numattempts > (core::Size)option[maxattempts]()) {
					perturbationfailed = true;
					break;
				}
			}

			if(is_cyclic && afterend == endres) {
				mypose_copy.append_residue_by_bond(dipeptide1.residue(2), false); //Append another copy of the end residue.
				mypose_copy.conformation().delete_residue_slow(afterend);
				//mypose_copy.conformation().detect_disulfides();
				beforestart--; start--; mid--; end--;
				//mypose_copy.dump_pdb("mypose_copy.pdb"); //DELETE ME
			}
		}

		if(perturbationfailed) {
			printf("Could not find clash-free perturbation for end residue (-maxattempts flag exceeded).  Trying again with new start, middle, and end positions.\n");
			continue;
		}

		//Kinematic closure:
		protocols::loops::loop_closure::kinematic_closure::KinematicMover kinmover;
		kinmover.set_temperature( 1.0 );
		kinmover.set_vary_bondangles( false );
		kinmover.set_sample_nonpivot_torsions( true ); //true
		kinmover.set_rama_check( true );
		kinmover.set_idealize_loop_first( true ); //true
		kinmover.set_sfxn(sfxn);
		kinmover.set_pivots(start, mid, end);

		kinmover.apply(mypose_copy);

		if(	kinmover.last_move_succeeded() &&
				good_geometry(mypose_copy,
					startres - ((is_cyclic && afterend == endres)?1:0),
					startres+loop_length,
					afterend
				)
		  )
		{
			//mypose.dump_pdb("premut.pdb"); //DELETE ME
			for(core::Size ir=beforestart; ir<=end; ir++) {
				//Mutate the mypose sequence to match mypose_copy:
				std::string rsdname = mypose_copy.residue(ir).name();
				protocols::simple_moves::MutateResidue mutres(ir+( (is_cyclic && afterend == endres) ? 1 : 0 ), rsdname.substr(0, rsdname.find(chemical::patch_linker)));
				mutres.apply(mypose);
				copy_residue_position(mypose_copy, ir, mypose, ir+( (is_cyclic && afterend == endres) ? 1 : 0 ));
			}
			copy_residue_position(mypose_copy, ( (is_cyclic && afterend == endres) ? mypose_copy.n_residue() : afterend), mypose, afterend);
			mypose.update_residue_neighbors();
			//mypose.dump_pdb("postmutmove.pdb"); //DELETE ME
			break;
		} else {
			//mypose_copy.dump_pdb("temp.pdb");  exit(1); //DELETE ME
			printf("Closure failed.  Randomizing and trying again.\n");
		}
	}

	return;
}

//Function to determine whether a residue is near the design chain:
bool near_chain (
	const core::pose::Pose &mypose,
	const core::Size position,
	const core::Size designchain,
	const core::Real radius
) {
	const numeric::xyzVector < core::Real > resCAxyz = mypose.residue(position).atom("CA").xyz();
	numeric::xyzVector < core::Real > chainCAxyz;
	core::Real distsq=0.0;
	const core::Real radiussq=radius*radius;

	for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++) {
		if((core::Size)mypose.chain(ir)!=designchain) continue; //Skip residues outside of the design chain
		if(!mypose.residue(ir).has("CA")) continue; //Skip if there's no alpha carbon
		chainCAxyz = mypose.residue(ir).atom("CA").xyz();
		distsq = (chainCAxyz - resCAxyz).length_squared();
		if(distsq < radiussq) return true;
	}

	return false;
}

//Function to design the new loop
void designloop(
	core::pose::Pose &mypose,
	const core::Size startres,
	const core::Size endres,
	const core::Size loop_length,
	const bool is_cyclic,
	core::scoring::ScoreFunctionOP sfxn,
	const utility::vector1 < core::Size > &no_repack_positions,
	utility::vector1 < core::Size > &can_repack_positions,
	const utility::vector1 < std::string > &allowed_L_aas,
	const utility::vector1 < std::string > &allowed_D_aas,
	const bool allow_gly
) {
	printf("Designing loop.\n"); fflush(stdout);

	core::Real farep_default = sfxn->get_weight(core::scoring::fa_rep);
	core::Real faatr_default = sfxn->get_weight(core::scoring::fa_atr);

	if(option[design_farep_mult].user()) { //If the user has specified a factor by which to reduce fa_rep during design:
		sfxn->set_weight(core::scoring::fa_rep, farep_default*option[design_farep_mult]());
	}

	if(option[design_faatr_mult].user()) { //If the user has specified a factor by which to reduce fa_rep during design:
		sfxn->set_weight(core::scoring::fa_rep, faatr_default*option[design_faatr_mult]());
	}


	const core::Size designchain = mypose.chain(mypose.n_residue()); //The last chain is the design chain

	utility::vector1 < core::Size > dpositions; //List of D-amino acids

	//mypose.dump_pdb("temp.pdb"); //DELETE ME

	for(core::Size ir=startres,lastres=(is_cyclic?mypose.n_residue():endres); ir<=lastres; ir++) {
		if(core::chemical::is_canonical_D_aa(mypose.residue(ir).aa())) {
			dpositions.push_back(ir);
			if(!is_in_list(mypose.residue(ir).name3(), allowed_D_aas)) { //If this is DALA and DALA isn't an allowed D-amino acid, mutate to glycine temporarily
				protocols::simple_moves::MutateResidue mutres(ir, "GLY");
				mutres.apply(mypose);
			}
		}
	}

	//DELETE THE FOLLOWING:
	//if(dpositions.size()>0){
	//	printf("DPOSITIONS: ");
	//	for(core::Size i=1, imax=dpositions.size(); i<=imax; i++) printf("%lu ", dpositions[i]);
	//	printf("\n"); fflush(stdout);
	//}

	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(mypose);

	for(core::Size ir=1,nres=mypose.n_residue(); ir<=nres; ir++) {
		if(!mypose.residue(ir).is_protein()) {
			task->nonconst_residue_task(ir).prevent_repacking();
			continue;
		}
		if((core::Size)mypose.chain(ir)==designchain && ir>startres && (is_cyclic?true:ir<endres)) { //If this is the design chain
			if(is_in_list(ir, dpositions)) { //If this is a D-amino acid position
				utility::vector1< bool > aas(20,false);
				if(allow_gly) aas[core::chemical::aa_gly] = true;
				task->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas);

				if(allowed_D_aas.size()>0) {
					for(core::Size i=1, imax = allowed_D_aas.size(); i<=imax; i++) { //Loop through the allowed D-amino acids
						task->nonconst_residue_task(ir).allow_aa( core::chemical::get_D_equivalent( core::chemical::aa_from_name( allowed_D_aas[i] )  ) );
					}
				}
			} else {
				utility::vector1< bool > aas(20,false);
				if(allow_gly) aas[core::chemical::aa_gly] = true;
				if(allowed_L_aas.size()>0 && !is_in_list("NONE", allowed_L_aas)) {
					for(core::Size i=1, imax=allowed_L_aas.size(); i<=imax; i++) {
						aas[ core::chemical::aa_from_name( allowed_L_aas[i] ) ] = true;
					}
				}

				task->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas);
			}
		} else { //If this is NOT the design chain (or is the design chain, but outside of the loop to be designed)
			if( is_in_list(ir, no_repack_positions) ) {
				task->nonconst_residue_task(ir).prevent_repacking();
			} else if( is_in_list(ir, can_repack_positions) ) {
				task->nonconst_residue_task(ir).restrict_to_repacking();
			} else if( option[repackradius].user() ){
				if( near_chain(mypose, ir, designchain, option[repackradius]()) ) {
					can_repack_positions.push_back (ir);
					task->nonconst_residue_task(ir).restrict_to_repacking();
				} else {
					task->nonconst_residue_task(ir).prevent_repacking();
				}
			} else {
				task->nonconst_residue_task(ir).restrict_to_repacking();
			}
		}
	}

	//Repack, using the PackerTask declared above:
	protocols::simple_moves::PackRotamersMover repack( sfxn, task );
	repack.apply(mypose);

	printf("Design complete.  Sequence is ");
	for(core::Size ir=startres, irmax=startres+loop_length; ir<=irmax; ir++) {
		char aaname = mypose.residue(ir).name1();
		if( core::chemical::is_canonical_D_aa(mypose.residue(ir).aa()) ) aaname=tolower(aaname);
		printf("%c", aaname);
		if(ir==irmax) {
			aaname=mypose.residue(endres).name1();
			if( core::chemical::is_canonical_D_aa(mypose.residue(endres).aa()) ) aaname=tolower(aaname);
			printf("%c", aaname);
		}
	}
	printf("\n"); fflush(stdout);

	if(option[design_farep_mult].user()) { //If the user has specified a factor by which to reduce fa_rep during design, reset it now.
		sfxn->set_weight(core::scoring::fa_rep, farep_default);
	}

	if(option[design_faatr_mult].user()) { //If the user has specified a factor by which to reduce fa_rep during design, reset it now.
		sfxn->set_weight(core::scoring::fa_rep, faatr_default);
	}


	return;
}

//Function to relax the new loop
void relaxloop(
	core::pose::Pose &mypose,
	const core::Size startres,
	const core::Size endres,
	const core::Size loop_length,
	protocols::relax::FastRelaxOP frlx,
	protocols::simple_moves::MinMoverOP minmove,
	core::kinematics::MoveMapOP mm,
	const utility::vector1 < core::Size > &no_repack_positions, //Positions whose sidechains can NOT repack or relax (takes precidence over can_repack positions -- entries in both lists can NOT repack or relax).
	const utility::vector1 < core::Size > &can_repack_positions //Positions whose sidechains CAN repack or relax.
) {
	printf("Attempting relax...\n"); fflush(stdout);

	mm->set_jump(option[relaxjumps]());
	mm->set_chi(option[repackradius].user()?false:true);
	mm->set_bb(option[relaxbb]());

	//Allow chi and bb minimization of loop positions:
	for (core::Size ir = startres, irmax = startres+loop_length; ir<=irmax; ir++) {
		if( !is_in_list( ir, no_repack_positions ) ) {
			mm->set_bb(ir, (ir==startres? option[relaxjumps]() :true) );
			mm->set_chi(ir, true);
		}
	}

	mm->set_bb(endres, option[relaxbb]()); //Redundant, but better to make this explicit
	mm->set_chi(endres, (is_in_list(endres, no_repack_positions)?false:true) );
	mm->set_bb(startres, option[relaxbb]()); //Redundant, but better to make this explicit
	mm->set_chi(startres, (is_in_list(endres, no_repack_positions)?false:true) );
	mm->set(core::id::TorsionID(startres, core::id::BB, 2), true); //Allow startres psi to be minimized
	mm->set(core::id::TorsionID(endres, core::id::BB, 1), false); //Do not allow endres phi to be minimized (TEMPORARILY -- TODO: allow jumps and fix this!)

	if(option[repackradius].user() && can_repack_positions.size()!=0) { //If the user has specified a relax radius, chi is false by default --> need to turn ON for entries in the can_repack_positions list.
		for(core::Size i=1, imax=can_repack_positions.size(); i<=imax; i++) {
			if( !is_in_list(can_repack_positions[i], no_repack_positions) ) mm->set_chi(can_repack_positions[i], true); //Turn on ONLY if not in the no_repack_positions list.
		}
	} else if(!option[repackradius].user()){ //If the user has NOT specified a relax radius, all residues can repack by default --> need to turn OFF for entries in the no_repack_positions list.
		for(core::Size ir=1, nres=mypose.n_residue(); ir<=nres; ir++) {
			if(is_in_list(ir, no_repack_positions)) mm->set_chi(ir, false);
		}
	}

	if(option[relaxjumps]()) {
		printf("Doing initial minimization without jumps.\n"); fflush(stdout);
		mm->set_jump(false);
		minmove->apply(mypose);
		mm->set_jump(true);
	}


	if(option[relaxjumps]()) {printf("Doing fastrelax with jumps.\n"); fflush(stdout);}
	frlx->apply(mypose);
	printf("Relax complete.\n"); fflush(stdout);

	return;
}

//Function to set the new loop conformation to an input conformation from a reference pose:
void set_loop_conformation(
	core::pose::Pose &mypose,
	const core::pose::Pose &refpose, //The reference pose (bestpose)
	const core::Size startres,
	const core::Size endres,
	const core::Size loop_length
) {
	core::Size irmax=startres+loop_length;
	const bool is_cyclic = (startres+loop_length != endres-1); //Is this a cyclic peptide?
	const core::Size jumpnumber = mypose.num_jump(); //The jump should be the last one

	if(option[relaxjumps]()) { //Set the jump if jumps have been relaxed:
		mypose.set_jump(jumpnumber, refpose.jump(jumpnumber));
	}

	if(!is_cyclic) irmax++; //If this is not a cyclic peptide, we want to set phi for the end residue

	for(core::Size ir=startres; ir<=irmax; ir++) {
		if (ir>startres) mypose.set_phi(ir, refpose.phi(ir));
		if(ir<irmax) {
			mypose.set_psi(ir, refpose.psi(ir));
			mypose.set_omega(ir, refpose.omega(ir));
		}
	}

	if(!is_cyclic) {
		//If this is not a cyclic peptide, all the backbone dihedral values are correct, now, BUT the residues after the loop
		//are hanging out in space somewhere.  We need to fix this now.
		core::pose::Pose afterresidues = mypose; //Make a copy of mypose.
		superimposebb(afterresidues, refpose, endres, afterresidues.n_residue(), endres, refpose.n_residue(), true); //Superimpose the whole last chain.
		for(core::Size ir=endres, irmax=afterresidues.n_residue(); ir<=irmax; ir++) copy_residue_position(afterresidues, ir, mypose, ir); //Copy the residue positions.
		return; //We're done if this is NOT a cyclic peptide.  Otherwise, proceed to set the last psi and after-end phi values.
	}

	core::Real lastpsi = numeric::dihedral_radians (
		refpose.residue(irmax).xyz("N"),
		refpose.residue(irmax).xyz("CA"),
		refpose.residue(irmax).xyz("C"),
		refpose.residue(irmax).xyz( "O")
	);
	mypose.conformation().set_torsion_angle(//Set psi of last residue
		core::id::AtomID(mypose.residue(irmax).atom_index("N"), irmax),
		core::id::AtomID(mypose.residue(irmax).atom_index("CA"), irmax),
		core::id::AtomID(mypose.residue(irmax).atom_index("C"), irmax),
		core::id::AtomID(mypose.residue(irmax).atom_index("O"), irmax),
		lastpsi);

	core::Real endphi = numeric::dihedral_radians (
		refpose.residue(endres).xyz( (refpose.residue(endres).has("H") ? "H" : "CD") ),
		refpose.residue(endres).xyz("N"),
		refpose.residue(endres).xyz("CA"),
		refpose.residue(endres).xyz( "C")
	);

	core::pose::Pose single_aa_pose;
	single_aa_pose.append_residue_by_jump(refpose.residue(endres),1,"","",true);
	single_aa_pose.conformation().set_torsion_angle(//Set phi of endres residue
		core::id::AtomID(single_aa_pose.residue(1).atom_index( (single_aa_pose.residue(1).has("H") ? "H" : "CD") ), 1),
		core::id::AtomID(single_aa_pose.residue(1).atom_index("N"), 1),
		core::id::AtomID(single_aa_pose.residue(1).atom_index("CA"), 1),
		core::id::AtomID(single_aa_pose.residue(1).atom_index("C"), 1),
		endphi);
	superimposebb(single_aa_pose, refpose, endres, endres, 1, 1, false); //Superimpose the residue following phi change.
	copy_residue_position(single_aa_pose, 1, mypose, endres); //Copy the phi-altered residue.  Blargh -- it shouldn't be this complicated.

	return;
}

//Function that applies the Metropolies criterion
bool metropolis_criterion (
	const core::Real newval,
	const core::Real oldval,
	const core::Real MCtemperature,
	numeric::random::RandomGenerator &RG
) {
	if( newval < oldval) return true;
	else {
		core::Real randval = RG.uniform();
		core::Real expval = exp( -1.0 * (newval-oldval) / MCtemperature );
		return (expval > randval);
	}

	return true;
}

//Function to add a small random value to all backbone dihedral angles:
void jitterloop (
	core::pose::Pose &mypose,
	const core::Size startres,
	const core::Size loop_length,
	const core::Real magnitude,
	numeric::random::RandomGenerator &RG
) {
	printf("Perturbing loop.\n"); fflush(stdout);
	core::Real accumulator = 0.0;
	utility::vector1 < core::Real > pertvect;
	pertvect.resize(loop_length*3, 0.0);
	for(core::Size i=1, imax=pertvect.size(); i<=imax; i++) {
		pertvect[i]=RG.gaussian();
		if(i%3==2) pertvect[i] *= 0.05; //Omega perturbations should be SMALL.
		accumulator+=(pertvect[i]*pertvect[i]);
	}
	accumulator=RG.gaussian()*magnitude/sqrt(accumulator);
	for(core::Size i=1, imax=pertvect.size(); i<=imax; i++) pertvect[i] *= accumulator;

	for(core::Size counter=1, ir=startres, nr=startres+loop_length; ir<=nr; ir++) {
		if(ir>startres) mypose.set_phi(ir, mypose.phi(ir) + pertvect[counter++]);
		if(ir<nr) {
			mypose.set_psi(ir, mypose.psi(ir) + pertvect[counter++]);
			mypose.set_omega(ir, mypose.omega(ir) + pertvect[counter++]);
		}
	}
	mypose.update_residue_neighbors();

	return;
}

static 	numeric::random::RandomGenerator RG( 847322 ); //Random generator and seed

//MAIN
int main( int argc, char * argv [] ) {

	using namespace protocols::moves;
	using namespace utility::file;
	using namespace std;

	printf("Starting design_truecycpeptide_fromstub.cc\nFile created 21 Oct 2013 by Vikram K. Mulligan.\n\n"); fflush(stdout);


	register_options();
	devel::init(argc, argv);

	//Do initial checks of the options provided by the user:
	bool failnow = do_initial_checks();
	if(failnow) { printf("Program will now terminate.\n"); fflush(stdout); exit(1); }

	//Scorefunction, movemap, and reusable movers:
	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::get_score_function();
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_jump(option[relaxjumps]());
	mm->set_bb(option[relaxbb]());
	mm->set_chi(true);
	// AMW: VKM prefers dfpmin here
	protocols::simple_moves::MinMoverOP minmove = new protocols::simple_moves::MinMover(mm, sfxn, "dfpmin_armijo_nonmonotone", 0.0000001, true, false, false);
	protocols::relax::FastRelaxOP frlx;
	if(option[OptionKeys::relax::script].user()) {
		frlx = new protocols::relax::FastRelax(sfxn, option[OptionKeys::relax::script]()); //A fastrelax mover initialized from a script
	} else {
		frlx = new protocols::relax::FastRelax(sfxn, (core::Size) option[relaxrounds]()); //A fastrelax mover
	}
	frlx->set_movemap(mm);

	//Allowed amino acids:
	utility::vector1 < std::string > allowed_L_aas = option[allowedLaas]();
	utility::vector1 < std::string > allowed_D_aas = option[allowedDaas]();
	const bool allowL = (allowed_L_aas.size()>0 && !is_in_list("NONE", allowed_L_aas));
	const bool allowD = (allowed_D_aas.size()>0);
	const bool allowGLY = ( is_in_list("GLY", allowed_L_aas) || is_in_list("GLY", allowed_D_aas) ); //Is glycine in one of the allowed amino acid lists?  If the other list is empty, the kinclose function will allow glys.
	if(allowL) {
		printf ("These L-amino acids are allowed in the designed loop/cyclic peptide:\t");
		for(core::Size i=1, imax=allowed_L_aas.size(); i<=imax; i++) { printf("%s ", allowed_L_aas[i].c_str());}
		printf("\n"); fflush(stdout);
	}
	if(allowD) {
		printf ("These D-amino acids are allowed in the designed loop/cyclic peptide:\t");
		for(core::Size i=1, imax=allowed_D_aas.size(); i<=imax; i++) { printf("%s ", allowed_D_aas[i].c_str());}
		printf("\n"); fflush(stdout);
	}

	//The master pose and positions that don't repack:
	core::pose::Pose masterpose;
	utility::vector1 < core::Size > no_repack_positions = option[norepackpositions]();
	core::Size startres=(core::Size)option[newloopstart](); //This will be updated, below, to point to the new startres position (given that the start_chain has moved).
	core::Size endres=(core::Size)option[newloopend](); //This will be updated, below, to point to the new endres position (given that the start_chain has moved, and that a loop has been inserted).
	set_up_pose(masterpose, no_repack_positions, startres, endres, !allowL); //Import and rearrange chains, updating the positions that can't repack and the starting and ending points of the chain0.
	if(option[norepackpositions].user()) set_sc_repacking(mm, no_repack_positions); //Set the sidechain positions that can't repack.

	//Do checks of options provided by the user that are dependent on the pose:
	failnow = do_postload_checks(masterpose);
	if(failnow) { printf("Program will now terminate.\n"); fflush(stdout); exit(1); }

	(*sfxn)(masterpose);
	//masterpose.dump_scored_pdb("temp.pdb", *sfxn); //DELETE ME

	//Containers for poses:
	core::pose::Pose trialpose, trialpose2;
	core::pose::Pose bestpose=masterpose;
	core::pose::Pose verybestpose=masterpose;

	//Input variables used henceforth:
	const core::Size loop_length=(core::Size)option[looplength]();
	const bool is_cyclic=option[cyclic]();

	//Main loop: iterations of kinematic closure
	for(core::Size iter=1, maxiter=(core::Size)option[maxiterations](); iter<=maxiter; iter++) {
		printf("\nMOVE %lu\n", iter); fflush(stdout);

		if(option[dab_sev_strength].user()) {
			sfxn->set_weight(core::scoring::dab_sev, 0.0); //Turn OFF dab_sev for the initial scoring
		}

		trialpose = bestpose;
		//trialpose=masterpose; //Always start fresh
		//trialpose.conformation().detect_disulfides();
		//TODO: to make this work with loops that are mixtures of L- and D-amino acids, mutate the loop SEQUENCE to the best found so far, right here.
		//set_loop_conformation(trialpose, bestpose, startres, endres, loop_length); //Set the loop conformation to the best one found so far (though the rest of the protein will be in the starting conformation).

		kinclose(trialpose, startres, endres, loop_length, is_cyclic, ((option[enumerate_only_mode]() || iter==1)?false:true), RG, sfxn, allowL, allowD, allowGLY); //Kinematic loop closure.
		(*sfxn)(trialpose);

		//A vector of positions that CAN move.  (The no_repack_positions vector overrides this -- if a position is in both lists, it does NOT move.)
		utility::vector1 < core::Size > can_repack_positions;

		//Inner loop: rounds of design, perturbation, and relaxation.
		{
			if(option[dab_sev_strength].user()) {
				sfxn->set_weight(core::scoring::dab_sev, option[dab_sev_strength]()); //Turn OFF dab_sev
			}

			core::pose::Pose bestrelaxedpose=trialpose;
			for(core::Size irounds=1, maxrounds=(core::Size)option[maxdesignrounds](); irounds<=maxrounds; irounds++) {
				trialpose2=trialpose;
				if(option[jittersize].user()) jitterloop(trialpose2, startres, loop_length, option[jittersize](),RG); //Add a small random value to all backbond dihedral angles
				if(!option[enumerate_only_mode]()) designloop(trialpose2, startres, endres, loop_length, is_cyclic, sfxn, no_repack_positions, can_repack_positions, allowed_L_aas, allowed_D_aas, allowGLY); //designloop will add to the can_repack_positions vector based on the cutoff radius, if specified.
				constrain_loop_ends(trialpose2, startres+loop_length, endres, sfxn); //Add constraints to the score function and pose to keep the termini together
				relaxloop(trialpose2, startres, endres, loop_length, frlx, minmove, mm, no_repack_positions, can_repack_positions); //relax the loop
				if(option[dab_sev_strength].user()) {
					sfxn->set_weight(core::scoring::dab_sev, 0.0); //Turn OFF dab_sev
					printf("Minimizing without dab_sev.\n"); fflush(stdout);
					minmove->apply(trialpose2); //Minimize without dab_sev
				}
				trialpose2.remove_constraints(); //Remove the constraints from the pose.
				reset_sfxn(sfxn); //Reset the scorefunction (no constraints).

				(*sfxn)(trialpose2);
				if(irounds==1 || metropolis_criterion(trialpose2.energies().total_energy(), trialpose.energies().total_energy(), option[MCtemp](), RG) ) {
					printf("\tAccepting design/relax round %lu.\n", irounds); fflush(stdout);
					trialpose=trialpose2;
					if(trialpose2.energies().total_energy() < bestrelaxedpose.energies().total_energy()) bestrelaxedpose = trialpose2; //Store this if it is the lowest-energy structure encountered so far.
				} else {
					printf("\tRejecting design/relax round %lu.\n", irounds); fflush(stdout);
				}
			}
			trialpose = bestrelaxedpose; //Make the trial pose the lowest-energy encountered during relaxation.
		}

		(*sfxn)(trialpose);
		if(option[enumerate_only_mode]() || (iter==1 && option[acceptfirstmove]()) || metropolis_criterion(trialpose.energies().total_energy(), bestpose.energies().total_energy(), option[MCtemp](), RG) ) {
			printf("Accepting move %lu.\n", iter); fflush(stdout);
			bestpose=trialpose;

			//Dump out a PDB for each move accepted
			char outfile [256];
			sprintf(outfile, "accepted_%04lu.pdb", iter);
			trialpose.dump_scored_pdb(outfile, *sfxn);
			printf("Wrote %s\n", outfile); fflush(stdout);

			if(!option[enumerate_only_mode]() && trialpose.energies().total_energy()<verybestpose.energies().total_energy()) { //If this is the lowest-energy structure found so far
				printf("New lowest-energy sequence found at move %lu.\n", iter); fflush(stdout);
				verybestpose=trialpose;
				sprintf(outfile, "lowestE_%04lu.pdb", iter);
				trialpose.dump_scored_pdb(outfile, *sfxn);
				printf("Wrote %s\n", outfile); fflush(stdout);
			}
		} else { //If we don't accept the move
			printf("Rejecting move %lu.\n", iter); fflush(stdout);
		}
		printf("END OF MOVE %lu\n", iter); fflush(stdout);
	}

	printf("\n***JOB COMPLETED***\n"); fflush(stdout);
	return 0;
}

