/****************************************************************************************************
	test_KIC_betapeptide.cc

	A pilot app to test kinematic closure with beta-3-amino acids (and possibly other nonstandard
	backbones?)

	History:
	--File created 22 August 2013 by Vikram K. Mulligan, Baker Laboratory.
****************************************************************************************************/

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/format.hh>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <protocols/relax/FastRelax.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/Energies.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>


using ObjexxFCL::fmt::F;
using namespace protocols::cluster;


//OPT_1GRP_KEY( FileVector, ccd, rmsd_matrix )
//OPT_KEY( Real, v_bb_perturbation ) //The average number of degrees by which backbone dihedrals are perturbed; default 3
//OPT_KEY ( Boolean, v_use_sidechains ) //Should I use side chain dihedrals or not?  Default true (i.e. use side chains).
//OPT_KEY ( String, v_sequence ) //The sequence of beta-amino acids.
OPT_KEY( Real, v_MCtemperature) //The value of k_b*T to use for the Monte Carlo search of backbone conformation space (default 1.5).
OPT_KEY( Integer, v_MCrounds) //The number of rounds of Monte Carlo searches of backbone conformation space (default 10).
OPT_KEY( Integer, v_MCsteps_per_round) //The number of steps in each of the rounds of Monte Carlo searches of backbone conformation space (default 10).
OPT_KEY( Real, v_MCanglepert) //The amount by which backbone dihedral angles (other than omega) are perturbed during the Monte Carlo search of backbone conformation space (default 10 degrees).

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//NEW_OPT (v_sequence, "", "AAAAAA");
	//NEW_OPT( v_use_sidechains, "", true);
	//NEW_OPT ( v_bb_perturbation, "", 3);
	NEW_OPT( v_MCtemperature, "", 1.5);
	NEW_OPT( v_MCrounds, "", 10);
	NEW_OPT( v_MCsteps_per_round, "", 10);
	NEW_OPT( v_MCanglepert, "", 2.0);
}

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

using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

//Function to return a random dihedral angle.
core::Real randangle (numeric::random::RandomGenerator &RG)
{	return (core::Real)(RG.uniform()*360.0-180.0);
}

//Function to randomly jitter backbone phi, psi, and theta.
void perturb_bb (
	core::pose::Pose& pose,
	numeric::random::RandomGenerator &RG
) {
	for (core::Size ir=1; ir<=pose.n_residue(); ir++) {//loop through all residues
		if(pose.residue(ir).has("CM")) { //If this is a beta-amino acid, perturb phi, psi, and theta.
			betapeptide_setphi(pose, ir, randangle(RG));
			betapeptide_settheta(pose, ir, randangle(RG));
			betapeptide_setpsi(pose, ir, randangle(RG));
		} else { //If this is an alpha-amino acid, perturb phi and psi.
			pose.set_phi(ir, randangle(RG));
			pose.set_psi(ir, randangle(RG));
		}//end if
	} //end loop through all residues.
	return;
}

//Funciton to do a local MC search of backbone conformation space to find the nearest low-energy conformation.
void jitter_and_minimize (
	core::pose::Pose &pose,
	numeric::random::RandomGenerator &RG,
	protocols::relax::FastRelax &frlx,
	core::scoring::ScoreFunctionOP &sfxn,
	core::Size MCrounds,
	core::Size MCsteps_per_round,
	core::Real anglepertsize,
	core::Real MCtemperature
) {
	using namespace core::id;

	core::pose::Pose trialpose, lastacceptpose, lowestEpose, verylowestEpose;
	verylowestEpose = pose;

	//Begin rounds of MC searches
	for (core::Size iround=1; iround<=MCrounds; iround++) {
		lowestEpose=pose;
		lastacceptpose=pose;

		//Begin a Monte Carlo search
		for (core::Size istep=1; istep<=MCsteps_per_round; istep++) {
			trialpose = lastacceptpose;
			//Loop through all residues and wiggle all backbone dihedral angles except omega:
			for(core::Size ir=1; ir<=trialpose.n_residue(); ir++) {
				if(trialpose.residue(ir).has("CM")) { //If this is a beta-amino acid, perturb phi, psi, and theta.
					betapeptide_setphi(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 1))+(core::Real)RG.gaussian()*anglepertsize );
					betapeptide_settheta(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 2))+(core::Real)RG.gaussian()*anglepertsize);
					betapeptide_setpsi(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 3))+(core::Real)RG.gaussian()*anglepertsize);
				} else { //If this is an alpha-amino acid, perturb phi and psi.
					trialpose.set_phi(ir, trialpose.phi(ir)+(core::Real)RG.gaussian()*anglepertsize);
					trialpose.set_psi(ir, trialpose.psi(ir)+(core::Real)RG.gaussian()*anglepertsize);
				}//end if
			}
			frlx.apply(trialpose); //Relax the jittered structure.
			(*sfxn)(trialpose); //Score the jittered structure.

			//Apply the Metropolis criterion to accept or reject the jittering:
			if( trialpose.energies().total_energy() < lastacceptpose.energies().total_energy() ) { //If the energy is lower, accept the move.
				lastacceptpose=trialpose;
				if (lastacceptpose.energies().total_energy()<lowestEpose.energies().total_energy()) lowestEpose = lastacceptpose; //If this is the lowest energy encountered this round.
			} else { //If the energy has gone up, accept or reject based on Metropolis criterion:
				if( (core::Real)RG.uniform() < (core::Real)exp( -1.0*( trialpose.energies().total_energy() - lastacceptpose.energies().total_energy() )/MCtemperature) ) {
					lastacceptpose=trialpose;
				} //no else required, here.
			}//end applying MC criterion	

		}//End this MC trajectory (one in a series of rounds)

		if (lowestEpose.energies().total_energy() < verylowestEpose.energies().total_energy() ) verylowestEpose=lowestEpose; //If this round has produced the lowest energy so far.

	} //End the rounds of MC searches.

	//printf("Initial energy: %f\nEnergy after MC search: %f\n", pose.energies().total_energy(), verylowestEpose.energies().total_energy()); fflush(stdout);

	pose=verylowestEpose; //Spit out the very lowest energy pose generated.

	return;
}

int main( int argc, char * argv [] ) {

	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace protocols::cluster;
	using namespace std;
	using namespace chemical;
	using namespace conformation;
	using namespace core::pose;

	printf("Starting test_KIC_betapeptide.cc\nFile created 22 August 2013 by Vikram K. Mulligan\n\n");
	fflush(stdout);

	numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

	register_options();
	devel::init(argc, argv);

	core::Size numstructs = 25;
	if(option[out::nstruct].user()) numstructs=(core::Size) option[out::nstruct]();

	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::getScoreFunction();

	protocols::relax::FastRelax frlx(sfxn, 1); //A fastrelax mover
	protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn); //Create the RepackSidechains mover and set the score function.

	core::pose::Pose mypose;
	const string sequence = "GA[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]G";
	//const string sequence = "GA[B3A]A[B3A]A[B3A]A[DALA]A[DALA]AG";

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( sequence, *rsd_set, false );
	mypose.clear();

	for ( core::Size i = 1, ie = requested_types.size(); i <= ie; ++i ) {
		chemical::ResidueType const & rsd_type = *requested_types[ i ];
		core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );
		if(i==1) {
			mypose.append_residue_by_jump( *new_rsd, 1, "", "", true );
			core::pose::add_lower_terminus_type_to_pose_residue(mypose, 1);
		} else {
			mypose.append_residue_by_bond( *new_rsd, true );
		}
	}

	core::pose::add_upper_terminus_type_to_pose_residue(mypose, mypose.n_residue());

	for(core::Size ir=1; ir<=mypose.n_residue(); ir++) {
		if(mypose.residue(ir).has("CM")) {
			betapeptide_setphi(mypose, ir, -122.6);
			betapeptide_settheta(mypose, ir, -58.4);
			betapeptide_setpsi(mypose, ir, -165.3);
			betapeptide_setomega(mypose, ir, 180.0);
		} else {
			mypose.set_phi(ir, -60.0);
			if(ir<mypose.n_residue()) {
				mypose.set_omega(ir, 180.0);
				mypose.set_psi(ir, 160.0);
			}
		}
	}

	betapeptide_setphi(mypose, 9, 45.0);
	betapeptide_setphi(mypose, 4, -26.3);
	betapeptide_settheta(mypose, 6, 58.4);
	mypose.update_residue_neighbors();

	mypose.dump_scored_pdb("first.pdb", *sfxn); //Pre-KIC pose

	//Kinematic closure mover:
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinmover = new protocols::loops::loop_closure::kinematic_closure::KinematicMover;
	kinmover->set_temperature( 1.0 );
	kinmover->set_vary_bondangles( false );
	kinmover->set_sample_nonpivot_torsions( true ); //true
	kinmover->set_rama_check( true );
	kinmover->set_idealize_loop_first(true); //true
	kinmover->set_sfxn(sfxn);
	kinmover->set_pivots(2, 5, mypose.n_residue()-1);

	for(core::Size i=1; i<=numstructs; i++) {
		core::pose::Pose temppose=mypose;

		kinmover->apply(temppose);

		betapeptide_setomega(temppose, 7, 100.0);
		temppose.update_residue_neighbors();

		char outfilename [256];
		sprintf(outfilename, "out_%04lu.pdb", i);
		temppose.dump_scored_pdb(outfilename, *sfxn);

		frlx.apply(temppose);
		sprintf(outfilename, "relaxed_%04lu.pdb", i);
		temppose.dump_scored_pdb(outfilename, *sfxn);

	}

	printf("\nJOB COMPLETED.\n"); fflush(stdout);
	return 0;
}

