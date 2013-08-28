/****************************************************************************************************
	make_ideal_betapeptide.cc

	A pilot app to build a beta-peptide, and to generate ideal secondary structure.

	History:
	--File created 4 January 2013 by Vikram K. Mulligan, Baker Laboratory.
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
	int resnumber,
	double angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 1 ), angle );
	return;
}

void betapeptide_settheta (
	core::pose::Pose &pose,
	int resnumber,
	double angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 2 ), angle );
	return;
}

void betapeptide_setpsi (
	core::pose::Pose &pose,
	int resnumber,
	double angle
) {
	using namespace core;
	using namespace core::id;
	pose.set_torsion( TorsionID( resnumber, id::BB, 3 ), angle );
	return;
}

void betapeptide_setomega (
	core::pose::Pose &pose,
	int resnumber,
	double angle
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
double randangle (numeric::random::RandomGenerator &RG)
{	return (double)(RG.uniform()*360.0-180.0);
}

//Function to randomly jitter backbone phi, psi, and theta.
void perturb_bb (
	core::pose::Pose& pose,
	numeric::random::RandomGenerator &RG
) {
	for (int ir=1; ir<=pose.n_residue(); ir++) {//loop through all residues
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
	int MCrounds,
	int MCsteps_per_round,
	double anglepertsize,
	double MCtemperature
) {
	using namespace core::id;

	core::pose::Pose trialpose, lastacceptpose, lowestEpose, verylowestEpose;
	verylowestEpose = pose;

	//Begin rounds of MC searches
	for (int iround=1; iround<=MCrounds; iround++) {
		lowestEpose=pose;
		lastacceptpose=pose;

		//Begin a Monte Carlo search
		for (int istep=1; istep<=MCsteps_per_round; istep++) {
			trialpose = lastacceptpose;
			//Loop through all residues and wiggle all backbone dihedral angles except omega:
			for(int ir=1; ir<=trialpose.n_residue(); ir++) {
				if(trialpose.residue(ir).has("CM")) { //If this is a beta-amino acid, perturb phi, psi, and theta.
					betapeptide_setphi(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 1))+(double)RG.gaussian()*anglepertsize );
					betapeptide_settheta(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 2))+(double)RG.gaussian()*anglepertsize);
					betapeptide_setpsi(trialpose, ir, trialpose.torsion(TorsionID(ir, id::BB, 3))+(double)RG.gaussian()*anglepertsize);
				} else { //If this is an alpha-amino acid, perturb phi and psi.
					trialpose.set_phi(ir, trialpose.phi(ir)+(double)RG.gaussian()*anglepertsize);
					trialpose.set_psi(ir, trialpose.psi(ir)+(double)RG.gaussian()*anglepertsize);
				}//end if
			}
			frlx.apply(trialpose); //Relax the jittered structure.
			(*sfxn)(trialpose); //Score the jittered structure.

			//Apply the Metropolis criterion to accept or reject the jittering:
			if( trialpose.energies().total_energy() < lastacceptpose.energies().total_energy() ) { //If the energy is lower, accept the move.
				lastacceptpose=trialpose;
				if (lastacceptpose.energies().total_energy()<lowestEpose.energies().total_energy()) lowestEpose = lastacceptpose; //If this is the lowest energy encountered this round.
			} else { //If the energy has gone up, accept or reject based on Metropolis criterion:
				if( (double)RG.uniform() < (double)exp( -1.0*( trialpose.energies().total_energy() - lastacceptpose.energies().total_energy() )/MCtemperature) ) {
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

	printf("Starting make_ideal_betapeptide.cc\nFile created 4 January 2013 by Vikram K. Mulligan\n\n");
	fflush(stdout);

	numeric::random::RandomGenerator RG( 923749 ); //Random generator and seed

	register_options();
	devel::init(argc, argv);

	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::getScoreFunction();

	protocols::relax::FastRelax frlx(sfxn, 1); //A fastrelax mover
	protocols::simple_moves::RepackSidechainsMover repack_sc(sfxn); //Create the RepackSidechains mover and set the score function.
	
	core::pose::Pose mypose;
	const string sequence = "GA[B3A]A[B3A]A[B3A]A[B3A]A[B3A]A[B3A]G";

	// grab residue types
	chemical::ResidueTypeCAPs requested_types = core::pose::residue_types_from_sequence( sequence, *( chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" )), true );
	mypose.clear();

	int n_struct = option[out::nstruct](); //Number of structures to generate
	string outprefix = "";
	if (option[out::file::o].user()) outprefix= option[out::file::o]()+"_"; //Output file name prefix
	string outfile = "";
	char charbuffer[24];

	for ( Size i = 1; i <= requested_types.size(); ++i )
	{
		// grab the new residue
		chemical::ResidueType const & rsd_type = *requested_types[ i ];
		core::conformation::ResidueOP new_rsd( NULL );
		new_rsd = conformation::ResidueFactory::create_residue ( rsd_type );
		if (i>1) mypose.append_residue_by_bond( *new_rsd, true );
		else {
			mypose.append_residue_by_jump(*new_rsd, 1, "", "", true);
			core::pose::add_lower_terminus_type_to_pose_residue(mypose, 1);
		}
	}

	core::pose::add_upper_terminus_type_to_pose_residue(mypose, mypose.n_residue());

	double philist[] = {	-122.6,	-71.2,	58.9,	-53.6,	50.4,	64.8,	-175.1,	-110.4,	-67.2,	59.5,	83.5,	-171.8,	-161.9,	63.7,	59.6	};
	double thetalist[] = {	-58.4,	134.0,	-121.4,	-47.7,	51.7,	63.3,	52.9,	60.5,	170.4,	166.2,	-53.0,	-70.6,	76.3,	76.7,	57.7	};
	double psilist[] = {	-165.3,	-64.8,	78.3,	106.5,	-109.6,	-172.6,	107.0,	21.8,	-177.6,	168.4,	-92.7,	31.9,	-36.4,	-54.4,	84.7	};

	printf("\nFile\tphi\ttheta\tpsi\tEnergy\n----\t---\t-----\t---\t------\n"); fflush(stdout);
	for (int iconf = 0; iconf<15; iconf++)
	{
		core::pose::Pose temppose=mypose;
		temppose.set_omega(1,180);
		for(int ir=2; ir<temppose.n_residue(); ir++) {
			betapeptide_setphi(temppose, ir, philist[iconf]);
			betapeptide_settheta(temppose, ir, thetalist[iconf]);
			betapeptide_setpsi(temppose, ir, psilist[iconf]);
			if (ir<temppose.n_residue()) betapeptide_setomega(temppose, ir, 180);
		}
		temppose.update_residue_neighbors();
		repack_sc.apply(temppose); //Repack side-chains.
		temppose.update_residue_neighbors();
		frlx.apply(temppose);
		jitter_and_minimize(temppose, RG, frlx, sfxn, option[v_MCrounds](), option[v_MCsteps_per_round](), option[v_MCanglepert](), option[v_MCtemperature]());
		(*sfxn)(temppose);
		char curnum[8];
		sprintf(curnum, "%02i", iconf+1);
		temppose.dump_scored_pdb( outprefix + "B" + curnum + ".pdb", *sfxn );
		printf("%s:\t%.1f\t%.1f\t%.1f\t%.6f\n", (outprefix + "B" + curnum + ".pdb").c_str(), philist[iconf], thetalist[iconf], psilist[iconf], temppose.energies().total_energy());
		fflush(stdout);
	}

	printf("\nJOB COMPLETED.\n"); fflush(stdout);
	return 0;
}

