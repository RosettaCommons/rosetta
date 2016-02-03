/*
	test_sidechainKIC.cc
	A pilot app to test kinematic closure through sidechains.

	File history:
		--Created 20 March 2014.
*/

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <devel/init.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <stdio.h>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/relax/FastRelax.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>

#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <core/id/NamedAtomID.hh>

#define PI 3.1415926535897932384626433832795

OPT_KEY( Boolean, buildideal)

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//utility::vector1 <core::Size> emptyintlist;

	NEW_OPT( buildideal, "Should we build ideal geometry?  Default true.", true);
}


/**********
  MAIN!!!
**********/
int main(int argc, char *argv[]) {
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options; 
	using namespace basic::options::OptionKeys; 
	using namespace core::scoring::constraints; 
	using namespace core::id;
	using namespace core::pose;
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace protocols::generalized_kinematic_closure;

	//using namespace protocols::toolbox::match_enzdes_util; //Florian's stuff


	printf("Starting test_sidechainKIC.cc\n");
	printf("Pilot app created 20 March 2014 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	if(!option[in::file::s].user()) {
		printf("Error!  The user must specify an input file with the -in:file:s flag.  Crashing.\n"); fflush(stdout);
		exit(1);
	}

	//Default scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.


	core::pose::Pose mypose;
	core::import_pose::pose_from_file (mypose, option[in::file::s]()[1]);
	//const core::Size nres = mypose.n_residue();
	printf("Import complete.\n"); fflush(stdout);

	mypose.dump_scored_pdb("1_initial.pdb", (*sfxn)); printf("Wrote 1_initial.pdb.\n"); fflush(stdout);

	for(core::Size counter=1; counter<=100; counter++) {
		//core::Real const dihed = 231.0+25.0*cos(3.141592654/50.0*(core::Real)counter);
		//core::Real const bondang = 124.6 + 10.0*cos(3.141592654/50.0*(core::Real)counter);
		core::Real const bondlen = 2.03055 + 1.0*sin(3.141592654/50.0*(core::Real)counter);

		core::pose::Pose mypose2=mypose;

		GeneralizedKICOP scKIC = new GeneralizedKIC;
		scKIC->set_mover_effect_on_bonded_geometry(1);
		scKIC->set_build_ideal_geometry( option[buildideal]() );

		for(core::Size i=18; i<=20; i++) {
			scKIC->add_loop_residue(i);
		}
		for(core::Size i=40; i>=38; i--) {
			scKIC->add_loop_residue(i);
		}

		scKIC->set_pivot_atoms( 18, "CA", 40, "CA" , 38, "CA" );
		scKIC->set_closure_attempts(1);

		//Create a set_dihedrals perturber:
		//scKIC->add_perturber("set_dihedral");
		//scKIC->add_perturber("set_bondangle");
		scKIC->add_perturber("set_bondlength");
		/*utility::vector1 < core::id::AtomID > dihedralatoms;
		dihedralatoms.push_back( AtomID(3, 18) );
		dihedralatoms.push_back( AtomID(1, 19) );
		dihedralatoms.push_back( AtomID(2, 19) );
		dihedralatoms.push_back( AtomID(3, 19) );*/
		//utility::vector1 < core::id::AtomID > bondangleatoms;
		//bondangleatoms.push_back( AtomID(2, 20) ); //CYS20 CA
		//bondangleatoms.push_back( AtomID(5, 20) ); //CYS20 CB
		//bondangleatoms.push_back( AtomID(6, 20) ); //CYS20 SG
		utility::vector1 < core::id::NamedAtomID > bondlenatoms;
		bondlenatoms.push_back( NamedAtomID("SG", 20) ); //CYS20 SG
		bondlenatoms.push_back( NamedAtomID("SG", 40) ); //CYS40 SG
		//scKIC->add_AtomIDset_to_perturber_AtomIDset_list(1, dihedralatoms);
		//scKIC->add_value_to_perturber_value_list(1, dihed);
		//scKIC->add_AtomIDset_to_perturber_AtomIDset_list(1, bondangleatoms);
		//scKIC->add_value_to_perturber_value_list(1, bondang);
		scKIC->add_atomset_to_perturber_atomset_list(1, bondlenatoms);
		scKIC->add_value_to_perturber_value_list(1, bondlen);

		scKIC->apply(mypose2);

		char filename[256];
		sprintf(filename, "out_%04lu.pdb", counter);

		mypose2.dump_pdb(filename); printf("Wrote %s.\n", filename); fflush(stdout);

		//mypose.dump_scored_pdb("2_final.pdb", (*sfxn)); printf("Wrote 2_final.pdb.\n"); fflush(stdout);
	}

	printf("*****JOB COMPLETED.*****\n");

	return 0;
}

