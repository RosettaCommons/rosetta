/*
	test_pose_copy.cc
	A pilot app to test whether poses are copied properly with automatic metal setup.
	Created by Vikram K. Mulligan, Baker Laboratory.

	File history:
		--Created 17 April 2014.
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
//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>

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

#define PI 3.1415926535897932384626433832795
#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

//OPT_KEY( Integer, peptidesize)
//OPT_KEY( Real, energycutoff)
//OPT_KEY( Integer, relaxrnds)
//OPT_KEY( IntegerVector, dpositions)

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//utility::vector1 <core::Size> emptyintlist;

	//NEW_OPT( peptidesize, "The number of residues in the peptide, including terminal cysteines.  Default 9.", 9);
	//NEW_OPT( energycutoff, "The energy above which structures are discarded.  Only used if a value is specified by the user.", 1000.0);
	//NEW_OPT( relaxrnds, "The number of rounds of relaxation with the fastrelax mover.  Default 3.", 3);
	//NEW_OPT( dpositions, "The positions at which a D-alanine should be placed.  Default empty list.", emptyintlist);
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


	printf("Starting test_pose_copy.cc\n");
	printf("Pilot app created 17 April 2014 by Vikram K. Mulligan, Baker Laboratory.\n");
	register_options(); //Parse the input flags and set global variables appropriately.
	devel::init(argc, argv); //Initialize Rosetta

	//Default scorefunction:
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//sfxn->set_weight(cart_bonded, 1.0); //Turn on cart_bonded.

	core::pose::Pose pose1;
	core::pose::Pose pose2;

	core::import_pose::pose_from_file (pose1, option[in::file::s]()[1]);
	printf("Imported pose1.\n"); fflush(stdout);

	pose2=pose1;

	pose1.dump_scored_pdb("pose1.pdb", *sfxn);
	printf("Wrote pose1.pdb.\n"); fflush(stdout);
	pose2.dump_scored_pdb("pose2.pdb", *sfxn);
	printf("Wrote pose2.pdb.\n"); fflush(stdout);

	printf("***JOB COMPLETED.***\n");
	return 0;
}

