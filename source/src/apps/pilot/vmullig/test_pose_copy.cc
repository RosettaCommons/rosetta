/*
	test_pose_copy.cc
	A pilot app to test whether poses are copied properly with automatic metal setup.
	Created by Vikram K. Mulligan, Baker Laboratory.

	File history:
		--Created 17 April 2014.

*/

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/conformation/Residue.fwd.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <sstream>
#include <stdio.h>
//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/minimization_packing/MinMover.hh>


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

