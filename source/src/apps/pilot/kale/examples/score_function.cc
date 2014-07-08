// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <devel/init.hh>

// C++ headers
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {

	devel::init(argc, argv);

	using core::pose::Pose;
	using core::import_pose::pose_from_pdb;

	Pose pose;
	pose_from_pdb(pose, "structures/kic/1srp.pdb");

	using core::scoring::ScoreFunctionOP;
	using core::scoring::get_score_function;

	ScoreFunctionOP function = get_score_function();

	cout << function->score(pose) << endl;

}
