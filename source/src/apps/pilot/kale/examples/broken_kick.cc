// Headers {{{1
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/balanced_kic/algorithms.hh>
#include <devel/init.hh>
#include <iomanip>

// Namespaces {{{1
using namespace std;
using namespace core;
using namespace devel;
using namespace devel::balanced_kic::algorithms;
// }}}1

void make_correct_perturbation( // {{{1
		ClosureProblem &problem) {

	problem.torsion_angles[ 1] = 29.2461;
	problem.torsion_angles[ 2] = 171.427;
	problem.torsion_angles[ 3] = 176.577;

	problem.torsion_angles[ 4] = 77.1789;
	problem.torsion_angles[ 5] = 158.531;
	problem.torsion_angles[ 6] = 178.858;

	problem.torsion_angles[ 7] = -19.0967;
	problem.torsion_angles[ 8] = -56.2946;
	problem.torsion_angles[ 9] = 178.923;

	problem.torsion_angles[10] = -104.399;
	problem.torsion_angles[11] = 58.4015;
	problem.torsion_angles[12] = 178.845;

	problem.torsion_angles[13] = 56.9483;
	problem.torsion_angles[14] = 45.5138;
	problem.torsion_angles[15] = 181.292;

	problem.torsion_angles[16] = 58.7624;
	problem.torsion_angles[17] = -140.242;
	problem.torsion_angles[18] = 179.689;

	problem.torsion_angles[19] = -72.8008;
	problem.torsion_angles[20] = 35.8359;
	problem.torsion_angles[21] = 177.923;

	problem.torsion_angles[22] = 131.783;
	problem.torsion_angles[23] = -111.547;
	problem.torsion_angles[24] = 177.051;

	problem.torsion_angles[25] = -172.692;
	problem.torsion_angles[26] = -93.694;
	problem.torsion_angles[27] = 176.571;

	problem.torsion_angles[28] = -64.8262;
	problem.torsion_angles[29] = 145.006;
	problem.torsion_angles[30] = 179.204;

	problem.torsion_angles[31] = 86.5967;
	problem.torsion_angles[32] = -61.3628;
	problem.torsion_angles[33] = 176.83;
                           
	problem.torsion_angles[34] = -74.2586;
	problem.torsion_angles[35] = -167.88;
	problem.torsion_angles[36] = 227.887;
                           
	problem.torsion_angles[37] = 318.961;
	problem.torsion_angles[38] = 64.3249;
	problem.torsion_angles[39] = 323.456;
}

void pick_correct_solution( // {{{1
		SolutionList const &solutions,
		ClosureSolution const *&solution) {

	solution = &solutions[4];
}
// }}}1

void make_pathological_perturbation(ClosureProblem &problem) { // {{{1

	/*
	 * torsion_angles     29.2461    171.427    176.577     147.79    12.8668    178.858   -80.5171    -143.65    178.923     14.635   -28.4338    178.845    39.2256   -6.25637    181.292    48.1442   -83.4468    179.689   -169.399   -42.8457    177.923    99.1484   -83.8437    177.051    157.905    93.0474    176.571    105.486    150.429    179.204   -127.218   -59.8358     176.83   -13.9848    74.3853    227.887    318.961    64.3249    323.456 
	 * bond_angles        163.122    110.949    111.247    114.693    110.981    114.701    120.399    110.977     113.26    120.354    110.976    115.985     123.07    110.942    113.978    120.437    110.995    113.837    120.808    110.998     116.09    122.525    110.966    114.171     120.44     110.96    114.442    120.881    110.995    115.729    121.958    111.023    114.901    121.872    111.011    121.441    142.219    110.958    126.139 
	 * bond_lengths       1.46079    1.51089    1.38058    1.45975    1.50969    1.35064    1.45991    1.50969    1.32866    1.45993    1.51056    1.33521    1.45958    1.51033    1.32498    1.45978    1.50988    1.33911    1.45966    1.50935     1.3387    1.45966     1.5105    1.34395    1.46054    1.50979     1.3149    1.45971    1.50918    1.33328    1.45981    1.50992    1.30402    1.45931    1.51015    1.26892     1.4601    1.51075    3.07789 
	 * origin               1.628      0.206      0.303 
	 * frame             0.403206 -0.0773553  -0.911834   0.631593  -0.697521    0.33846  -0.662205  -0.712377  -0.232387 
	 * solution        6 of 6
	 */

	problem.torsion_angles[ 1] = 29.2461;
	problem.torsion_angles[ 2] = 171.427;
	problem.torsion_angles[ 3] = 176.577;
                               
	problem.torsion_angles[ 4] = 147.79;
	problem.torsion_angles[ 5] = 12.8668;
	problem.torsion_angles[ 6] = 178.858;
                               
	problem.torsion_angles[ 7] = -80.5171;
	problem.torsion_angles[ 8] = -143.65;
	problem.torsion_angles[ 9] = 178.923;
                               
	problem.torsion_angles[10] = 14.635;
	problem.torsion_angles[11] = -28.4338;
	problem.torsion_angles[12] = 178.845;
                               
	problem.torsion_angles[13] = 39.2256;
	problem.torsion_angles[14] = -6.25637;
	problem.torsion_angles[15] = 181.292;
                               
	problem.torsion_angles[16] = 48.1442;
	problem.torsion_angles[17] = -83.4468;
	problem.torsion_angles[18] = 179.689;
                               
	problem.torsion_angles[19] = -169.399;
	problem.torsion_angles[20] = -42.8457;
	problem.torsion_angles[21] = 177.923;
                               
	problem.torsion_angles[22] = 99.1484;
	problem.torsion_angles[23] = -83.8437;
	problem.torsion_angles[24] = 177.051;
                               
	problem.torsion_angles[25] = 157.905;
	problem.torsion_angles[26] = 93.0474;
	problem.torsion_angles[27] = 176.571;
                               
	problem.torsion_angles[28] = 105.486;
	problem.torsion_angles[29] = 150.429;
	problem.torsion_angles[30] = 179.204;
                               
	problem.torsion_angles[31] = -127.218;
	problem.torsion_angles[32] = -59.8358;
	problem.torsion_angles[33] = 176.83;
                               
	problem.torsion_angles[34] = -13.9848;
	problem.torsion_angles[35] = 74.3853;
	problem.torsion_angles[36] = 227.887;
                               
	problem.torsion_angles[37] = 318.961;
	problem.torsion_angles[38] = 64.3249;
	problem.torsion_angles[39] = 323.456;

	/*
	problem.torsion_angles[ 1] = 29.2461;
	problem.torsion_angles[ 2] = 171.427;
	problem.torsion_angles[ 3] = 176.577;

	problem.torsion_angles[ 4] = 77.1789;
	problem.torsion_angles[ 5] = 158.531;
	problem.torsion_angles[ 6] = 178.858;

	problem.torsion_angles[ 7] = -19.0967;
	problem.torsion_angles[ 8] = -56.2946;
	problem.torsion_angles[ 9] = 178.923;

	problem.torsion_angles[10] = -104.399;
	problem.torsion_angles[11] = 58.4015;
	problem.torsion_angles[12] = 178.845;

	problem.torsion_angles[13] = 56.9483;
	problem.torsion_angles[14] = 45.5138;
	problem.torsion_angles[15] = 181.292;

	problem.torsion_angles[16] = 58.7624;
	problem.torsion_angles[17] = -140.242;
	problem.torsion_angles[18] = 179.689;

	problem.torsion_angles[19] = -72.8008;
	problem.torsion_angles[20] = 35.8359;
	problem.torsion_angles[21] = 177.923;

	problem.torsion_angles[22] = 131.783;
	problem.torsion_angles[23] = -111.547;
	problem.torsion_angles[24] = 177.051;

	problem.torsion_angles[25] = -172.692;
	problem.torsion_angles[26] = -93.694;
	problem.torsion_angles[27] = 176.571;

	problem.torsion_angles[28] = -64.8262;
	problem.torsion_angles[29] = 145.006;
	problem.torsion_angles[30] = 179.204;

	problem.torsion_angles[31] = 86.5967;
	problem.torsion_angles[32] = -61.3628;
	problem.torsion_angles[33] = 176.83;
                           
	problem.torsion_angles[34] = -74.2586;
	problem.torsion_angles[35] = -167.88;
	problem.torsion_angles[36] = 227.887;
                           
	problem.torsion_angles[37] = 318.961;
	problem.torsion_angles[38] = 64.3249;
	problem.torsion_angles[39] = 323.456;
	*/
}

void pick_pathological_solution( // {{{1
		SolutionList const &solutions,
		ClosureSolution const *&solution) {

	solution = &solutions[6];
}
// }}}1

// This program demonstrates that torsions generated by the KIC solver are 
// faithfully translated into the final PDB structure.  Consequently, the bug 
// is most likely in the KIC solver itself.  To prove that this is the case, I 
// will try to elicit the bug using the original protocols::KinematicMover.

int main(int argc, char** argv) {
	devel::init(argc, argv);

	pose::Pose pose;
	string input_pdb = "structures/cyclic/14.marked.pdb";
	import_pose::pose_from_file(pose, input_pdb, core::import_pose::PDB_file);

	ClosureProblem problem;
	ClosureSolution const *solution = NULL;
	SolutionList solutions;

	pose.dump_pdb("problem.pdb");

	define_closure_problem(problem, pose, 2, 7, 12);
	make_pathological_perturbation(problem);
	solve_closure_problem(problem, solutions);

	for (int i = 1; i <= solutions.size(); i++) {
		cout << "Torsions (" << i << "/6)" << endl; 
		cout << "==============" << endl;
		for (int j = 1; j <= solutions[i].torsion_angles.size(); j++) {
			Real angle = solutions[i].torsion_angles[j];
			angle = (angle < 0) ? angle + 360 : angle;
			cout << fixed << setprecision(4) << setw(8) << angle << endl;
		}
		cout << endl;

		cout << "Angles (" << i << "/6)" << endl; 
		cout << "============" << endl;
		for (int j = 1; j <= solutions[i].bond_angles.size(); j++) {
			Real angle = solutions[i].bond_angles[j];
			angle = (angle < 0) ? angle + 360 : angle;
			cout << fixed << setprecision(4) << setw(8) << angle << endl;
		}
		cout << endl;

		cout << "Lengths (" << i << "/6)" << endl; 
		cout << "=============" << endl;
		for (int j = 1; j <= solutions[i].bond_lengths.size(); j++) {
			Real length = solutions[i].bond_lengths[j];
			cout << fixed << setprecision(4) << setw(8) << length << endl;
		}
		cout << endl;
	}

	apply_all_closure_solutions(solutions);
	pick_pathological_solution(solutions, solution);

	solutions[1].pose.dump_pdb("solution.1.pdb");
	solutions[2].pose.dump_pdb("solution.2.pdb");
	solutions[3].pose.dump_pdb("solution.3.pdb");
	solutions[4].pose.dump_pdb("solution.4.pdb");
	solutions[5].pose.dump_pdb("solution.5.pdb");
	solutions[6].pose.dump_pdb("solution.6.pdb");
}
