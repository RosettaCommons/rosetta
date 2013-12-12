// Headers {{{1
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>

#include <devel/balanced_kic/algorithms.hh>
#include <devel/init.hh>

#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <iomanip>
#include <unistd.h>
#include <boost/lexical_cast.hpp>


// Namespaces {{{1
using namespace std;
using namespace core;
using namespace devel;
using namespace devel::balanced_kic::algorithms;
using namespace protocols::loops::loop_closure::kinematic_closure;
using namespace basic::options;

using core::Size;
using core::Real;

// Global Variables {{{1
static numeric::random::RandomGenerator RG(431375);

OPT_1GRP_KEY(Boolean, kale, terminus)

// }}}1

class NullKinematicPerturber : public KinematicPerturber { // {{{1
	
public:
	
	std::string perturber_type() const {
		return "NullKinematicPerturber";
	}
	
	void perturb_chain(
			core::pose::Pose const & pose,
			utility::vector1< core::Real> & torsion_angles,
			utility::vector1< core::Real> & bond_angles,
			utility::vector1< core::Real> & bond_lengths) {

		/// Don't do anything.
		return;
	}
};
// }}}1

void make_random_perturbation( // {{{1
		ClosureProblem &problem) {

	Size first_torsion = 4;
	Size last_torsion = problem.torsion_angles.size() - 3;

	for(Size i = first_torsion; i <= last_torsion; /* increment in loop */) {
			Real phi = 360 * RG.uniform() - 180;
			Real psi = 360 * RG.uniform() - 180;
			Real omega = 360 * RG.uniform() - 180;

			problem.torsion_angles[i++] = phi;
			problem.torsion_angles[i++] = psi;

			if (i != last_torsion) problem.torsion_angles[i++] = omega;
			else i++;
	}
}

void pick_random_solution( // {{{1
		SolutionList const &solutions,
		ClosureSolution const *&solution) {

	Size index = RG.random_range(1, solutions.size());
	solution = &solutions[index];
}
// }}}1

void devel_main(int argc, char** argv) { // {{{1
	chdir("../src/apps/pilot/kale/examples");
	NEW_OPT(kale::terminus, "Attempt to close the N-terminus.", false);

	devel::init(argc, argv);
	RG.set_seed(0);

	pose::Pose pose;
	string input_pdb, output_pdb;
	Size first_pivot, second_pivot, third_pivot;

	ClosureProblem problem;
	ClosureSolution const *solution = NULL;
	SolutionList solutions;

	if (option[OptionKeys::kale::terminus]() == true) {
		input_pdb = "structures/linear/4.helix.pdb";
		output_pdb = "terminus_closure.";
		first_pivot = 1; second_pivot = 2; third_pivot = 3;
	}
	else {
		input_pdb = "structures/linear/5.helix.pdb";
		output_pdb = "standard_closure.";
		first_pivot = 2; second_pivot = 3; third_pivot = 4;
	}

	cout << "Input: " << input_pdb << endl;
	cout << "Outputs: " << output_pdb << "*.pdb" << endl;
	cout << "Pivots: " << first_pivot << "/"
		                 << second_pivot << "/"
		                 << third_pivot << endl << endl;

	import_pose::pose_from_pdb(pose, input_pdb);

	cout << "Before Anything" << endl;
	cout << "===============" << endl;

	cout << "Pose XYZs" << endl; 
	cout << "---------" << endl;
	for (int j = 1; j <= pose.total_residue(); j++) {
		id::NamedAtomID n_id("N", j);
		id::NamedAtomID ca_id("CA", j);
		id::NamedAtomID c_id("C", j);

		core::PointPosition n_xyz  = pose.xyz(n_id);
		core::PointPosition ca_xyz = pose.xyz(ca_id);
		core::PointPosition c_xyz  = pose.xyz(c_id);

		cout << fixed << setprecision(4) << setw(8) << n_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << n_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << n_xyz(3) << endl;

		cout << fixed << setprecision(4) << setw(8) << ca_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << ca_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << ca_xyz(3) << endl;

		cout << fixed << setprecision(4) << setw(8) << c_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << c_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << c_xyz(3) << endl;

		cout << endl;
	}
	cout << endl;

	define_closure_problem(problem, pose, 
			first_pivot, second_pivot, third_pivot);

	cout << "Problem" << endl;
	cout << "=======" << endl;

	cout << "KIC XYZs" << endl; 
	cout << "--------" << endl;
	for (int j = 1; j <= problem.atom_xyzs.size(); j++) {
		Coordinate xyz = problem.atom_xyzs[j];
		cout << fixed << setprecision(4) << setw(8) << xyz[1] << " ";
		cout << fixed << setprecision(4) << setw(8) << xyz[2] << " ";
		cout << fixed << setprecision(4) << setw(8) << xyz[3] << endl;
		if (j % 3 == 0) cout << endl;
	}
	cout << endl;

	cout << "Pose XYZs" << endl; 
	cout << "---------" << endl;
	for (int j = 1; j <= problem.pose.total_residue(); j++) {
		id::NamedAtomID n_id("N", j);
		id::NamedAtomID ca_id("CA", j);
		id::NamedAtomID c_id("C", j);

		core::PointPosition n_xyz  = problem.pose.xyz(n_id);
		core::PointPosition ca_xyz = problem.pose.xyz(ca_id);
		core::PointPosition c_xyz  = problem.pose.xyz(c_id);

		cout << fixed << setprecision(4) << setw(8) << n_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << n_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << n_xyz(3) << endl;

		cout << fixed << setprecision(4) << setw(8) << ca_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << ca_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << ca_xyz(3) << endl;

		cout << fixed << setprecision(4) << setw(8) << c_xyz(1) << " ";
		cout << fixed << setprecision(4) << setw(8) << c_xyz(2) << " ";
		cout << fixed << setprecision(4) << setw(8) << c_xyz(3) << endl;

		cout << endl;
	}
	cout << endl;

	//make_random_perturbation(problem);
	solve_closure_problem(problem, solutions);
	apply_all_closure_solutions(solutions);
	//pick_random_solution(solutions, solution);

	for (int i = 1; i <= solutions.size(); i++) {
		cout << "Solution " << i << endl;
		cout << "==========" << endl;

		cout << "KIC XYZs" << endl; 
		cout << "--------" << endl;
		for (int j = 1; j <= solutions[i].atom_xyzs.size(); j++) {
			Coordinate xyz = solutions[i].atom_xyzs[j];
			cout << fixed << setprecision(4) << setw(8) << xyz[1] << " ";
			cout << fixed << setprecision(4) << setw(8) << xyz[2] << " ";
			cout << fixed << setprecision(4) << setw(8) << xyz[3] << endl;
			if (j % 3 == 0) cout << endl;
		}
		cout << endl;

		cout << "Pose XYZs" << endl; 
		cout << "---------" << endl;
		for (int j = 1; j <= solutions[i].pose.total_residue(); j++) {
			id::NamedAtomID n_id("N", j);
			id::NamedAtomID ca_id("CA", j);
			id::NamedAtomID c_id("C", j);

			core::PointPosition n_xyz  = solutions[i].pose.xyz(n_id);
			core::PointPosition ca_xyz = solutions[i].pose.xyz(ca_id);
			core::PointPosition c_xyz  = solutions[i].pose.xyz(c_id);

			cout << fixed << setprecision(4) << setw(8) << n_xyz(1) << " ";
			cout << fixed << setprecision(4) << setw(8) << n_xyz(2) << " ";
			cout << fixed << setprecision(4) << setw(8) << n_xyz(3) << endl;

			cout << fixed << setprecision(4) << setw(8) << ca_xyz(1) << " ";
			cout << fixed << setprecision(4) << setw(8) << ca_xyz(2) << " ";
			cout << fixed << setprecision(4) << setw(8) << ca_xyz(3) << endl;

			cout << fixed << setprecision(4) << setw(8) << c_xyz(1) << " ";
			cout << fixed << setprecision(4) << setw(8) << c_xyz(2) << " ";
			cout << fixed << setprecision(4) << setw(8) << c_xyz(3) << endl;

			cout << endl;
		}
		cout << endl;

		/*
		cout << "Solved Torsions" << endl; 
		cout << "===============" << endl;
		for (int j = 1; j <= solutions[i].torsion_angles.size(); j++) {
			Real angle = solutions[i].torsion_angles[j];
			angle = (angle < 0) ? angle + 360 : angle;
			cout << fixed << setprecision(4) << setw(8) << angle << endl;
			if (j % 3 == 0) cout << endl;
		}
		cout << endl;

		cout << "Applied Torsions" << endl; 
		cout << "================" << endl;
		for (int j = 1; j <= solutions[i].pose.total_residue(); j++) {
			Real phi = solutions[i].pose.phi(j);
			Real psi = solutions[i].pose.psi(j);
			Real omega = solutions[i].pose.omega(j);

			phi = (phi < 0) ? phi + 360 : phi;
			psi = (psi < 0) ? psi + 360 : psi;
			omega = (omega < 0) ? omega + 360 : omega;

			cout << fixed << setprecision(4) << setw(8) << phi << endl;
			cout << fixed << setprecision(4) << setw(8) << psi << endl;
			cout << fixed << setprecision(4) << setw(8) << omega << endl;
			cout << endl;
		}
		cout << endl;
		*/

		/*
		cout << "Angles  " << endl; 
		cout << "========" << endl;
		for (int j = 1; j <= solutions[i].bond_angles.size(); j++) {
			Real angle = solutions[i].bond_angles[j];
			angle = (angle < 0) ? angle + 360 : angle;
			cout << fixed << setprecision(4) << setw(8) << angle << endl;
			if (j % 3 == 0) cout << endl;
		}
		cout << endl;

		cout << "Lengths " << endl; 
		cout << "========" << endl;
		for (int j = 1; j <= solutions[i].bond_lengths.size(); j++) {
			Real length = solutions[i].bond_lengths[j];
			cout << fixed << setprecision(4) << setw(8) << length << endl;
			if (j % 3 == 0) cout << endl;
		}
		cout << endl;
		*/
	}

	for (int i = 1; i <= solutions.size(); i++) {
		string index = boost::lexical_cast<std::string>(i);
		solutions[i].pose.dump_pdb(output_pdb + index + ".pdb");
	}
}

void protocols_main(int argc, char** argv) { // {{{1
	chdir("../src/apps/pilot/kale/examples");
	NEW_OPT(kale::terminus, "Attempt to close the N-terminus.", false);

	devel::init(argc, argv);
	RG.set_seed(0);

	pose::Pose pose;
	string input_pdb, output_pdb;
	Size first_pivot, second_pivot, third_pivot;

	KinematicMoverOP mover = new KinematicMover();
	KinematicPerturberOP perturber = new NullKinematicPerturber();
	mover->set_perturber(perturber);

	if (option[OptionKeys::kale::terminus]() == true) {
		input_pdb = "structures/linear/4.helix.pdb";
		output_pdb = "terminus_closure.prc.pdb";
		first_pivot = 1; second_pivot = 2; third_pivot = 3;
	}
	else {
		input_pdb = "structures/linear/5.helix.pdb";
		output_pdb = "standard_closure.prc.pdb";
		first_pivot = 2; second_pivot = 3; third_pivot = 4;
	}

	cout << "Input: " << input_pdb << endl;
	cout << "Outputs: " << output_pdb << endl;
	cout << "Pivots: " << first_pivot << "/"
		                 << second_pivot << "/"
		                 << third_pivot << endl << endl;

	import_pose::pose_from_pdb(pose, input_pdb);

	mover->apply(pose);

	cout << "Applied Torsions" << endl; 
	cout << "================" << endl;
	for (int j = 1; j <= pose.total_residue(); j++) {
		Real phi = pose.phi(j);
		Real psi = pose.psi(j);
		Real omega = pose.omega(j);

		phi = (phi < 0) ? phi + 360 : phi;
		psi = (psi < 0) ? psi + 360 : psi;
		omega = (omega < 0) ? omega + 360 : omega;

		cout << fixed << setprecision(4) << setw(8) << phi << endl;
		cout << fixed << setprecision(4) << setw(8) << psi << endl;
		cout << fixed << setprecision(4) << setw(8) << omega << endl;
		cout << endl;
	}
	cout << endl;

	pose.dump_pdb(output_pdb);
}
// }}}1

int main(int argc, char** argv) {
	devel_main(argc, argv);
	//protocols_main(argc, argv);
}
