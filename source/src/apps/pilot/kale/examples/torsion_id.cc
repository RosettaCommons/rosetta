#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

using namespace std;
using core::pose::Pose;
using core::id::TorsionID;
using core::id::BB;
using core::import_pose::pose_from_file;

int main(int argc, char** argv) {
	devel::init(argc, argv);
	Pose pose;

	pose_from_file(pose, "structures/linear/6.symmetry.pdb", core::import_pose::PDB_file);

	cout << "Phi:   " << pose.phi(1) << endl;
	cout << "Psi:   " << pose.psi(1) << endl;
	cout << "Omega: " << pose.omega(1) << endl;

	cout << endl;

	cout << "BB 1:  " << pose.torsion(TorsionID(1, BB, 1)) << endl;
	cout << "BB 2:  " << pose.torsion(TorsionID(1, BB, 2)) << endl;
	cout << "BB 3:  " << pose.torsion(TorsionID(1, BB, 3)) << endl;
}
