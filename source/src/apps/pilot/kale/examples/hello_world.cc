#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

using namespace std;
using core::pose::Pose;
using core::import_pose::pose_from_file;

int main(int argc, char** argv) {
	devel::init(argc, argv);
	Pose pose;

	pose_from_file(pose, "structures/linear/6.symmetry.pdb", core::import_pose::PDB_file);
	pose.dump_pdb("pdb_from_pose.pdb");
}
