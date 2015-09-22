#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale.examples" );

using namespace core;

int main(int argc, char **argv) {

	devel::init(argc, argv);

	pose::Pose pose;
	import_pose::pose_from_pdb(pose, "structures/marked_loop.8.pdb");

	for (Size i = 1; i <= pose.n_residue(); i++) {
		conformation::Residue residue = pose.residue(i);

		for (Size j = 1; j <= residue.natoms(); j++) {
			TR << "Atom ID: " << j << std::endl;
			residue.atom_type(j).print(TR);
			TR << std::endl;
		}
	}

}
