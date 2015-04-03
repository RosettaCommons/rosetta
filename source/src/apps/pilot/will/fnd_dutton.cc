#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/abinitio/BrokerMain.hh>


int main(int argc, char *argv[]) {

	try {

	protocols::abinitio::register_options_broker();

	devel::init(argc,argv);
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::simple_moves::symmetry;
	using namespace protocols::symmetric_docking;

	SetupForSymmetryMoverOP setup_mover = new SetupForSymmetryMover;
	SymDockProtocolOP dock_mover = new SymDockProtocol;
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
	seq_mover->add_mover( setup_mover );
	seq_mover->add_mover( dock_mover );
	protocols::jd2::JobDistributor::get_instance()->go( seq_mover );

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

