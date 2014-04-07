/// @file
/// @brief


#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <devel/init.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <utility/excn/Exceptions.hh>



///////////////////////////////////////////////////////////////////////////////

class CustomMover : public protocols::moves::Mover {
public:
	CustomMover(){}
	void apply( Pose & pose) {
		numeric::xyzVector< core::Real> X = core::pose::symmetry::get_symm_axis( pose );
		std::cerr << "SYMM AXIS : " << X[0] << "," << X[1] << "," << X[2] << std::endl;
	}

	virtual std::string get_name() const {
		return "CustomMover";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( new SetupForSymmetryMover() );
	seq->add_mover( new CustomMover() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
    try {
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
    return 0;
    }
