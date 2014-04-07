/// @file
/// @brief

#include <devel/init.hh>

#include <protocols/comparative_modeling/hybridize/HybridizeProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>



///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try{
		protocols::jd2::JobDistributor::get_instance()->go(
			new protocols::comparative_modeling::hybridize::HybridizeProtocol( )
		);
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
    devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
    }


