/// @file
/// @brief

#include <devel/init.hh>

#include <protocols/comparative_modeling/hybridize/HybridizeProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	try{
		protocols::jd2::JobDistributor::get_instance()->go( new protocols::comparative_modeling::hybridize::HybridizeProtocol() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}


