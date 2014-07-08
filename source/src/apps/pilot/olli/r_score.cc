// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/moves/Mover.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>


// Utility headers
#include <basic/options/option_macros.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/excn/EXCN_Base.hh>
#include <utility/io/mpistream.hh>



static basic::Tracer tr("main");

using namespace core;
using namespace protocols;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;




void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  OPT( in::file::s );
}

// Forward
class ScoreMover;

// Types
typedef  utility::pointer::owning_ptr< ScoreMover >  ScoreMoverOP;
typedef  utility::pointer::owning_ptr< ScoreMover const >  ScoreMoverCOP;

class ScoreMover : public moves::Mover {
public:
	ScoreMover( scoring::ScoreFunctionOP score ) : score_( score ){};
	virtual void apply( core::pose::Pose& );

private:
	scoring::ScoreFunctionOP score_;
};

void ScoreMover::apply( core::pose::Pose &pose ) {
	std::string fname( jd2::current_output_name() );
	Real score ( (* score_)( pose ));
	tr.Info << fname << " " << score << std::endl;
}



void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	ScoreMoverOP tool =  new ScoreMover( scoring::get_score_function() );
	protocols::jd2::JobDistributor::get_instance()->go( tool );


	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	register_options();
	devel::init( argc, argv );

	try{
		run();
	} catch ( utility::excn::EXCN_Base& excn ) {
		excn.show( std::cerr );
	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


