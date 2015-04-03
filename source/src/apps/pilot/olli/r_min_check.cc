// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

// #include <protocols/abinitio/Templates.hh>
// #include <protocols/abinitio/TemplateJumpSetup.hh>
// #include <protocols/abinitio/PairingStatistics.hh>
// #include <protocols/abinitio/StrandConstraints.hh>
// #include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
//#include <core/pose/util.hh>
#include <devel/init.hh>
//#include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
//#include <core/scoring/constraints/BoundConstraint.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/util.hh>

//for derivative check
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/loops/Loops.hh>
#include <core/scoring/constraints/util.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
//#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/exit.hh>

//Auto Headers
#include <utility/excn/EXCN_Base.hh>
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>


static thread_local basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
//using namespace abinitio;
//using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
using namespace ObjexxFCL::format;

OPT_1GRP_KEY( Boolean, score_app, linmin )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  //  Templates::register_options();
  OPT( in::file::s );
  NEW_OPT( score_app::linmin, "Do a quick round of linmin before reporting the score", false );
}

// Forward
class MinToolMover;

// Types
typedef  utility::pointer::owning_ptr< MinToolMover >  MinToolMoverOP;
typedef  utility::pointer::owning_ptr< MinToolMover const >  MinToolMoverCOP;

class MinToolMover : public moves::Mover {
public:
  virtual std::string get_name() const {return "mintool mover";}
  MinToolMover() {};
  virtual void apply( core::pose::Pose& );
private:
  Size nstruct_;
};

void MinToolMover::apply( core::pose::Pose &pose ) {

  scoring::ScoreFunction score;
  std::string fname( jd2::current_output_name() );

  core::scoring::constraints::add_constraints_from_cmdline( pose, score);
  score.show( tr.Info, pose );
  tr.Info << fname << " " << score( pose ) << std::endl;
  ++nstruct_;

  core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  movemap->set_bb( true ); movemap->set_chi( true );
  protocols::simple_moves::MinMoverOP minmover =
    new protocols::simple_moves::MinMover(
        movemap, score.clone(), "linmin", 1e-4,
	false /*use_nblist*/, true /*deriv_check*/, true /*verbose driv check*/
    );
  minmover->apply( pose );

  tr.Info << fname << " " << score( pose ) << std::endl;

}


void run() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  MinToolMoverOP cst_tool =  new MinToolMover;
  protocols::jd2::JobDistributor::get_instance()->go( cst_tool, new jd2::NoOutputJobOutputter );

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

  //	try{
  run();
  //	} catch ( utility::excn::EXCN_Base& anExcn ) {
  //		anExcn.show( std::cerr );
  //	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}


