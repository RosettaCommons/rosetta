// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/doug/PeptoidDihedralGrabber.cc
/// @brief Simply prints backbone and side chain dihedral angles of a peptoid
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// protocols header
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueType.hh>

// devel headers
#include <devel/init.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cyclization.OptionKeys.gen.hh>

// c++
#include <sstream>

// tracer
static thread_local basic::Tracer TR( "PeptoidDihedralGrabber" );

// local options
basic::options::BooleanOptionKey const cyclic( "cyclic" );

// super simple class to grab and print stuff
class PeptoidDihedralGrabber : public protocols::moves::Mover {
public:
// ctor
PeptoidDihedralGrabber( bool cyclic );

//dtor
virtual ~PeptoidDihedralGrabber(){}

// mover interface
virtual void apply( core::pose::Pose & pose );
virtual std::string get_name() const { return "PeptoidDihedralGrabber"; }
virtual protocols::moves::MoverOP clone() const { return new PeptoidDihedralGrabber( *this ); }
virtual protocols::moves::MoverOP fresh_instance() const { return clone(); }

private:
bool cyclic_;

};

PeptoidDihedralGrabber::PeptoidDihedralGrabber( bool cyclic ) :
  cyclic_( cyclic )
{}

void
PeptoidDihedralGrabber::apply( core::pose::Pose & pose )
{
	using namespace core;
	using namespace pose;
	using namespace conformation;
	using namespace chemical;

	for( Size i(1); i <= pose.total_residue(); ++i ) {
		// print name
		TR << "aa: " << pose.residue( i ).type().name3() << ", ";

		// print preceding omg, phi psi
		// first residue
		if ( i == 1 ) {
			if ( cyclic_ ) {
				TR << "omg: " << pose.omega( pose.total_residue() ) << ", phi: " << pose.phi( i ) << ", psi: " << pose.psi( i ) << ", ";
			} else {
				TR << "omg: " << 0.0 << ", phi: " << pose.phi( i ) << ", psi: " << pose.psi( i ) << ", ";
			}
		// last residue
		} else if ( i == pose.total_residue() ) {
			TR << "omg: " << pose.omega( i - 1 ) << ", phi: " << pose.phi( i - 1 ) << ", psi: " << pose.psi( i - 1 ) << ", ";
		// everything else
		} else {
			TR << "omg: " << pose.omega( i - 1 ) << ", phi: " << pose.phi( i - 1 ) << ", psi: " << pose.psi( i - 1 ) << ", ";
		}

		// print sidechain info
		for( Size j(1); j <= pose.residue( i ).type().nchi(); ++j ) {
			std::stringstream chi_string;
			chi_string << "x" << j << ": ";
			TR << chi_string.str() << pose.residue( i ).chi( j );
		}
		TR << std::endl;
	}
}

// typedefs
typedef utility::pointer::owning_ptr< PeptoidDihedralGrabber > PeptoidDihedralGrabberOP;

int
main( int argc, char * argv [] )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys::cyclization;
  using namespace protocols::simple_moves;
  using namespace protocols::moves;

  // add local options
 	option.add( cyclic, "cyclic" ).def("True");

  // init
  devel::init( argc, argv );

	// setup sequence mover
	SequenceMoverOP sm( new SequenceMover() );

  // setup the cyclization mover(s) ( just add patches and constraints don't minimize )
	if ( option[chains_to_cyclize].user() ) {
		core::Size num_cyclic_chains( option[chains_to_cyclize].value().size() );
		for ( core::Size i(1); i <= num_cyclic_chains; ++i ) {
			sm->add_mover( new CyclizationMover( option[chains_to_cyclize].value()[i], true, false, 0 ) );
		}
	}

	// setup peptoid dihedral grabber mover
	PeptoidDihedralGrabberOP pdg( new PeptoidDihedralGrabber( option[cyclic].value() ) );

	// setup a task factory
	core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
	tf->push_back( new core::pack::task::operation::InitializeFromCommandline() );
	tf->push_back( new core::pack::task::operation::RestrictResidueToRepacking() );
	//tf->push_back( new core::pack::task::operation::ReadResfile() );

	// setup score function
	core::scoring::ScoreFunctionOP scrfxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );

	// setup pack rotamers movie
	PackRotamersMoverOP pm( new protocols::simple_moves::PackRotamersMover() );
	pm->task_factory( tf );
	pm->score_function( scrfxn );

	// add movers
	sm->add_mover( pdg );
	sm->add_mover( pm );
	sm->add_mover( pdg );

  // go go go
  protocols::jd2::JobDistributor::get_instance()->go( sm );

  TR << "\n+-----------------------------------------------------------------+\n"
     <<   "|                              DONE                               |\n"
     <<   "+-----------------------------------------------------------------+" << std::endl;

  return 0;
}
