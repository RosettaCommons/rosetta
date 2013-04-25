// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
//     vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file SurfaceOrientMover.cc
/// @author Robin A Thottungal (rathottungal@gmail.com)
/**

**/

// Unit Headers
#include <protocols/surface_docking/SlideIntoSurface.hh>
#include <protocols/surface_docking/SurfaceOrientMover.fwd.hh>

// Package Headers

// Project headers
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// for adding data to pose
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

//Utility Headers
#include <utility/exit.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.SurfaceDocking.SlideIntoSurface");

namespace protocols {
namespace surface_docking {

using namespace protocols;
using namespace core;
using namespace core::scoring;
using namespace protocols::moves;


// default constructor
SlideIntoSurface::SlideIntoSurface() : Mover()
{
	Mover::type( "SlideIntoSurface" );
	rb_jump_ = 1;
	scorefxn_ = ScoreFunctionFactory::create_score_function
                                            ( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//constructor
SlideIntoSurface::SlideIntoSurface(
	core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump)
{
	Mover::type( "DockingSlideIntoContact" );
	scorefxn_ = ScoreFunctionFactory::create_score_function
                                            ( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//destructor
SlideIntoSurface::~SlideIntoSurface() {}


void SlideIntoSurface::apply( core::pose::Pose & pose )
{
	using namespace moves;
	/***
	// Extracting surfaceVector info
	std::string SurfVectors[3];
	for(Size i=0; i<pose.pdb_info()->remarks().size(); i++) {
		TR << pose.pdb_info()->remarks().at(i).num << " " << pose.pdb_info()->remarks().at(i).value << std::endl;
		SurfVectors[i]=pose.pdb_info()->remarks().at(i).value;
	}
    surface_docking::SurfaceParametersOP surfaceParams=
            new surface_docking::SurfaceParameters(SurfVectors[0],SurfVectors[1],SurfVectors[2]);
	**/
	surface_docking::SurfaceParameters & surfaceParams=
				*( static_cast< surface_docking::SurfaceParameters * >
							( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SURFACE_PARAMS )() ));
	rigid::RigidBodyTransMoverOP mover( new rigid::RigidBodyTransMover(pose,rb_jump_));
	//TR<<"SlideAxis Vector:"<<surfaceParams.slideaxis<<std::endl;
	mover->trans_axis(-1*surfaceParams.slideaxis);
	( *scorefxn_ )( pose );
	TR<<"Initial Score:"<<pose.energies().total_energies()
                                        [scoring::interchain_vdw ]<<std::endl;
	TR << "Moving away" << std::endl;
	core::Size const counter_breakpoint( 500 );
	core::Size counter( 0 );
	// first try moving away from each other
	while ( pose.energies().total_energies()[ scoring::interchain_vdw ] > 0.1
                                            && counter <= counter_breakpoint ) {
		mover->apply( pose );
		( *scorefxn_ )( pose );
		TR<<"current score:"<<pose.energies().total_energies()
                                    [ scoring::interchain_vdw ]<<std::endl;
		++counter;
	}
	if( counter > counter_breakpoint ){
		TR<<"failed moving away with original vector.Aborting SlideIntoSurface"
                                                                    <<std::endl;
		set_current_tag( "fail" );
		return;
	}
	counter = 0;
	// then try moving towards each other
	TR << "Moving together" << std::endl;
	//mover->trans_axis().negate();
	while ( counter <= counter_breakpoint && pose.energies().total_energies()
                                        [ scoring::interchain_vdw ] < 0.1 ) {
		mover->apply( pose );
		//TR<<"current score:"<<pose.energies().total_energies()
                                       //[ scoring::interchain_vdw ]<<std::endl;
		( *scorefxn_ )( pose );
		TR<<"Initial Score:"<<pose.energies().total_energies()<<std::endl;
		++counter;
	}
	if( counter > counter_breakpoint ){
		TR<<"moving together failed. Aborting SlideIntoSurface"<<std::endl;
		set_current_tag( "fail" );
		return;
	}
	// move away again until just touching
	mover->trans_axis().negate();
	mover->apply( pose );
}

std::string SlideIntoSurface::get_name() const {
	return "SlideIntoSurface";
}

/////////////////////////////// FaSlideIntoSurface/////////////////////////////
// default constructor
FaSlideIntoSurface::FaSlideIntoSurface()
{

	Mover::type( "FaSlideIntoSurface" );
	rb_jump_ = 1;
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}


//constructor
FaSlideIntoSurface::FaSlideIntoSurface(
	core::Size const rb_jump
	//protocols::surface_docking::SurfaceParametersOP surfParams
) : Mover(), rb_jump_(rb_jump), tolerance_(0.2)
{
	Mover::type( "FaSlideIntoSurface" );
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

//destructor
FaSlideIntoSurface::~FaSlideIntoSurface() {}

void FaSlideIntoSurface::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	// Reading datacache for surface Vectors
	//scoring::SurfaceParameters & surfaceParams =
	//		*( static_cast< core::scoring::SurfaceParameters * >
	//					( pose.data().get_ptr( SURFACE_PARAMS )() ));
	/***
	// Extracting surfaceVector info
		std::string SurfVectors[3];
		for(Size i=0; i<pose.pdb_info()->remarks().size(); i++) {
			TR << pose.pdb_info()->remarks().at(i).num << " " << pose.pdb_info()->remarks().at(i).value << std::endl;
			SurfVectors[i]=pose.pdb_info()->remarks().at(i).value;
		}
	surface_docking::SurfaceParametersOP surfaceParams=
	          new surface_docking::SurfaceParameters(SurfVectors[0],SurfVectors[1],SurfVectors[2]);
	**///does not work!
	surface_docking::SurfaceParameters & surfaceParams=
					*( static_cast< surface_docking::SurfaceParameters * >
								( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SURFACE_PARAMS )() ));
	rigid::RigidBodyTransMoverOP trans_mover( new rigid::RigidBodyTransMover(pose,rb_jump_));
	TR<<"SlideAxis Vector:"<<surfaceParams.slideaxis<<std::endl;
	trans_mover->trans_axis(-1*surfaceParams.slideaxis);

	(*scorefxn_)( pose );
	core::Real const initial_fa_rep = pose.energies().total_energies()[ fa_rep ];
	bool are_touching = false;
	//moves::RigidBodyTransMover trans_mover( pose, rb_jump_ );

	//int i=1;
	// Take 2A steps till clash, then back apart one step.  Now you're within 2A of touching.
	// Repeat with 1A steps, 0.5A steps, 0.25A steps, etc until you're as close are you want.
	for( core::Real stepsize = 2.0; stepsize > tolerance_; stepsize /= 2.0 ) {
		trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move together
		trans_mover->step_size(stepsize);
		core::Size const counter_breakpoint( 500 );
		core::Size counter( 0 );
		do
		{
			trans_mover->apply( pose );
			(*scorefxn_)( pose );
			core::Real const push_together_fa_rep = pose.energies().total_energies()[ fa_rep ];
			//std::cout << "fa_rep = " << push_together_fa_rep << std::endl;
			are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
			//std::ostringstream s;
			//s << "snapshot" << i << ".pdb";
			//pose.dump_pdb(s.str());
			//i += 1;
			++counter;
		} while( counter <= counter_breakpoint && !are_touching );
		if( counter > counter_breakpoint ){
			TR<<"Failed Fadocking Slide Together. Aborting."<<std::endl;
			set_current_tag( "fail" );
		}
		trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move apart
		trans_mover->apply( pose );
	}
}

std::string
FaSlideIntoSurface::get_name() const {
	return "FaSlideIntoSurface";
}



/////////////////////////////// FaSlideAwayFromSurface/////////////////////////////
// default constructor
FaSlideAwayFromSurface::FaSlideAwayFromSurface()
{

	Mover::type( "FaSlideIntoSurface" );
	rb_jump_ = 1;
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}


//constructor
FaSlideAwayFromSurface::FaSlideAwayFromSurface(
	core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump), tolerance_(0.2)
{
	Mover::type( "FaSlideIntoSurface" );
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

//destructor
FaSlideAwayFromSurface::~FaSlideAwayFromSurface() {}

void FaSlideAwayFromSurface::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	// Reading datacache for surface Vectors
	surface_docking::SurfaceParameters & surfaceParams=
					*( static_cast< surface_docking::SurfaceParameters * >
								( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SURFACE_PARAMS )() ));
    rigid::RigidBodyTransMoverOP trans_mover( new rigid::RigidBodyTransMover(pose,rb_jump_));
	TR<<"SlideAxis Vector:"<<surfaceParams.slideaxis<<std::endl;
	trans_mover->trans_axis(-1*surfaceParams.slideaxis);

	(*scorefxn_)( pose );
	//core::Real const initial_fa_rep = pose.energies().total_energies()[ fa_rep ];
	//bool are_touching = false;
	//moves::RigidBodyTransMover trans_mover( pose, rb_jump_ );

	//Moving away
	TR<<"Moving Away"<<std::endl;
	trans_mover->trans_axis( trans_mover->trans_axis() ); // now move away
	trans_mover->step_size(100);
	trans_mover->apply( pose );

	/***
	//int i=1;
	// Take 2A steps till clash, then back apart one step.  Now you're within 2A of touching.
	// Repeat with 1A steps, 0.5A steps, 0.25A steps, etc until you're as close are you want.
	for( core::Real stepsize = 2.0; stepsize > tolerance_; stepsize /= 2.0 ) {
		trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move together
		trans_mover->step_size(stepsize);
		core::Size const counter_breakpoint( 500 );
		core::Size counter( 0 );
		do
		{
			trans_mover->apply( pose );
			(*scorefxn_)( pose );
			core::Real const push_together_fa_rep = pose.energies().total_energies()[ fa_rep ];
			//std::cout << "fa_rep = " << push_together_fa_rep << std::endl;
			are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
			//std::ostringstream s;
			//s << "snapshot" << i << ".pdb";
			//pose.dump_pdb(s.str());
			//i += 1;
			++counter;
		} while( counter <= counter_breakpoint && !are_touching );
		if( counter > counter_breakpoint ){
			TR<<"Failed Fadocking Slide Together. Aborting."<<std::endl;
			set_current_tag( "fail" );
		}
		trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move apart
		trans_mover->apply( pose );
	}
	**/
}

std::string
FaSlideAwayFromSurface::get_name() const {
	return "FaSlideAwayFromSurface";
}

}	//surfaceDockingProtocol

}	//protocol
