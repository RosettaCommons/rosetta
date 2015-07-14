// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Optimizes the membrane position given the high-res score function
/// @details	Optimizes the membrane position given the smooth high-res score
///				function; scans the center along the normal around the initial center
///				in 0.1A steps; scans the normal in 0.2degree steps along arches
///				over the x-axis, y-axis, xy-direction, -xy-direction; outcome is
///				deterministic
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeMembranePositionMover_cc
#define INCLUDED_protocols_membrane_OptimizeMembranePositionMover_cc

// Unit Headers
#include <protocols/membrane/OptimizeMembranePositionMover.hh> 
#include <protocols/membrane/OptimizeMembranePositionMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh> 
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/conformation/membrane/types.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.OptimizeMembranePositionMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
	
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: scorefxn = smooth2012
OptimizeMembranePositionMover::OptimizeMembranePositionMover() : protocols::moves::Mover()
{
	set_defaults();
	register_options();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
OptimizeMembranePositionMover::OptimizeMembranePositionMover( OptimizeMembranePositionMover const & src ) : protocols::moves::Mover( src ),
	ft_( src.ft_ ),
	sfxn_( src.sfxn_ ),
	score_best_( src.score_best_ ),
	starting_z_( src.starting_z_ ),
	best_z_( src.best_z_ ),
	stepsize_z_( src.stepsize_z_ ),
	best_center_( src.best_center_ ),
	stepsize_angle_( src.stepsize_angle_ ),
	best_normal_( src.best_normal_ )
{}

/// @brief Assignment Operator
OptimizeMembranePositionMover & OptimizeMembranePositionMover::operator = ( OptimizeMembranePositionMover const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
		
	// Otherwise, create a new object
	return *( new OptimizeMembranePositionMover( *this ) );
}

/// @brief Destructor
OptimizeMembranePositionMover::~OptimizeMembranePositionMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
OptimizeMembranePositionMover::clone() const {
	return ( protocols::moves::MoverOP( new OptimizeMembranePositionMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
OptimizeMembranePositionMover::fresh_instance() const {
	return protocols::moves::MoverOP( new OptimizeMembranePositionMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
OptimizeMembranePositionMover::parse_my_tag(
	 utility::tag::TagCOP /*tag*/,
	 basic::datacache::DataMap &,
	 protocols::filters::Filters_map const &,
	 protocols::moves::Movers_map const &,
	 core::pose::Pose const &
	 ) {

	// TODO: implement this
	
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
OptimizeMembranePositionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new OptimizeMembranePositionMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
OptimizeMembranePositionMoverCreator::keyname() const {
	return OptimizeMembranePositionMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
OptimizeMembranePositionMoverCreator::mover_name() {
	return "OptimizeMembranePositionMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (OptimizeMembranePositionMover)
std::string
OptimizeMembranePositionMover::get_name() const {
	return "OptimizeMembranePositionMover";
}

/// @brief Flip the downstream partner in the membrane
void OptimizeMembranePositionMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::conformation::membrane;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;
	
	TR << "Optimizing membrane position using scorefunction " << sfxn_ << std::endl;
	
	// remember original foldtree
	ft_ = pose.fold_tree();
	
	// reorder foldtree with residue 1 as root of the foldtree
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	foldtree.reorder( 1 );
	pose.fold_tree( foldtree );
	
	// show foldtree
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);

	// compute initial membrane position from structure
	MembranePositionFromTopologyMoverOP initmem( new MembranePositionFromTopologyMover( true ) );
	initmem->apply( pose );
	
	// show foldtree
	TR << "foldtree after MembranePositionFromTopologyMover" << std::endl;
	pose.fold_tree().show(std::cout);
	
	// optimize membrane center
	optimize_membrane_center( pose );
	
	// optimize membrane normal
	optimize_membrane_normal( pose );
	
	// reset foldtree to original one
	pose.fold_tree( ft_ );
	TR << "foldtree reset to original one" << std::endl;
	pose.fold_tree().show(std::cout);
	
}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void OptimizeMembranePositionMover::register_options() {
	
	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
	
}

/// @brief Set default values
void OptimizeMembranePositionMover::set_defaults() {

	using namespace core::scoring;

	// create scorefxn
	sfxn_ = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
	
	// starting z and angle and their stepsizes
	starting_z_ = -10.0;
	stepsize_z_ = 0.1;
	stepsize_angle_ = 0.5;
	
	// set best score
	score_best_ = 999999;
	best_z_ = starting_z_;
	
}// set_defaults

////////////////////////////////////////////////////////////////////////////////

/// @brief Optimize membrane center
void OptimizeMembranePositionMover::optimize_membrane_center( Pose & pose ) {

	using namespace numeric;

	// starting membrane center
	TR << "Optimizing membrane center" << std::endl;
	TR << "Starting MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
	
	// get components of center vector
	core::Vector center( pose.conformation().membrane_info()->membrane_center() );
	
	// initialize new center
	core::Vector new_center;
	core::Real new_z( static_cast< core::Real >( starting_z_ ) );
	
	// set new starting center
	new_center.assign( 0, 0, new_z );
	
	// SetMembraneCenterMover
	SetMembraneCenterMoverOP set_center( new SetMembraneCenterMover( new_center ) );
	set_center->apply( pose );
	
	// score the pose
	core::Real score_old = sfxn_->score( pose );
	core::Real score_new( score_old );
	score_best_ = score_old;
	
	// get iterations
	Size iterations = nearest_size( abs_difference( starting_z_, 0.0 ) * 2 / stepsize_z_ );
	
	// drag along the z axis and score
	for ( Size i = 1; i <= iterations; ++i ){
		
		// drag z through the membrane
		new_z += stepsize_z_;
		
		// set new center
		new_center.assign( 0, 0, new_z );
		
		// SetMembraneCenterMover
		SetMembraneCenterMoverOP set_center1( new SetMembraneCenterMover( new_center ) );
		set_center1->apply( pose );
		
		// rescore pose
		score_new = sfxn_->score( pose );
		
		// save best center
		if ( score_new < score_best_ ){
			
			best_z_ = new_z;
			score_best_ = score_new;
//			TR << "best center: " << new_center.to_string() << " with score " << score_new << std::endl;
		}
		
		// set old score to this one
		score_old = score_new;
	}
	
	// Set center to best one that we found
	best_center_.assign( 0, 0, best_z_ );
	TR << "Best center: " << best_center_.to_string() << std::endl;
	
	// SetMembraneCenterMover
	SetMembraneCenterMoverOP set_center2( new SetMembraneCenterMover( best_center_ ) );
	set_center2->apply( pose );
	
	// final score from center search
	score_best_ = sfxn_->score( pose );
	TR << "Final score from center search: " << score_best_ << std::endl;
	
	// final membrane center
	best_center_ = pose.conformation().membrane_info()->membrane_center();
	TR << "Final MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;

} // optimize membrane center

////////////////////////////////////////////////////////////////////////////////
	
/// @brief Optimize membrane normal
void OptimizeMembranePositionMover::optimize_membrane_normal( Pose & pose ) {
	
	// starting membrane normal
	TR << "Optimizing membrane normal" << std::endl;
	TR << "Starting MemInfo normal: " << pose.conformation().membrane_info()->membrane_normal().to_string() << std::endl;
	
	// initialize best normal
	best_normal_ = pose.conformation().membrane_info()->membrane_normal();
	core::Real angle( 0 );
	core::Real new_x, new_y, new_z;
	core::Vector new_normal;

	// score the pose
//	core::Real score_old = score_best_;
	core::Real score_new = score_best_;
	
	// sample an arch over several different axes / directions
	for ( Size j = 1; j <= 4; ++j ){
		
		angle = 45;
		
		// sample angles
		for ( Size i = 1; i <= 180; ++i ){
			
			// increase angle by 1 degree
			angle += stepsize_angle_;
			
			// for first iteration, sample arch over x-axis
			if ( j == 1 ) {
				new_x = 0;
				new_y = 15 * cos( numeric::conversions::radians( angle ) );
				new_z = 15 * sin( numeric::conversions::radians( angle ) );
			}
			
			// for second iteration, sample arch over y-axis
			else if ( j == 2 ) {
				new_x = 15 * cos( numeric::conversions::radians( angle ) );
				new_y = 0;
				new_z = 15 * sin( numeric::conversions::radians( angle ) );
			}
			
			// for third iteration, sample arch over x=y
			else if ( j == 3 ) {
				if ( i == 1 ) { angle = -45; }
				new_x = 15 * sin( numeric::conversions::radians( angle ) );
				new_y = 15 * sin( numeric::conversions::radians( angle ) );
				new_z = 15 * cos( numeric::conversions::radians( angle ) );
			}
			
			// for fourth iteration, sample arch over x=-y
			else if ( j == 4 ) {
				if ( i == 1 ) { angle = -45; }
				new_x = 15 * sin( numeric::conversions::radians( angle ) );
				new_y = -15 * sin( numeric::conversions::radians( angle ) );
				new_z = 15 * cos( numeric::conversions::radians( angle ) );
			}
			
			// set new normal
			new_normal.assign( new_x, new_y, new_z );
			new_normal.normalize( 15 );
			
			// SetMembraneNormalMover
			SetMembraneNormalMoverOP set_normal( new SetMembraneNormalMover( new_normal ) );
			set_normal->apply( pose );
			
			// rescore pose
			score_new = sfxn_->score( pose );
			
			// save best normal
			if ( score_new < score_best_ ){
				
				best_normal_ = new_normal;
				score_best_ = score_new;
//				TR << "best normal: " << best_normal_.to_string() << " with score " << score_best_ << std::endl;
			}
			
			// set old score to this one
//			score_old = score_new;
		}
	}
	
	// Set normal to best one that we found
	TR << "Best normal: " << best_normal_.to_string() << std::endl;
	
	// SetMembraneNormalMover
	SetMembranePositionMoverOP set_final( new SetMembranePositionMover( best_center_, best_normal_ ) );
	set_final->apply( pose );
	
	// final score from center search
	score_new = sfxn_->score( pose );
	TR << "Final score from normal search: " << score_new << std::endl;
	
	// final membrane info
	TR << "Final MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
	TR << "Final MemInfo normal: " << pose.conformation().membrane_info()->membrane_normal().to_string() << std::endl;
	
} // optimize membrane normal


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeMembranePositionMover_cc
