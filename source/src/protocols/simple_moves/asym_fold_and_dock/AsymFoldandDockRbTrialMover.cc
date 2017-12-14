// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockRbTrialMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/MonteCarlo.hh>

// Package headers
// options
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <utility>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/Tag.hh>
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockRbTrialMoverCreator.hh>


namespace protocols {
namespace simple_moves {
namespace asym_fold_and_dock {

static basic::Tracer TR( "protocols.simple_moves.asym_fold_and_dock.AsymFoldandDockRbTrialMover" );

AsymFoldandDockRbTrialMover::AsymFoldandDockRbTrialMover() :
	protocols::moves::Mover( "AsymFoldandDockRbTrialMover" ),
	smooth_move_(false), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

AsymFoldandDockRbTrialMover::AsymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn
) :
	protocols::moves::Mover( "AsymFoldandDockRbTrialMover" ),
	scorefxn_(std::move(scorefxn)), smooth_move_(false), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{}

AsymFoldandDockRbTrialMover::AsymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move
) :
	protocols::moves::Mover( "AsymFoldandDockRbTrialMover" ),
	scorefxn_(std::move(scorefxn)), smooth_move_(smooth_move), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{}

AsymFoldandDockRbTrialMover::AsymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move,
	core::Real rot_mag,
	core::Real trans_mag
) :
	protocols::moves::Mover( "AsymFoldandDockRbTrialMover" ),
	scorefxn_(std::move( scorefxn)), smooth_move_(smooth_move), rot_mag_(rot_mag), trans_mag_(trans_mag)
{}

void
AsymFoldandDockRbTrialMover::apply( core::pose::Pose & pose )
{

	using namespace basic::options;
	using namespace protocols::moves;

	TR.Debug << "Rb move applied..." << std::endl;

	core::Real trans_mag_smooth = 0.1;
	core::Real rot_mag_smooth = 1.0;

	// Docking options
	if ( smooth_move_ ) {
		if ( option[ OptionKeys::fold_and_dock::trans_mag_smooth ].user() ) {
			trans_mag_smooth = option[ OptionKeys::fold_and_dock::trans_mag_smooth ];
		}
		if ( option[ OptionKeys::fold_and_dock::rot_mag_smooth ].user() ) {
			rot_mag_smooth = option[ OptionKeys::fold_and_dock::rot_mag_smooth ];
		}
	}

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rb_rot_magnitude ].user() ) {
		rot_mag_ = ( option[ OptionKeys::fold_and_dock::rb_rot_magnitude ] );
	}

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rb_trans_magnitude ].user() ) {
		trans_mag_ = ( option[ OptionKeys::fold_and_dock::rb_trans_magnitude ] );
	}

	core::Real rot_mag_trial = smooth_move_ ? rot_mag_smooth : rot_mag_;
	core::Real trans_mag_trial = smooth_move_ ? trans_mag_smooth : trans_mag_;

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rigid_body_cycles ].user() ) {
		rigid_body_cycles_ = ( option[ OptionKeys::fold_and_dock::rigid_body_cycles ] );
	}

	if ( option[ OptionKeys::fold_and_dock::rigid_body_disable_mc ].user() ) {
		mc_filter_ = false;
	}

	// Setup Monte Carlo object
	protocols::moves::MonteCarloOP monteCarlo_( new protocols::moves::MonteCarlo(pose, *scorefxn_, 2.0 ) );


	// Setup the movemap
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_bb( false );
	movemap->set_chi( false );
	movemap->set_jump( true );

	//set up mover for docking
	protocols::rigid::RigidBodyPerturbMover rb_perturb =
		protocols::rigid::RigidBodyPerturbMover(
		pose, *movemap, rot_mag_trial, trans_mag_trial,
		rigid::partner_downstream, false );

	for ( Size i = 1; i <= rigid_body_cycles_; ++i ) {
		rb_perturb.apply( pose );
		if ( mc_filter_ ) monteCarlo_->boltzmann( pose );
	}
	//monteCarlo_->recover_low(pose);
}

moves::MoverOP
AsymFoldandDockRbTrialMover::clone() const
{
	return moves::MoverOP( new AsymFoldandDockRbTrialMover( *this ) );
}

moves::MoverOP
AsymFoldandDockRbTrialMover::fresh_instance() const
{
	return moves::MoverOP( new AsymFoldandDockRbTrialMover  );
}

void
AsymFoldandDockRbTrialMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	rot_mag_ = tag->getOption< core::Real >( "rot_mag", 8.0 );
	trans_mag_ = tag->getOption< core::Real >( "trans_mag", 3.0 ) ;
	rigid_body_cycles_ = tag->getOption< core::Size >( "rigid_body_cycles", 50 ) ;
	smooth_move_  = tag->getOption< bool >( "smooth_move", false ) ;
	mc_filter_ = tag->getOption< bool >( "mc_filter_", true ) ;
}

std::string
AsymFoldandDockRbTrialMover::get_name() const {
	return "AsymFoldandDockRbTrialMover";
}

std::string
AsymFoldandDockRbTrialMover::mover_name()
{
	return "AsymFoldandDockRbTrialMover";
}

void AsymFoldandDockRbTrialMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "rot_mag", xsct_real, "rotamer magnitude", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "trans_mag", xsct_real, "translation magnitude", "3.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rigid_body_cycles", xsct_non_negative_integer, "num of rigid body cycles", "50" )
		+ XMLSchemaAttribute::attribute_w_default( "smooth_move", xsct_rosetta_bool, "to use smooth move?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "mc_filter", xsct_rosetta_bool, "use Monte Carlo filter?", "true" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "a symmetric fold and dock", attlist  );
}

std::string AsymFoldandDockRbTrialMoverCreator::keyname() const {
	return AsymFoldandDockRbTrialMover::mover_name();
}

protocols::moves::MoverOP
AsymFoldandDockRbTrialMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AsymFoldandDockRbTrialMover );
}

void AsymFoldandDockRbTrialMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AsymFoldandDockRbTrialMover::provide_xml_schema( xsd );
}

} // asym_fold_and_dock
} // moves
} // protocols
