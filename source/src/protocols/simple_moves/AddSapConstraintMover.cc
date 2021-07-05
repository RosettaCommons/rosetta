// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddSapConstraintMover.cc
/// @brief Mover that adds the SapConstraint to the pose
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <protocols/simple_moves/AddSapConstraintMover.hh>
#include <protocols/simple_moves/AddSapConstraintMoverCreator.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraint.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/util.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/pointer/memory.hh>

// C++ Headers
#include <iostream>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace core::pack::guidance_scoreterms::sap;

// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.simple_moves.AddSapConstraintMover" );


// AddSapConstraintMover; based on the protocols::moves::Mover basis class
AddSapConstraintMover::AddSapConstraintMover() :
	protocols::moves::Mover("AddSapConstraintMover"),
	options_( utility::pointer::make_shared< core::pack::guidance_scoreterms::sap::SapConstraintOptions >() )
{
	set_speed( "fast" );
}


AddSapConstraintMover::AddSapConstraintMover( AddSapConstraintMover const & ot ) :
	protocols::moves::Mover( ot )
{
	*this = ot;
}

AddSapConstraintMover &
AddSapConstraintMover::operator=( AddSapConstraintMover const & ot ) {
	protocols::moves::Mover::operator=( ot );
	runtime_assert( ot.options_ );
	options_ = ot.options_->clone();
	return *this;
}

void
AddSapConstraintMover::apply( core::pose::Pose & pose )
{
	options_->sanity_check();
	SapConstraintOP cst = utility::pointer::make_shared<SapConstraint>( options_ );
	pose.add_constraint( cst );
}

moves::MoverOP
AddSapConstraintMover::clone() const
{
	return utility::pointer::make_shared< AddSapConstraintMover >( *this );
}

moves::MoverOP
AddSapConstraintMover::fresh_instance() const
{
	return utility::pointer::make_shared< AddSapConstraintMover >( );
}

void AddSapConstraintMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	if ( tag->hasOption( "name" ) ) {
		set_name( tag->getOption< std::string >( "name" ) );
	}

	if ( tag->hasOption( "speed" ) ) {
		set_speed( tag->getOption< std::string >( "speed" ) );
	}

	if ( tag->hasOption( "sap_goal" ) ) {
		set_sap_goal( tag->getOption< core::Real >( "sap_goal" ) );
	}
	if ( tag->hasOption( "sap_lb_goal" ) ) {
		set_sap_lb_goal( tag->getOption< core::Real >( "sap_lb_goal" ) );
	}
	if ( tag->hasOption( "packing_correction" ) ) {
		set_packing_correction( tag->getOption< core::Real >( "packing_correction" ) );
	}
	if ( tag->hasOption( "penalty_per_sap" ) ) {
		set_penalty_per_sap( tag->getOption< core::Real >( "penalty_per_sap" ) );
	}
	if ( tag->hasOption( "full_accuracy_when_scoring" ) ) {
		set_full_accuracy_when_scoring( tag->getOption< bool >( "full_accuracy_when_scoring" ) );
	}

	if ( tag->hasOption( "score_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data, "score_selector" ) );
		if ( selector != nullptr ) {
			options_->score_selector( selector );
		}
	}

	if ( tag->hasOption( "sap_calculate_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data, "sap_calculate_selector" ) );
		options_->sap_calculate_selector( selector );
	}

	if ( tag->hasOption( "sasa_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data, "sap_calculate_selector" ) );
		options_->sasa_selector( selector );
	}
}


void
AddSapConstraintMover::set_speed( std::string const & speed ) {
	if ( speed == "slow" ) {
		options_->fast(false);
		options_->lightning(false);
	} else if ( speed == "fast" ) {
		options_->fast(true);
		options_->lightning(false);
	} else if ( speed == "lightning" ) {
		options_->fast(false);
		options_->lightning(true);
	} else {
		utility_exit_with_message("AddSapConstraintMover: Unknown speed setting: " + speed);
	}
}

void
AddSapConstraintMover::set_score_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in ) {
	options_->score_selector( selector_in );
}

void
AddSapConstraintMover::set_sap_calculate_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in ) {
	options_->sap_calculate_selector( selector_in );
}

void
AddSapConstraintMover::set_sasa_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in ) {
	options_->sasa_selector( selector_in );
}

void
AddSapConstraintMover::set_sap_goal( core::Real goal ) {
	options_->sap_goal( goal );
}

void
AddSapConstraintMover::set_sap_lb_goal( core::Real lb_goal ) {
	options_->sap_lb_goal( lb_goal );
}

void
AddSapConstraintMover::set_penalty_per_sap( core::Real penalty ) {
	options_->penalty_per_sap( penalty );
}

void
AddSapConstraintMover::set_packing_correction( core::Real correction ) {
	options_->packing_correction( correction );
}

void
AddSapConstraintMover::set_full_accuracy_when_scoring( bool full_accuracy ) {
	options_->full_accuracy_when_scoring( full_accuracy );
}

void
AddSapConstraintMover::set_name( std::string const & name ) {
	options_->name( name );
}

std::string AddSapConstraintMover::get_name() const {
	return mover_name();
}

std::string AddSapConstraintMover::mover_name() {
	return "AddSapConstraintMover";
}

void AddSapConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		//cslist of variant types to add
		+ XMLSchemaAttribute::attribute_w_default(
		"sap_goal", xsct_real,
		"The final sap value you wish for your selection to achieve. There will be a penalty if the sap is above this number.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"penalty_per_sap", xsct_real,
		"If the current sap value is outside the sap_goal, what penalty shall be applied for each unit of sap?",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"sap_lb_goal", xsct_real,
		"The lowerbound on the final sap value you wish for your selection to achieve.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"packing_correction", xsct_real,
		"This method approximates sap during packing and it is often the case that it is wrong by a constant factor. This option allows"
		" you to adjust for this fact. After each packing trajectory, the SapConstraintHelper will output a message indicating"
		" what the algorithm thought the sap was, what it actually was, and what packing_correction you should use. Running a few trajectories"
		" should get you a good idea of what number you should place here. Setting it wrong will result in worse scoring poses.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"full_accuracy_when_scoring", xsct_rosetta_bool,
		"Debugging: Instead of using the full sap calculation when reporting scores, use instead the selected \"speed\" method.",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"speed", xs_string,
		"How fast should this calculation be during packing; faster is less precise. Choices: slow, fast, lightning. Slow uses full-atom"
		" calculations for sasa (but still uses approximate sasa). Fast and lightning assume that designable residues look like CYS for"
		" sasa calcualtions and only calculate sasa once at the start. Lightning approximates even further and treats each rotamer of an aa-type"
		" at the same position the same. Further, lightning only does the 0-clip on the whole-residue rather than at the per-atom level.",
		"fast" )
		+ XMLSchemaAttribute::attribute_w_default(
		"score_selector", xs_string,
		"Which residues should be included in the sap score? Optional, will default to full-pose.",
		"true_selector" )
		+ XMLSchemaAttribute(
		"sap_calculate_selector", xs_string,
		"Which residues should be present during the sap calculation? Only residues in the score_selector will have their values reported, but residues"
		" in this selector will still be assigned atom-saps which can affect the residues in score_selector. Optional, will default to score_selector." )
		+ XMLSchemaAttribute(
		"sasa_selector", xs_string,
		"Which residues should be present during the sasa calculation? Optional, will default to sap_calculate_selector." );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"A mover that adds the SapConstraint to the pose. The sap_constraint scoreterm then uses these constraints to try to match your"
		" sap_goal that you specified. A very important consideration is that all the speeds to choose from are approximations, so while"
		" the slower speeds are more precise, they are not necessarily more accurate. See \"packing_correction\" for more information."
		" See this paper for more info on sap: Developability index: a rapid in silico tool for the screening of antibody aggregation propensity. "
		" Lauer, et. al. J Pharm Sci 2012",
		attlist );
}

std::string AddSapConstraintMoverCreator::keyname() const {
	return AddSapConstraintMover::mover_name();
}

protocols::moves::MoverOP
AddSapConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddSapConstraintMover );
}

void AddSapConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddSapConstraintMover::provide_xml_schema( xsd );
}


} // moves
} // protocols

