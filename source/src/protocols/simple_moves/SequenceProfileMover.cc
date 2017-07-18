// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SequenceProfileMover.cc
/// @brief  moverization of code that was embedded in the parser
/// @author Brian Weitzner brian.weitzner@gmail.com, Steven Lewis smlewi@gmail.com
/// @date   Rebased to next year.

// Unit Headers
#include <protocols/simple_moves/SequenceProfileMover.hh>
#include <protocols/simple_moves/SequenceProfileMoverCreator.hh>


// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/BasicFilters.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jd2/util.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/sequence/SequenceProfile.hh>
#include <utility/file/FileName.hh>

#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SequenceProfileMover" );

namespace protocols {
namespace simple_moves {

SequenceProfileMover::SequenceProfileMover() :
	protocols::moves::Mover( "SequenceProfileMover" )
{
}

SequenceProfileMover::~SequenceProfileMover() = default;

void
SequenceProfileMover::apply( core::pose::Pose & pose )
{
	using namespace core::sequence;

	SequenceProfileOP profile( new SequenceProfile );
	profile->read_from_checkpoint( cst_file_name_ );
	for ( core::Size seqpos( 1 ), end( pose.size() ); seqpos <= end; ++seqpos ) {
		pose.add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( pose, seqpos, profile ) ) ) );
	}

	TR << "Added sequence profile constraints specified in file " << cst_file_name_ << "." << std::endl;
}

// XRW TEMP std::string
// XRW TEMP SequenceProfileMover::get_name() const {
// XRW TEMP  return SequenceProfileMover::mover_name();
// XRW TEMP }

protocols::moves::MoverOP
SequenceProfileMover::clone() const{
	return protocols::moves::MoverOP( new SequenceProfileMover( *this ) );
}

protocols::moves::MoverOP
SequenceProfileMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SequenceProfileMover );
}

void
SequenceProfileMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	//handle cst_file_name
	std::string const input_file_name( protocols::jd2::current_input_tag() );
	core::Size const wheres_period( input_file_name.find_first_of( "." ) );
	std::string const dflt_cst_file_name( input_file_name.substr(0, wheres_period ) + ".cst" );
	set_cst_file_name( tag->getOption< std::string >( "file_name", dflt_cst_file_name ) );

	//handle profile_wgt
	set_profile_wgt( tag->getOption< core::Real >( "weight", 0.25 ) );

	if ( profile_wgt_ ) {
		using namespace utility::pointer;
		for ( std::map< std::string, ReferenceCountOP >::const_iterator it=data[ "scorefxns" ].begin(); it!=data[ "scorefxns" ].end(); ++it ) {
			ScoreFunctionOP scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", it->first ) );
			scorefxn->set_weight( res_type_constraint, profile_wgt_ );
			TR << "setting " << it->first << " res_type_constraint to " << profile_wgt_ << "\n";
		}
	}
	TR << "Changed all scorefxns to have profile weights of " << profile_wgt_ << std::endl;
} //end parse_my_tag

// XRW TEMP std::string
// XRW TEMP SequenceProfileMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SequenceProfileMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SequenceProfileMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SequenceProfileMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SequenceProfileMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "profile";
// XRW TEMP }

std::string SequenceProfileMover::get_name() const {
	return mover_name();
}

std::string SequenceProfileMover::mover_name() {
	return "profile";
}

void SequenceProfileMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "file_name", xs_string, "XRW_TODO" )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "XRW TO DO", "0.25" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW_TODO", attlist );
}

std::string SequenceProfileMoverCreator::keyname() const {
	return SequenceProfileMover::mover_name();
}

protocols::moves::MoverOP
SequenceProfileMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SequenceProfileMover );
}

void SequenceProfileMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SequenceProfileMover::provide_xml_schema( xsd );
}



} // simple_moves
} // protocols

