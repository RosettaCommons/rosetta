// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ContingentAcceptMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/ContingentAcceptMover.hh>
#include <protocols/simple_moves/ContingentAcceptMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.simple_moves.ContingentAcceptMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP ContingentAcceptMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ContingentAcceptMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ContingentAcceptMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ContingentAcceptMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ContingentAcceptMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ContingentAccept";
// XRW TEMP }

ContingentAcceptMover::ContingentAcceptMover()
: moves::Mover("ContingentAccept"),
	filter_(/* NULL */),
	mover_(/* NULL */),
	delta_(5)
{
}

protocols::filters::FilterOP
ContingentAcceptMover::filter() const{
	return filter_;
}

void
ContingentAcceptMover::filter( protocols::filters::FilterOP f ){
	filter_ = f;
}

void
ContingentAcceptMover::apply( Pose & pose )
{
	core::pose::Pose const old_pose = pose;
	core::Real preMoveScore;
	preMoveScore = filter()->report_sm( pose );
	mover()->apply(pose);
	core::Real postMoveScore;
	postMoveScore = filter()->report_sm( pose );

	if ( postMoveScore <= preMoveScore + delta() ) {
		TR<<"postMoveScore is "<<postMoveScore<<" preMoveScore is "<<preMoveScore<<" . Accepting new pose."<<std::endl;
		return;
	} else {
		TR<<"postMoveScore is "<<postMoveScore<<" preMoveScore is "<<preMoveScore<<" . Rejecting new pose."<<std::endl;
		pose = old_pose;
	}
}

// XRW TEMP std::string
// XRW TEMP ContingentAcceptMover::get_name() const {
// XRW TEMP  return ContingentAcceptMover::mover_name();
// XRW TEMP }

moves::MoverOP
ContingentAcceptMover::clone() const
{
	return moves::MoverOP( new ContingentAcceptMover( *this ) );
}

moves::MoverOP
ContingentAcceptMover::fresh_instance() const
{
	return moves::MoverOP( new ContingentAcceptMover );
}

/*
returns the contains mover
*/
protocols::moves::MoverOP ContingentAcceptMover::mover() const {
	return mover_;
}

/*
Setting the internal mover to point to m.
*/
void ContingentAcceptMover::mover( protocols::moves::MoverOP m ){
	mover_ = m;
}

void
ContingentAcceptMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & )
{
	filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "mover" ), movers ) );
	delta( tag->getOption< core::Real >( "delta" , 5.0 )); // do I need anything like ( "MaxDeltaFilterVal", 5 )?? What is the integer for?
}

std::string ContingentAcceptMover::get_name() const {
	return mover_name();
}

std::string ContingentAcceptMover::mover_name() {
	return "ContingentAccept";
}

void ContingentAcceptMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "filter", xs_string, "Filter whose value is compared before and after applying the Mover; looked up from DataMap.")
		+ XMLSchemaAttribute::required_attribute( "mover", xs_string, "Mover whose move is contingently accepted.")
		+ XMLSchemaAttribute::attribute_w_default( "delta", xsct_real, "Amount Filter value can increase before rejection (units depend on Filter)", "5.0");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "checks filter, runs mover, checks filter.  Resets pose if filter value increases by more than delta. NOT Monte Carlo.", attlist );
}

std::string ContingentAcceptMoverCreator::keyname() const {
	return ContingentAcceptMover::mover_name();
}

protocols::moves::MoverOP
ContingentAcceptMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ContingentAcceptMover );
}

void ContingentAcceptMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ContingentAcceptMover::provide_xml_schema( xsd );
}




} // simple_moves
} // protocols
