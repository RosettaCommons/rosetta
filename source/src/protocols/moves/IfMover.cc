// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/IfMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Christopher Miles (cmiles@uw.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) Added logic.

// Unit headers
#include <protocols/moves/IfMover.hh>
#include <protocols/moves/IfMoverCreator.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/pointer/memory.hh>

namespace protocols {
namespace moves {

static basic::Tracer TR( "protocols.moves.IfMover" );

void IfMover::apply(core::pose::Pose& pose) {
	moves::MoverStatus status;


	if ( filter_->apply(pose) ) {
		true_mover_->apply(pose);
		status = true_mover_->get_last_move_status();
	} else {
		false_mover_->apply(pose);
		status = false_mover_->get_last_move_status();
	}

	// update this mover's status
	protocols::moves::Mover::set_last_move_status(status);
}

core::pose::PoseOP IfMover::get_additional_output_true_mover() {
	return true_mover_->get_additional_output();
}

core::pose::PoseOP IfMover::get_additional_output_false_mover() {
	return false_mover_->get_additional_output();
}

// backwards compatibility
core::pose::PoseOP IfMover::get_additional_output() {
	return get_additional_output_true_mover();
}


void IfMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	using namespace protocols::filters;
	//TR<<"If mover" << std::endl;
	true_mover_ = find_mover_or_die("null", tag, data);
	false_mover_ = find_mover_or_die("null", tag, data);
	if ( tag->hasOption("logic") ) {

		std::string logic = tag->getOption< std::string >("logic");
		//TR << "Parsing logic " << logic << std::endl;
		utility::vector1<std::string> logicSP = utility::string_split_simple(logic);
		if ( logicSP.size() < 4 ) {
			utility_exit_with_message("Logic string must be at least - if x : y");
		}
		core::Size filter_token = 2;
		core::Size true_token = 4;
		core::Size false_token = 6;
		bool flip=false;
		if ( utility::uppercased(logicSP[2]) == "NOT" ) {
			filter_token+=1;
			true_token+=1;
			false_token+=1;
			flip=true;
			//TR << "Flipping" <<std::endl;
		}
		filter_ = find_filter_or_die(logicSP[filter_token], tag, data);

		if ( flip ) {
			//TR << "Setting False " << logicSP[true_token] << std::endl;
			false_mover_ = find_mover_or_die(logicSP[true_token], tag, data);
		} else {
			//TR << "Setting True " << logicSP[true_token] << std::endl;
			true_mover_ = find_mover_or_die(logicSP[true_token], tag, data);
		}

		if ( logicSP.size() >= 6 ) {
			//TR << "False" << std::endl;
			if ( flip ) {
				true_mover_ = find_mover_or_die(logicSP[false_token], tag, data);
			} else {
				false_mover_ = find_mover_or_die(logicSP[false_token], tag, data);
			}
		}
		//TR <<"Parsed logic" << std::endl;
	} else {

		std::string filter_name = "";
		std::string const true_mover_name( tag->getOption< std::string >( "true_mover_name", "null" ));
		std::string const false_mover_name( tag->getOption< std::string >( "false_mover_name", "null" ));


		if ( true_mover_name == "null" && false_mover_name=="null" ) {
			utility_exit_with_message("Both movers cannot be null!");
		}

		/// see: protocols/moves/util.hh

		true_mover_ = find_mover_or_die(true_mover_name, tag, data);
		false_mover_ = find_mover_or_die(false_mover_name, tag, data);

		//Parse the filter either by value or identity.
		if ( tag->hasOption("value") ) {
			if ( tag->getOption< bool >("value") == true ) {
				filter_name= "TrueFilter";
				filter_ = utility::pointer::make_shared< filters::TrueFilter >();
			} else {
				filter_name="FalseFilter";
				filter_ = utility::pointer::make_shared< filters::FalseFilter >();
			}
		} else {
			filter_name = tag->getOption< std::string >( "filter_name" ) ;
			filter_ = find_filter_or_die(filter_name, tag, data);
		}

		TR << "with true_mover \"" << true_mover_name
			<< "\" and false_mover \"" << false_mover_name
			<< "\" filter \"" << filter_name
			<< std::endl;
	}
}

std::string IfMover::get_name() const {
	return mover_name();
}

std::string IfMover::mover_name() {
	return "If";
}

void IfMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"logic", xs_string,
		"Parse logic as `if x : y else z` This corresponds to filter_name,true_mover_name,false_mover_name. Use null for do nothing. not is also accepted so - if not x : y else z");
	attlist + XMLSchemaAttribute(
		"value",
		xs_boolean,
		"Alternative to filter_name/logic.  Allows control-flow through RS.  Set as 1 or 0");

	attlist + XMLSchemaAttribute(
		"filter_name", xs_string,
		"Filter used for the if else statement");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"true_mover_name", xs_string,
		"Mover to be execuated when filter returns true",
		"null");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"false_mover_name", xs_string,
		"Mover to be executed when filter returns false",
		"null");
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Implements a simple "
		"IF (filter(pose)) THEN true_mover(pose) ELSE false_mover(pose). "
		"If not using the logic option, true_mover is required, false_mover is not. Both movers default to null.",
		attlist );
}

std::string IfMoverCreator::keyname() const {
	return IfMover::mover_name();
}

protocols::moves::MoverOP
IfMoverCreator::create_mover() const {
	return utility::pointer::make_shared< IfMover >();
}

void IfMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IfMover::provide_xml_schema( xsd );
}


} //moves
} //protocols
