// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/MultipleOutputWrapper.cc
/// @brief  This mover wraps another mover or a ROSETTASCRIPTS block to obtain additional,
//    related poses. A new instance of the mover or protocol is created for each iteration.
/// @author Luki Goldschmidt (lugo@uw.edu)

// Unit headers
#include <protocols/rosetta_scripts/MultipleOutputWrapper.hh>
#include <protocols/rosetta_scripts/MultipleOutputWrapperCreator.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Boost headers

// C++ headers
#include <iomanip>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.MultipleOutputWrapper" );

using namespace protocols::moves;

////////////////////////////////////////////////////////////////////////

// XRW TEMP std::string MultipleOutputWrapper::mover_name()
// XRW TEMP {
// XRW TEMP  return "MultipleOutputWrapper";
// XRW TEMP }

// XRW TEMP std::string MultipleOutputWrapperCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return mover_name();
// XRW TEMP }

// XRW TEMP MoverOP MultipleOutputWrapperCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return MoverOP( new MultipleOutputWrapper() );
// XRW TEMP }

////////////////////////////////////////////////////////////////////////

MultipleOutputWrapper::MultipleOutputWrapper() :
	Mover( "MultipleOutputWrapper" ),
	mover_tag_(/* NULL */),
	rosetta_scripts_tag_(/* NULL */),
	reference_pose_(/* NULL */),
	max_poses_(0),
	max_attempts_(10),
	n_poses_(0),
	keep_mover_state_(false),
	mover_()
{
}

// XRW TEMP std::string MultipleOutputWrapper::get_name() const
// XRW TEMP {
// XRW TEMP  return MultipleOutputWrapper::mover_name();
// XRW TEMP }

/// @brief Process all input poses (provided pose and from previous mover)
void MultipleOutputWrapper::apply(core::pose::Pose& pose)
{
	reference_pose_ = core::pose::PoseOP( new core::pose::Pose(pose) );
	generate_pose(pose);
}

/// @brief Hook for multiple pose putput to JD2 or another mover
core::pose::PoseOP MultipleOutputWrapper::get_additional_output()
{
	if ( !reference_pose_ ) {
		return nullptr;
	}
	if ( (max_poses_ > 0) && (n_poses_ >= max_poses_) ) {
		return nullptr;
	}

	core::pose::PoseOP new_pose( new core::pose::Pose(*reference_pose_) );
	if ( !generate_pose(*new_pose) ) {
		return nullptr;
	}

	return new_pose;
}

bool MultipleOutputWrapper::generate_pose(core::pose::Pose & pose)
{
	// Empty objects... may not work...
	basic::datacache::DataMap data;
	protocols::filters::Filters_map filters;
	protocols::moves::Movers_map movers;
	if ( !keep_mover_state_ ) {
		mover_ = nullptr;
	}

	if ( !mover_ && rosetta_scripts_tag_ ) {
		protocols::rosetta_scripts::RosettaScriptsParser parser;
		mover_ = parser.parse_protocol_tag( pose, rosetta_scripts_tag_, basic::options::option );
	}

	if ( !mover_ && mover_tag_ ) {
		mover_ = MoverFactory::get_instance()->newMover(mover_tag_, data, filters, movers, pose );
	}

	runtime_assert( mover_ != nullptr );

	core::Size attempts;
	for ( attempts = 1; attempts <= max_attempts_; ++attempts ) {

		mover_->apply(pose);

		protocols::moves::MoverStatus status = mover_->get_last_move_status();
		set_last_move_status(status);
		if ( status != protocols::moves::MS_SUCCESS ) {
			TR << "Sub-mover or protocol reported failure on attempt " << attempts << " of " << max_attempts_ << std::endl;
			if ( reference_pose_ ) {
				pose = *reference_pose_;
			}
			continue;
		}

		if ( ! pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
			// Add new output tag
			using basic::datacache::DataCache_CacheableData;
			std::ostringstream tag;
			tag << name_ << "_" << std::setw(4) << std::setfill('0') << (n_poses_+1);
			pose.data().set(
				core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
				DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( tag.str() ) )
			);
		}

		++n_poses_;
		return true;
	}

	if ( !attempts ) {
		set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
	}

	TR << "Could not generate a pose" << std::endl;
	return false;
}

void MultipleOutputWrapper::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {

	if ( tag->hasOption("name") ) {
		name_ = tag->getOption<std::string>("name");
	}

	if ( tag->hasOption("max_output_poses") ) {
		max_poses_ = tag->getOption<int>("max_output_poses", 0);
	}

	if ( tag->hasOption("max_attempts") ) {
		max_attempts_ = tag->getOption<int>("max_attempts", 0);
	}

	if ( tag->hasOption("keep_mover_state") ) {
		keep_mover_state_ = tag->getOption<bool>("keep_mover_state");
	}

	try {

		// Children of tag are movers
		for ( utility::tag::TagCOP const curr_tag : tag->getTags() ) {
			// Try instantiating first mover from tag to test parsing
			if ( curr_tag->getName() == "ROSETTASCRIPTS" ) {
				// Treat subtag as a ROSETTASCRIPTS protocol
				protocols::rosetta_scripts::RosettaScriptsParser parser;
				protocols::moves::MoverOP mover( parser.parse_protocol_tag( curr_tag, basic::options::option ) );
				rosetta_scripts_tag_ = curr_tag;
			} else {
				// Treat subtag as a regular mover tag
				std::string name = curr_tag->getOption<std::string>("name");
				protocols::moves::MoverOP new_mover(
					protocols::moves::MoverFactory::get_instance()->
					newMover(curr_tag, data, filters, movers, pose)
				);
				mover_tag_ = curr_tag;
			}
			// Only first mover used -- add warning when multiple defined?
			break;
		}

		if ( !mover_tag_ && !rosetta_scripts_tag_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception, "No mover or ROSETTASCRIPTS tag found.");
		}

	} catch( utility::excn::Exception const & e ) {
		std::string my_name( tag->getOption<std::string>("name") );
		throw CREATE_EXCEPTION(utility::excn::Exception, "Exception in MultipleOutputWrapper with name \"" + my_name + "\": " + e.msg());
	}
}

std::string MultipleOutputWrapper::get_name() const {
	return mover_name();
}

std::string MultipleOutputWrapper::mover_name() {
	return "MultipleOutputWrapper";
}

void MultipleOutputWrapper::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	RosettaScriptsParser::write_ROSETTASCRIPTS_complex_type( xsd );

	typedef XMLSchemaAttribute Attr;
	AttributeList attlist;
	attlist
		+ optional_name_attribute()
		+ Attr( "max_output_poses", xsct_non_negative_integer, "XRW TO DO" )
		+ Attr( "max_attempts", xsct_non_negative_integer, "XRW TO DO" )
		+ Attr( "keep_mover_state", xsct_rosetta_bool, "XRW TO DO" );

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement(
		RosettaScriptsParser::rosetta_scripts_element_name(),
		& RosettaScriptsParser::rosetta_scripts_complex_type_naming_func );
	subelements.add_group_subelement( & moves::MoverFactory::mover_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.set_subelements_pick_one( subelements )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );

}

std::string MultipleOutputWrapperCreator::keyname() const {
	return MultipleOutputWrapper::mover_name();
}

protocols::moves::MoverOP
MultipleOutputWrapperCreator::create_mover() const {
	return protocols::moves::MoverOP( new MultipleOutputWrapper );
}

void MultipleOutputWrapperCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MultipleOutputWrapper::provide_xml_schema( xsd );
}


} //rosetta_scripts
} //protocols
