// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/MultiplePoseMover.cc
/// @brief This mover accepts multiple poses from a previous mover,
///   performs selection using a provided pose selector,
///   applies contained ROSETTASCRIPTS protocol (ParsedProtocol),
///   and output multiple poses to the next mover of JD2.
/// @author Luki Goldschmidt (lugo@uw.edu)

// Unit headers
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>
#include <protocols/rosetta_scripts/MultiplePoseMoverCreator.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <protocols/rosetta_scripts/PoseSelectorFactory.hh>

// C/C++ headers
#include <iostream>
#include <list>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/options/option.hh>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.MultiplePoseMover" );

////////////////////////////////////////////////////////////////////////

// XRW TEMP std::string MultiplePoseMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "MultiplePoseMover";
// XRW TEMP }

// XRW TEMP std::string MultiplePoseMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP MultiplePoseMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new MultiplePoseMover() );
// XRW TEMP }

////////////////////////////////////////////////////////////////////////

MultiplePoseMover::MultiplePoseMover() :
	protocols::moves::Mover( "MultiplePoseMover" ),
	cached_(false),
	max_input_poses_(0),
	max_output_poses_(0),
	rosetta_scripts_tag_(/* NULL */),
	selector_tag_(/* NULL */),
	previous_mover_(/* NULL */),
	selectors_(),
	pose_input_cache_(),
	pose_output_cache_(),
	poses_input_(0),
	poses_output_(0),
	selector_filters_()
{
}

// XRW TEMP std::string MultiplePoseMover::get_name() const
// XRW TEMP {
// XRW TEMP  return MultiplePoseMover::mover_name();
// XRW TEMP }

/// @brief Process input pose
void MultiplePoseMover::apply(core::pose::Pose& pose)
{
	protocols::moves::Mover::set_last_move_status(protocols::moves::FAIL_RETRY);

	// Reset state
	pose_input_cache_.clear();
	pose_output_cache_.clear();
	poses_input_ = 0;
	poses_output_ = 0;

	// Clone provided pose since we may be returning a different based on Selector
	pose_input_cache_.push_back(pose.clone());
	++poses_input_;

	if ( cached_ ) {
		// Collect all remaining poses from previous mover
		while ( fill_input_cache() ) {}

		TR << "Collected input poses: " << pose_input_cache_.size() << std::endl;
		for ( core::pose::PoseOP p : pose_input_cache_ ) {
			TR << "\t" << p->sequence() << std::endl;
		}
	}

	core::pose::PoseOP output_pose( generate_pose() );
	if ( output_pose ) {
		pose = *output_pose;
		protocols::moves::Mover::set_last_move_status(protocols::moves::MS_SUCCESS);
	} else {
		TR << "No poses selected or processed succesfully by sub-mover. Setting status to FAIL_RETRY." << std::endl;
	}
}

/// @brief
core::pose::PoseOP MultiplePoseMover::get_additional_output()
{
	return generate_pose();
}

////////////////////////////////////////////////////////////////////////////

bool MultiplePoseMover::fill_input_cache()
{
	if ( !previous_mover_ ) {
		// No previous mover (source)
		return false;
	}

	if ( max_input_poses_ && poses_input_ >= max_input_poses_ ) {
		// Input limit reached
		return false;
	}

	TR << "Obtaining additional pose from previous mover: " << previous_mover_->get_name() << std::endl;
	core::pose::PoseOP another_pose = previous_mover_->get_additional_output();
	if ( !another_pose ) {
		// Previous mover out of poses
		return false;
	}

	pose_input_cache_.push_back(another_pose);
	++poses_input_;
	return true;
}

core::pose::PoseOP MultiplePoseMover::generate_pose()
{
	while ( !max_output_poses_ || poses_output_ < max_output_poses_ ) {

		// Grab a pose that we already processed from the pose_output_cache_
		if ( !pose_output_cache_.empty() ) {
			core::pose::PoseOP pose = pose_output_cache_.front();
			pose_output_cache_.pop_front();
			++poses_output_;
			if ( pose ) {
				return pose;
			}
		}

		// Get and process more poses, and put them into the pose_output_cache_

		// 1. Obtain input pose
		if ( !cached_ && pose_input_cache_.empty() ) {
			// Pull another pose from previous mover when not running in cached mode
			fill_input_cache(); // one iteration using get_additional_output()
		}

		if ( pose_input_cache_.empty() ) {
			// No more input poses
			return nullptr;
		}

		// 2. Select poses
		std::deque < core::pose::PoseOP > selected_poses = select_poses( pose_input_cache_ );

		// 3. Process selected poses and put them into pose_output_cache_
		if ( !selected_poses.empty() ) {
			std::deque < core::pose::PoseOP > poses = process_poses(selected_poses);
			for ( core::pose::PoseOP pose : poses ) {
				pose_output_cache_.push_back(pose);
			}
		}

		// We've processed these input poses, so we're done with them
		pose_input_cache_.clear();

	}

	return nullptr;
}

/// @brief Select poses from set using specified selectors
std::deque < core::pose::PoseOP > MultiplePoseMover::select_poses( std::deque < core::pose::PoseOP > & poses)
{
	// Process selection criteria
	if ( selectors_.empty() ) {
		// No selector specified -- select all poses
		return poses;
	}

	utility::vector1 < bool > selected_poses_by_selectors;
	std::deque < core::pose::PoseOP > selected_poses_for_processing;

	// Temp work around that shouldn't be a big performance hit
	utility::vector1 < core::pose::PoseOP > pose_vector;
	for ( core::pose::PoseOP p : poses ) {
		pose_vector.push_back(p);
	}

	// Apply selectors
	for ( PoseSelectorOP selector : selectors_ ) {
		selected_poses_by_selectors = selector->select_poses(pose_vector);
		// TODO: How do we handle multiple selectors? AND? OR?
		break;
	}

	if ( selectors_.size() > 1 ) {
		TR << selectors_.size() << " selectors defined; only the first selector is currently used!" << std::endl;
	}

	// Make a new vector of PoseOP's for selected poses for easier handling
	core::Size i = 1;
	for ( auto & pose : poses ) {
		if ( selected_poses_by_selectors[i] ) {
			selected_poses_for_processing.push_back(pose);
		}
		++i;
	}

	// Apply movers to selected poses
	TR << "Selected poses for processing: " << selected_poses_for_processing.size() << std::endl;
	for ( core::pose::PoseOP p : selected_poses_for_processing ) {
		TR << "\t" << p->sequence() << std::endl;
	}

	return selected_poses_for_processing;
}

/// @brief Rub sub-protocol on set of poses
std::deque < core::pose::PoseOP > MultiplePoseMover::process_poses( std::deque < core::pose::PoseOP > & poses )
{
	// Collect additional poses first rather than adding them to selected_poses_ right away
	// as it may invalidate the iterator pointer when memory location changes
	std::deque < core::pose::PoseOP > selected_poses;
	utility::vector1 < core::pose::PoseOP > additional_poses;

	for ( core::pose::PoseOP p : poses ) {
		TR << "Applying mover to pose: " << p->sequence() << std::endl;
		if ( process_pose(*p, additional_poses) ) {
			// Processing successful (from sub-protocol and its filters)
			selected_poses.push_back( p );
		}
	}

	if ( !additional_poses.empty() ) {
		for ( core::pose::PoseOP p : additional_poses ) {
			selected_poses.push_back( p );
		}
		TR << additional_poses.size() << " additional poses obtained; total output poses: " << selected_poses.size() << std::endl;
	}

	return selected_poses;
}

/// @brief Process a single input pose by the RosettaScripts mover
bool MultiplePoseMover::process_pose( core::pose::Pose & pose, utility::vector1 < core::pose::PoseOP > & additional_poses )
{
	if ( !rosetta_scripts_tag_ ) {
		return true;
	}

	protocols::rosetta_scripts::RosettaScriptsParser parser;

	// rosetta_scripts_tag_ has been pre-parsed in parse_my_tag() so no parsing exception should be thrown here
	// No longer true, in fact, as the original parse_my_tag call has been removed, since the Pose needed for that
	// call is not present until we get to the MPM's apply call.
	protocols::moves::MoverOP mover;
	try {
		rosetta_scripts_tag_->reset_accessed_options();
		mover = parser.parse_protocol_tag( pose, rosetta_scripts_tag_, basic::options::option );
		if ( !mover ) {
			TR << "Failed to parse protocol? This should not happen. Not applying protocol to pose." << std::endl;
			return false;
		}
	} catch ( utility::excn::Exception const & e ) {
		std::ostringstream oss;
		oss << "MultiplePoseMover could not create the inner parsed protocol; error message generated from parser.parse_protocol_tag:\n" << e.msg();
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}

	mover->apply(pose);

	if ( mover->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
		TR << "Sub-mover reported failure; passing status on." << std::endl;
		return false;
	}

	// Get additional poses from protocol mover;
	// these are selected for output by this mover
	core::pose::PoseOP p;
	while ( ( p = mover->get_additional_output() ) ) {
		additional_poses.push_back(p);
	}

	return true;
}

/// @brief Recurse through the tree of Tags and for each one, read each option, so that every option is set as "read"
///
/// @details As a safety check to make sure that no options are misspelled, at the conclusion of parsing, all options
/// contained in all tags are examined to determine if they have been accessed or not. If an option has not been
/// accessed, then an error message is written to the screen, and the program exits. However, it is not possible
/// to parse all movers in the absence of the Pose that is going to be used, and that Pose isn't available to the MPM
/// until execution. So we work around this issue by accessing each of the options in the input tag and all of its
/// children

void
recursively_access_all_attributes( utility::tag::TagCOP tag )
{
	auto const & options = tag->getOptions();
	for ( auto option_iter : options ) {
		tag->getOption< std::string >( option_iter.first );
	}
	for ( auto child : tag->getTags() ) {
		recursively_access_all_attributes( child );
	}
}

/// @brief Parse settings in tag
void MultiplePoseMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {

	if ( tag->hasOption("max_input_poses") ) {
		max_input_poses_ = tag->getOption<int>("max_input_poses", 0);
	}
	if ( tag->hasOption("max_output_poses") ) {
		max_output_poses_ = tag->getOption<int>("max_output_poses", 0);
	}
	if ( tag->hasOption("cached") ) {
		cached_ = tag->getOption<bool>("cached");
	}

	try {

		// ROSETTASCRIPTS tag (optional)
		if ( tag->hasTag("ROSETTASCRIPTS") ) {
			utility::tag::TagCOP rs_tag = tag->getTag("ROSETTASCRIPTS");
			set_rosetta_scripts_tag( rs_tag );
			// Try parsing the ROSETTASCRIPTS block to avoid tripping die_for_unaccessed_options()
			// protocols::rosetta_scripts::RosettaScriptsParser parser;
			// protocols::moves::MoverOP mover( parser.parse_protocol_tag( rosetta_scripts_tag_ ) );
			recursively_access_all_attributes( rs_tag );
		}

		// SELECT tag (optional)
		if ( tag->hasTag("SELECT") ) {
			selectors_.clear();
			TagCOP select_tag( tag->getTag("SELECT") );
			for ( TagCOP curr_tag : select_tag->getTags() ) {
				PoseSelectorOP new_selector(
					PoseSelectorFactory::get_instance()->
					newPoseSelector( curr_tag, data, selector_filters_, movers, pose )
				);
				runtime_assert( new_selector != nullptr );
				selectors_.push_back( new_selector );
				TR << "Defined selector of type " << curr_tag->getName() << std::endl;
			}
		}

		// Warn if no ROSETTASCRIPTS protocol and no SELECTOR (i.e. null mover)
		if ( selectors_.size() < 1 && !rosetta_scripts_tag_ ) {
			std::string my_name( tag->getOption<std::string>("name") );
			TR.Warning << "Neither a ROSETTASCRIPTS protocol nor a SELECT statement specified in MultiplePoseMover with name \"" << my_name << "\". This mover has no effect. Are you sure this is what you intended?" << std::endl;
		}

		// TODO: Should we complain here is there are tags specified that we don't understand?

	} catch( utility::excn::Exception const & e ) {
		std::string my_name( tag->getOption<std::string>("name") );
		throw CREATE_EXCEPTION(utility::excn::Exception, "Exception in MultiplePoseMover with name \"" + my_name + "\": " + e.msg());
	}

	// Obtain flags from selector
	PoseSelectorFlags flags = PSF_NONE;
	for ( PoseSelectorOP selector : selectors_ ) {
		// flags |= selector->get_flags();
		flags = (PoseSelectorFlags)( flags | selector->get_flags() );
	}

	if ( (flags & PSF_NEED_FULL_POSE_SET) && !cached_ ) {
		cached_ = true;
		if ( tag->hasOption("cached") ) {
			TR.Warning << "Overriding \"cached\" option due to requirements in specified Selectors." << std::endl;
		}
	}

	if ( max_input_poses_ ) {
		TR << "Max input poses: " << max_input_poses_ << std::endl;
	}
	if ( max_output_poses_ ) {
		TR << "Max output poses: " << max_output_poses_ << std::endl;
	}
	if ( cached_ ) {
		TR << "Pose input caching enabled" << std::endl;
	}
}

/// @brief Used by RosettaScripts to set the previous mover to pull poses from
void MultiplePoseMover::set_previous_mover( protocols::moves::MoverOP const m ) { previous_mover_ = m; }

/// @brief sets rosettascripts tag
void MultiplePoseMover::set_rosetta_scripts_tag( utility::tag::TagCOP tag ) { rosetta_scripts_tag_ = tag; }

std::string MultiplePoseMover::get_name() const {
	return mover_name();
}

std::string MultiplePoseMover::mover_name() {
	return "MultiplePoseMover";
}

std::string mpm_mangler( std::string const & element_name )
{ return "multiple_pose_mover_" + element_name + "_type"; }

void MultiplePoseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	RosettaScriptsParser::write_ROSETTASCRIPTS_complex_type( xsd );

	using Attr = XMLSchemaAttribute;
	AttributeList attlist;
	attlist
		+ Attr( "max_input_poses", xsct_non_negative_integer, "XSD TO DO" )
		+ Attr( "max_output_poses", xsct_non_negative_integer, "XSD TO DO" )
		+ Attr( "cached", xsct_rosetta_bool, "XSD TO DO" );

	PoseSelectorFactory::get_instance()->define_pose_selector_group( xsd );

	XMLSchemaSimpleSubelementList pose_selector_subelement;
	pose_selector_subelement.add_group_subelement( & PoseSelectorFactory::pose_selector_group_name );
	XMLSchemaComplexTypeGenerator select_ct;
	select_ct.element_name( "SELECT" )
		.complex_type_naming_func( & mpm_mangler )
		.set_subelements_repeatable( pose_selector_subelement )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList rs_element, select_element;
	rs_element.add_already_defined_subelement(
		RosettaScriptsParser::rosetta_scripts_element_name(),
		& RosettaScriptsParser::rosetta_scripts_complex_type_naming_func );
	select_element.add_already_defined_subelement(
		"SELECT", & mpm_mangler );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( mover_name() )
		.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.add_ordered_subelement_set_as_optional( select_element )
		.add_ordered_subelement_set_as_optional( rs_element )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );
}

std::string MultiplePoseMoverCreator::keyname() const {
	return MultiplePoseMover::mover_name();
}

protocols::moves::MoverOP
MultiplePoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MultiplePoseMover );
}

void MultiplePoseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MultiplePoseMover::provide_xml_schema( xsd );
}

} //rosetta_scripts
} //protocols
