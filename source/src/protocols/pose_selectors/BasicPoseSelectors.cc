// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/pose_selectors/BasicPoseSelectors.hh
/// @brief  Collection of simple pose selectors
/// @author Luki Goldschmidt <lugo@uw.edu>


#ifndef INCLUDED_protocols_pose_selectors_BasicPoseSelectors_cc
#define INCLUDED_protocols_pose_selectors_BasicPoseSelectors_cc

// Unit Headers
#include <protocols/pose_selectors/BasicPoseSelectors.hh>
#include <protocols/pose_selectors/BasicPoseSelectorCreators.hh>
#include <protocols/rosetta_scripts/PoseSelectorFactory.hh>
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/filters/FilterFactory.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sort_predicates.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

static basic::Tracer TR( "protocols.pose_selectors.BasicPoseSelectors" );

namespace protocols {
namespace pose_selectors {

////////////////////////////////////////////////////////////////////////
// LogicalSelector base class

// Selector
LogicalSelector::LogicalSelector()
{
}


LogicalSelector::~LogicalSelector() = default;


utility::tag::XMLSchemaComplexTypeGeneratorOP
LogicalSelector::complex_type_generator_for_logical_selector( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	rosetta_scripts::PoseSelectorFactory::get_instance()->define_pose_selector_group( xsd );
	//AttributeList attlist;
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & rosetta_scripts::PoseSelectorFactory::pose_selector_group_name );
	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen
		->set_subelements_repeatable( subelements )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PoseSelectorFactory::complex_type_name_for_pose_selector );

	return ct_gen;
	//The children of this tag are other pose selectors

}

protocols::rosetta_scripts::PoseSelectorFlags LogicalSelector::get_flags() const
{
	protocols::rosetta_scripts::PoseSelectorFlags flags(protocols::rosetta_scripts::PSF_NONE);
	for ( protocols::rosetta_scripts::PoseSelectorOP selector : selectors_ ) {
		flags = (protocols::rosetta_scripts::PoseSelectorFlags)( flags | selector->get_flags() );
	}
	return flags;
}

void LogicalSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	// Children of tag are selectors
	for ( utility::tag::TagCOP const curr_tag : tag->getTags() ) {
		protocols::rosetta_scripts::PoseSelectorOP new_selector(
			protocols::rosetta_scripts::PoseSelectorFactory::get_instance()->
			newPoseSelector( curr_tag, data, filters, movers, pose )
		);
		runtime_assert( new_selector != nullptr );
		selectors_.push_back( new_selector );
		TR << "Defined pose selector of type " << curr_tag->getName() << std::endl;
	}
}

utility::vector1<bool> LogicalSelector::select_poses(
	utility::vector1< core::pose::PoseOP > poses
)
{
	utility::vector1<bool> selected_poses;

	TR << "Applying selector " << get_name() << std::endl;

	selected_poses.resize( poses.size(), get_default() );

	for ( protocols::rosetta_scripts::PoseSelectorOP selector : selectors_ ) {
		utility::vector1<bool> selector_selected_poses( selector->select_poses( poses ) );

		// Merge sets using logical operator
		for (
				auto
				i = selected_poses.begin(),
				j = selector_selected_poses.begin();
				i != selected_poses.end() &&
				j != selector_selected_poses.end();
				++i, ++j ) {

			(*i) = selection_operation(*i, *j);
		}

		TR.Debug << "Pose selections for " << get_name() << " after " << selector->get_name() << ": ";
		for ( bool const selection : selected_poses ) {
			TR.Debug << selection << " ";
		}
		TR.Debug << std::endl;
	}

	TR.Debug << "Final pose selections for " << get_name() << ": ";
	for ( bool const selection : selected_poses ) {
		TR.Debug << selection << " ";
	}
	TR.Debug << std::endl;

	return selected_poses;
}


////////////////////////////////////////////////////////////////////////
// AndSelector

// Creator
protocols::rosetta_scripts::PoseSelectorOP AndSelectorCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new AndSelector() );
}

void
AndSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AndSelector::provide_xml_schema( xsd );
}

// Selector
void
AndSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen = LogicalSelector::complex_type_generator_for_logical_selector( xsd );
	ct_gen->element_name( name() )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );
}


////////////////////////////////////////////////////////////////////////
// OrSelector

// Creator
protocols::rosetta_scripts::PoseSelectorOP OrSelectorCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new OrSelector() );
}

void
OrSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	OrSelector::provide_xml_schema( xsd );
}


void
OrSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){

	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen = LogicalSelector::complex_type_generator_for_logical_selector( xsd );
	ct_gen->element_name( name() )
		.description( "XRW TO DO" )
		.write_complex_type_to_schema( xsd );
}
////////////////////////////////////////////////////////////////////////
// TopNByProperty

// Creator
protocols::rosetta_scripts::PoseSelectorOP TopNByPropertyCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new TopNByProperty() );
}


void
TopNByPropertyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	TopNByProperty::provide_xml_schema( xsd );
}

// Selector
TopNByProperty::TopNByProperty() :
	reporter_(/* NULL */),
	order_(1)
{
}
void
TopNByProperty::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;


	rosetta_scripts::PosePropertyReporterFactory::get_instance()->define_pose_reporter_group( xsd );

	XMLSchemaRestriction order_string;
	order_string.name( "order_string");
	order_string.base_type( xs_string );
	order_string.add_restriction( xsr_enumeration, "asc" );
	order_string.add_restriction( xsr_enumeration, "ascending" );
	order_string.add_restriction( xsr_enumeration, "desc" );
	order_string.add_restriction( xsr_enumeration, "descending" );
	xsd.add_top_level_element( order_string );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "n", xsct_non_negative_integer, "Selection limit" )
		+ XMLSchemaAttribute( "order", "order_string", "Ascending or descending order?" );
	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_group_subelement( & rosetta_scripts::PosePropertyReporterFactory::pose_reporter_group_name );


	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen
		.element_name( name() )
		.set_subelements_repeatable( subelements, 1, 1 )
		.add_attributes( attlist )
		.description( "XRW TO DO" )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PoseSelectorFactory::complex_type_name_for_pose_selector )
		.write_complex_type_to_schema( xsd );

}

void TopNByProperty::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	// n = selection limit
	n_ = tag->getOption<core::Size>("n");

	if ( tag->hasOption("order") ) {
		std::string order( tag->getOption<std::string>("order") );
		if ( order == "asc" || order == "ascending" ) {
			order_ = 1;
		} else if ( order == "desc" || order == "descending" ) {
			order_ = -1;
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("Unknown order: " + order);
		}
	}

	// Children of tag are reporters
	for ( utility::tag::TagCOP const curr_tag : tag->getTags() ) {
		protocols::rosetta_scripts::PosePropertyReporterOP new_reporter(
			protocols::rosetta_scripts::PosePropertyReporterFactory::get_instance()->
			newPosePropertyReporter( curr_tag, data, filters, movers, pose )
		);
		runtime_assert( new_reporter != nullptr );
		reporter_ = new_reporter;
		TR << "Defined pose property reporter of type " << curr_tag->getName() << std::endl;
		// Only first reporter used -- add warning when multiple defined?
		break;
	}
}

utility::vector1<bool> TopNByProperty::select_poses(
	utility::vector1< core::pose::PoseOP > poses
)
{
	utility::vector1<bool> selected_poses;

	typedef std::pair< core::Size, core::Real > Pose_Property;
	utility::vector1 < Pose_Property > pose_properties;

	TR << "Applying selector " << get_name() << std::endl;

	// Obtain properties of all poses
	{
		core::Size i = 1;
		for ( core::pose::PoseOP pose : poses ) {
			core::Real r = reporter_->report_property( *pose );
			pose_properties.push_back( Pose_Property(i, r) );
			++i;
		}
	}

	// Sort properties vector by reported property (second)
	// order_ =  1: ascending (default)
	// order_ = -1: descending
	// order_ =  0: no sorting
	if ( order_ != 0 ) {
		std::sort(pose_properties.begin(), pose_properties.end(), utility::SortSecond< core::Size, core::Real >());
		if ( order_ < 1 ) {
			std::reverse(pose_properties.begin(), pose_properties.end());
		}
	}

	// Debug
	TR.Debug << "Sorted poses:" << std::endl;
	for ( auto & p : pose_properties ) {
		TR.Debug << p.first << " = " << p.second << std::endl;
	}

	// Create selected poses vector
	selected_poses.resize(poses.size(), false);
	core::Size n = poses.size() < n_ ? poses.size() : n_;
	for ( core::Size i = 1; i <= n; ++i ) {
		selected_poses[ pose_properties[i].first ] = true;
	}

	return selected_poses;
}

////////////////////////////////////////////////////////////////////////
// Filter

// Creator
protocols::rosetta_scripts::PoseSelectorOP FilterCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new Filter() );
}


void
FilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	Filter::provide_xml_schema( xsd );
}

// Selector
Filter::Filter() :
	filter_(/* NULL */)
{
}

void
Filter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	protocols::filters::FilterFactory::get_instance()->define_filter_xml_schema( xsd );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "filter", xs_string, "Name attribute of a filter that was defined earlier in the RosettaScript.");
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & filters::FilterFactory::filter_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen
		.element_name( name() )
		.set_subelements_repeatable( subelements, 0, 1 )
		.add_attributes( attlist )
		.description( "XRW TO DO" )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PoseSelectorFactory::complex_type_name_for_pose_selector )
		.write_complex_type_to_schema( xsd );

}

void Filter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	using namespace utility::tag;
	using namespace protocols::rosetta_scripts;

	TagCOP filter_tag(nullptr);
	std::string filter_name;

	if ( tag->hasOption("filter") ) {
		// Find a filter by name defined somewhere upstream in the script
		RosettaScriptsParser parser;
		filter_name = tag->getOption<std::string>("filter");
		filter_tag = parser.find_rosettascript_tag(
			tag,
			"FILTERS",
			"name",
			filter_name
		);

		if ( !filter_tag ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Cannot find filter named \"" + filter_name + "\"");
		}

	} else {
		// Filter is defined inline (first child tag)
		utility::vector0< TagCOP > tags( tag->getTags() );
		for ( utility::vector0< TagCOP >::const_iterator it = tags.begin(); it != tags.end(); ++it ) {
			filter_tag = *it;
			break;
		}
	}

	if ( filter_tag ) {
		if ( filter_name.empty() && filter_tag->hasOption("name") ) {
			filter_name = filter_tag->getOption<std::string>("name");
		}
		filter_  = protocols::filters::FilterFactory::get_instance()->newFilter( filter_tag, data, filters, movers, pose );
	}

	if ( !filter_ ) {
		std::ostringstream s;
		s << "Cannot create filter from script tag: " << tag;
		throw utility::excn::EXCN_RosettaScriptsOption(s.str());
	}

	if ( filters.find(filter_name) != filters.end() ) {
		TR.Warning << "Filter named \"" << filter_name << "\" already defined. Not adding this filter instance to the map." << std::endl;
	} else {
		filters.insert( std::make_pair( filter_name, filter_ ) );
		TR << "Defined filter named \"" << filter_name << "\" of type " << filter_tag->getName() << std::endl;
	}
}

utility::vector1<bool> Filter::select_poses(
	utility::vector1< core::pose::PoseOP > poses
)
{
	utility::vector1<bool> selected_poses(poses.size(), true);

	if ( !filter_ ) {
		TR << "No filter instance!" << std::endl;
		return selected_poses;
	}

	TR << "Applying selector " << get_name() << ": " << filter_->get_user_defined_name() << std::endl;

	core::Size i = 1;
	for ( core::pose::PoseOP pose : poses ) {
		TR.Debug << "Pose " << i << "..." << std::endl;

		bool ok = filter_->apply( *pose );
		core::Real const filter_value( filter_->report_sm( *pose ) );
		setPoseExtraScore( *pose, filter_->get_user_defined_name(), (float)filter_value );

		TR.Debug << "Pose " << i << ": " << (ok ? "Pass" : "Fail") << std::endl;
		selected_poses[i] = ok;
		++i;
	}

	return selected_poses;
}

////////////////////////////////////////////////////////////////////////

} // pose_selectors
} // protocols

#endif //INCLUDED_protocols_pose_selectors_BasicPoseSelectors_cc
