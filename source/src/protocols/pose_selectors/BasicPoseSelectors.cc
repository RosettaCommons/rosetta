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
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/filters/FilterFactory.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/sort_predicates.hh>
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

// C++ Headers
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.pose_selectors.BasicPoseSelectors" );

namespace protocols {
namespace pose_selectors {

////////////////////////////////////////////////////////////////////////
// LogicalSelector base class

// Selector
LogicalSelector::LogicalSelector()
{
}

LogicalSelector::~LogicalSelector()
{
}

protocols::rosetta_scripts::PoseSelectorFlags LogicalSelector::get_flags() const
{
	protocols::rosetta_scripts::PoseSelectorFlags flags(protocols::rosetta_scripts::PSF_NONE);
	BOOST_FOREACH ( protocols::rosetta_scripts::PoseSelectorOP selector, selectors_ ) {
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
	BOOST_FOREACH ( utility::tag::TagCOP const curr_tag, tag->getTags() ) {
		protocols::rosetta_scripts::PoseSelectorOP new_selector(
			protocols::rosetta_scripts::PoseSelectorFactory::get_instance()->
			newPoseSelector( curr_tag, data, filters, movers, pose )
		);
		runtime_assert( new_selector != 0 );
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

	BOOST_FOREACH ( protocols::rosetta_scripts::PoseSelectorOP selector, selectors_ ) {
		utility::vector1<bool> selector_selected_poses( selector->select_poses( poses ) );

		// Merge sets using logical operator
		for (
				utility::vector1<bool>::iterator
				i = selected_poses.begin(),
				j = selector_selected_poses.begin();
				i != selected_poses.end() &&
				j != selector_selected_poses.end();
				++i, ++j ) {

			(*i) = selection_operation(*i, *j);
		}

		TR.Debug << "Pose selections for " << get_name() << " after " << selector->get_name() << ": ";
		BOOST_FOREACH ( bool selection, selected_poses ) {
			TR.Debug << selection << " ";
		}
		TR.Debug << std::endl;
	}

	TR.Debug << "Final pose selections for " << get_name() << ": ";
	BOOST_FOREACH ( bool selection, selected_poses ) {
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

////////////////////////////////////////////////////////////////////////
// OrSelector

// Creator
protocols::rosetta_scripts::PoseSelectorOP OrSelectorCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new OrSelector() );
}

////////////////////////////////////////////////////////////////////////
// TopNByProperty

// Creator
protocols::rosetta_scripts::PoseSelectorOP TopNByPropertyCreator::create_selector() const {
	return protocols::rosetta_scripts::PoseSelectorOP( new TopNByProperty() );
}

// Selector
TopNByProperty::TopNByProperty() :
	reporter_(/* NULL */),
	order_(1)
{
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
	BOOST_FOREACH ( utility::tag::TagCOP const curr_tag, tag->getTags() ) {
		protocols::rosetta_scripts::PosePropertyReporterOP new_reporter(
			protocols::rosetta_scripts::PosePropertyReporterFactory::get_instance()->
			newPosePropertyReporter( curr_tag, data, filters, movers, pose )
		);
		runtime_assert( new_reporter != 0 );
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
		BOOST_FOREACH ( core::pose::PoseOP pose, poses ) {
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
	for ( utility::vector1<Pose_Property>::iterator it=pose_properties.begin(); it!=pose_properties.end(); ++it ) {
		Pose_Property &p = *it;
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

// Selector
Filter::Filter() :
	filter_(/* NULL */)
{
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

	TagCOP filter_tag(NULL);
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
	BOOST_FOREACH ( core::pose::PoseOP pose, poses ) {
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
