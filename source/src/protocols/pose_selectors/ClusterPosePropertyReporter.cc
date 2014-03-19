// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/rosetta_scripts/ComparisonPoseReporters.hh
/// @brief  Collection of simple pose property reporters
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_ComparisonPoseReporters_CC
#define INCLUDED_protocols_rosetta_scripts_ComparisonPoseReporters_CC

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>
#include <protocols/rosetta_scripts/BasicPosePropertyReporters.hh>
#include <protocols/rosetta_scripts/BasicPosePropertyReporterCreators.hh>

// Package headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/filters/FilterFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <string>

static basic::Tracer TR( "protocols.rosetta_scripts.BasicPosePropertyReporters" );

namespace protocols {
namespace rosetta_scripts {

////////////////////////////////////////////////////////////////////////
// EnergyReporter

// Creator
PosePropertyReporterOP EnergyReporterCreator::create_reporter() const {
  return new EnergyReporter();
}

// Reporter
EnergyReporter::EnergyReporter() :
	scorefxn_(NULL),
	term_("")
{
}

core::Real EnergyReporter::report_property( core::pose::Pose & p ) const
{
	core::Real r = 0;

	if(term_ != "")	{
		core::scoring::ScoreType t( core::scoring::ScoreTypeManager::score_type_from_name(term_) );
		r = scorefxn_->score_by_scoretype(p, t);
	} else
		r = scorefxn_->score(p);

	return r;
}

void EnergyReporter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	scorefxn_ = tag->hasOption("scorefunction") ?
		core::scoring::ScoreFunctionFactory::create_score_function( tag->getOption<std::string>("scorefunction") ) :
		core::scoring::getScoreFunction();

	term_ = tag->hasOption("term") ?
		tag->getOption<std::string>("term") :
		"";
}


////////////////////////////////////////////////////////////////////////
// FilterReporter

// Creator
PosePropertyReporterOP FilterReporterCreator::create_reporter() const {
  return new FilterReporter();
}

// Reporter
FilterReporter::FilterReporter() :
	filter_(NULL)
{
}

core::Real FilterReporter::report_property( core::pose::Pose & p ) const
{
	core::Real r = 0;

	if(filter_) {
		r = filter_->report_sm(p);
	} else
		TR << "No filter instance; cannot score pose!" << std::endl;
	
	return r;
}

void FilterReporter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	using namespace utility::tag;
	using namespace protocols::rosetta_scripts;

	TagCOP filter_tag(NULL);

	if(tag->hasOption("filter")) {
		// Find a filter by name defined somewhere upstream in the script
		std::string filter_name( tag->getOption<std::string>("filter") );
		RosettaScriptsParser parser;
		filter_tag = parser.find_rosettascript_tag(
				tag,
				"FILTERS",
				"name",
				filter_name
		);

		if(!filter_tag)
			throw utility::excn::EXCN_RosettaScriptsOption("Cannot find filter named \"" + filter_name + "\"");

	} else {
		// Filter is defined inline (first child tag)
		utility::vector0< TagCOP > tags( tag->getTags() );
		for(utility::vector0< TagCOP >::const_iterator it = tags.begin(); it != tags.end(); ++it) {
			filter_tag = *it;
			break;
		}
	}

	if(filter_tag)
		filter_  = protocols::filters::FilterFactory::get_instance()->newFilter( filter_tag, data, filters, movers, pose );

	if(!filter_)
		std::ostringstream s;
		s << "Cannot create filter from script tag: " << tag;
		throw utility::excn::EXCN_RosettaScriptsOption(s.str());
	}
}

////////////////////////////////////////////////////////////////////////

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_ComparisonPoseReporters_CC
