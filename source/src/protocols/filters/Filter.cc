// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/Filter.cc
/// @brief
/// @details
///   Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <basic/Tracer.hh>

#include <utility>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.filters.Filter" );

namespace protocols {
namespace filters {

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
typedef utility::tag::TagCOP TagCOP;

FilterCollection::~FilterCollection() = default;

bool
FilterCollection::apply( core::pose::Pose const & pose ) const
{
	for ( protocols::filters::FilterCOP filter : filters_ ) {
		if ( ! filter->apply( pose ) ) {
			return false;
		}
	}

	return true;
}

void
FilterCollection::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	for ( protocols::filters::FilterCOP filter : filters_ ) {
		filter->report( out, pose );
	}
}

Filter::Filter()
: utility::pointer::ReferenceCount(),
	type_( "UNDEFINED TYPE" ),
	scorename_("defaultscorename")
{}

Filter::Filter( std::string  type )
: utility::pointer::ReferenceCount(),
	type_(std::move( type )),
	scorename_("defaultscorename")
{}

Filter::Filter( Filter const & init )
: utility::pointer::ReferenceCount(),
	type_( init.type_ ),
	user_defined_name_( init.user_defined_name_ ),
	scorename_("defaultscorename")

{}

Filter::~Filter() = default;

void
Filter::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & )
{}

core::Real Filter::score( core::pose::Pose & pose ) {
	core::Real score = report_sm( pose );
	core::pose::setPoseExtraScore( pose, scorename_, score );
	return score;
}

} // filters
} // protocols
