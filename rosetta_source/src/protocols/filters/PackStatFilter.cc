// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/PackStatFilter.cc
/// @brief filter structures by packstat score
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/filters/PackStatFilter.hh>
#include <protocols/filters/PackStatFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

//Auto Headers


//// C++ headers
static basic::Tracer tr("protocols.filters.PackStatFilter");

namespace protocols {
namespace filters {

// @brief default constructor
PackStatFilter::PackStatFilter():
	Filter( "PackStat" ),
	filtered_score_( 0.58 )  // ideally, ~0.65 is required for good packing
{}


// @brief constructor with arguments
PackStatFilter::PackStatFilter( Real const & score ):
	Filter( "PackStat" ),
	filtered_score_( score )
{}

// @brief copy constructor
PackStatFilter::PackStatFilter( PackStatFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_score_( rval.filtered_score_ )
{}

// @brief set filtered secondary structure
void PackStatFilter::filtered_score( Real const & score )
{
	filtered_score_ = score;
}

/// @brief
PackStatFilter::Real
PackStatFilter::compute( Pose const & pose ) const
{
	// calc packstat
	core::Real packscore;
	packscore = core::scoring::packstat::compute_packing_score( pose );
	return packscore;
}

/// @brief
PackStatFilter::Real
PackStatFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool PackStatFilter::apply( Pose const & pose ) const
{
	Real score = compute( pose );
	if( score > filtered_score_ ){
		tr << "Successfully filtered: " << score << std::endl;
		return true;
	}else{
		tr << "Filter failed current/threshold=" << score << "/" << filtered_score_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
PackStatFilter::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
 	filtered_score_ = tag->getOption<Real>( "threshold", 0.58 ); // ideally, ~0.65 is required for good packing
	tr << "Structures with packstat score " << filtered_score_ << " will be filtred." << std::endl;
}

FilterOP
PackStatFilterCreator::create_filter() const { return new PackStatFilter; }

std::string
PackStatFilterCreator::keyname() const { return "PackStat"; }


} // filters
} // protocols
