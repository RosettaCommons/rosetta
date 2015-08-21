// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/LeastNativeLike9merFilter
/// @brief filter structures by IntraRepeatContacts
/// @details
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/LeastNativeLike9merFilter.hh>
#include <protocols/simple_filters/LeastNativeLike9merFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoringManager.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

// vall_lookback
#include <core/scoring/methods/vall_lookback/VallLookbackPotential.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackData.hh>

//// C++ headers
static thread_local basic::Tracer tr("protocols.filters.LeastNativeLike9merFilter");

namespace protocols {
namespace simple_filters {

// @brief default constructor
LeastNativeLike9merFilter::LeastNativeLike9merFilter():
	Filter( "worst9mer" ),
	filtered_value_( 99 )
{}

// @brief copy constructor
LeastNativeLike9merFilter::LeastNativeLike9merFilter( LeastNativeLike9merFilter const & rval ):
	Super( rval ),
	filtered_value_( rval.filtered_value_ )
{}

// @brief destructor
LeastNativeLike9merFilter::~LeastNativeLike9merFilter() {}

// @brief set filtered value
void LeastNativeLike9merFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
LeastNativeLike9merFilter::Real
LeastNativeLike9merFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
}

/// @brief
void
LeastNativeLike9merFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "worst_9mer: " <<  compute( pose ) << std::endl;
}

/// @brief
LeastNativeLike9merFilter::Real
LeastNativeLike9merFilter::compute( const Pose & pose ) const
{
	using namespace core::scoring::methods;
	using namespace core::scoring;

	VallLookbackPotential const & potential_( ScoringManager::get_instance()->get_vallLookbackPotential());
	Real worst_9mer = potential_.lookback(pose);
	return(worst_9mer);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool LeastNativeLike9merFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	std::cout << "value" << value << "filtered_value_" << filtered_value_ << std::endl;
	if ( value < filtered_value_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
LeastNativeLike9merFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", 99.0 );
	tr << "Structures which has the best fragment RMSD at the worst position greater than " << filtered_value_ << " will be filtered." << std::endl;
}

filters::FilterOP
LeastNativeLike9merFilterCreator::create_filter() const { return protocols::filters::FilterOP(new LeastNativeLike9merFilter); }

std::string
LeastNativeLike9merFilterCreator::keyname() const { return "worst9mer"; }
} // filters
} // protocols
