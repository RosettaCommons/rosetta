// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/FragQualFilter.cc
/// @brief filter structures by packstat score
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/FragQualFilter.hh>
#include <protocols/fldsgn/filters/FragQualFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <core/scoring/packstat/compute_sasa.hh>
#include <basic/MetricValue.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/toolbox/pose_metric_calculators/FragQualCalculator.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//// C++ headers
static basic::Tracer tr("protocols.fldsgn.filters.FragQualFilter");

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
FragQualFilter::FragQualFilter():
	Filter( "FragQual" ),
	filtered_type_( "num_goodfrag" ),
	filtered_value_( 0.0 )
	// rmsd_cutoff_( 1.0 )
{}

// @brief copy constructor
FragQualFilter::FragQualFilter( FragQualFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_type_( rval.filtered_type_ ),
	filtered_value_( rval.filtered_value_ )
	// rmsd_cutoff_( 1.0 )
{}

// @brief set filtered value
void FragQualFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

// @brief set report type
void FragQualFilter::filtered_type( String const & value )
{
	filtered_type_ = value;
}

/// @brief
FragQualFilter::Real
FragQualFilter::report_sm( Pose const & pose ) const
{
	return 	compute( pose );
}

/// @brief
void
FragQualFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "FragQual: " <<  compute( pose ) << std::endl;
}

/// @brief
FragQualFilter::Real
FragQualFilter::compute( Pose const & pose ) const
{
	basic::MetricValue< Real > score;
	pose.metric( "FragQual", filtered_type_, score );
	return score.value();
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool FragQualFilter::apply( Pose const & pose ) const
{
	Real value = compute( pose );
	if( value > filtered_value_ ){
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	}else{
		tr << "Filter failed current/threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
FragQualFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose )
{
	using core::pose::metrics::CalculatorFactory;
	using protocols::toolbox::pose_metric_calculators::FragQualCalculator;

	// set filtered type
 	filtered_type_ = tag->getOption<String>( "type", "num_goodfrag" );
	if( filtered_type_ != "num_goodfrag" && filtered_type_ != "coverage" ) {
		tr << "Filter type, " << filtered_type_ <<  " is not defined." << std::endl;
		runtime_assert( false );
	}

	// set threshold
 	filtered_value_ = tag->getOption<Real>( "threshold", 0.0 );
	tr << "Structures with fragqual value, " << filtered_type_ << " above " << filtered_value_ << " will be filtred." << std::endl;

	// set FragQual
	CalculatorFactory::Instance().remove_calculator( "FragQual" );
	FragQualCalculator calculator;
	calculator.parse_my_tag( tag, data, filters, movers, pose  );
	CalculatorFactory::Instance().register_calculator( "FragQual", calculator.clone() );

	//calculator.begin( tag->getOption<Size>( "begin", 1 ) );
	//calculator.end( tag->getOption<Size>( "end", pose.total_residue() ) );
	//
	//rmsd_cutoff_ = tag->getOption<Real>( "rmsd_cutoff", 1.0 );
	//calculator.rmsd_cutoff( rmsd_cutoff_ );

}

protocols::filters::FilterOP
FragQualFilterCreator::create_filter() const { return new FragQualFilter; }

std::string
FragQualFilterCreator::keyname() const { return "FragQual"; }


} // filters
} // fldsgn
} // protocols
