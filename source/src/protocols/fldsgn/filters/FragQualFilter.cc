// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/FragQualFilter.cc
/// @brief filter structures by packstat score
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/FragQualFilter.hh>
#include <protocols/fldsgn/filters/FragQualFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/pose_metric_calculators/FragQualCalculator.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer tr( "protocols.fldsgn.filters.FragQualFilter" );

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
FragQualFilter::FragQualFilter( FragQualFilter const & /*rval*/ ) = default;

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
	return  compute( pose );
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
	if ( value > filtered_value_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
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
	using protocols::pose_metric_calculators::FragQualCalculator;

	// set filtered type
	filtered_type_ = tag->getOption<String>( "type", "num_goodfrag" );
	if ( filtered_type_ != "num_goodfrag" && filtered_type_ != "coverage" ) {
		tr.Fatal << "Filter type, " << filtered_type_ <<  " is not defined." << std::endl;
		utility_exit_with_message("Filter type not defined.");
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
	//calculator.end( tag->getOption<Size>( "end", pose.size() ) );
	//
	//rmsd_cutoff_ = tag->getOption<Real>( "rmsd_cutoff", 1.0 );
	//calculator.rmsd_cutoff( rmsd_cutoff_ );

}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP FragQualFilterCreator::create_filter() const { return protocols::filters::FilterOP( new FragQualFilter ); }

// XRW TEMP std::string
// XRW TEMP FragQualFilterCreator::keyname() const { return "FragQual"; }

std::string FragQualFilter::name() const {
	return class_name();
}

std::string FragQualFilter::class_name() {
	return "FragQual";
}

void FragQualFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction frag_qual_filter_type;
	frag_qual_filter_type.name( "frag_qual_filter_type");
	frag_qual_filter_type.base_type( xs_string );
	frag_qual_filter_type.add_restriction( xsr_enumeration, "num_goodfrag" );
	frag_qual_filter_type.add_restriction( xsr_enumeration, "coverage" );
	xsd.add_top_level_element( frag_qual_filter_type );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "type", "frag_qual_filter_type", "XRW TO DO", "num_goodfrag" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "XRW TO DO", "0.0" )
		//add all attributes detected in FragQUalCalculator's parse_my_tag
		+ XMLSchemaAttribute::attribute_w_default( "ratio_cutoff", xsct_real, "XRW TO DO", "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "begin", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute( "end", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "rmsd_cutoff", xsct_real, "XRW TO DO", "1.0" )
		+ XMLSchemaAttribute::required_attribute( "frag", xs_string, "XRW TO DO" ) //refers to a fragset in the data map
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "XRW TO DO", "false" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string FragQualFilterCreator::keyname() const {
	return FragQualFilter::class_name();
}

protocols::filters::FilterOP
FragQualFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new FragQualFilter );
}

void FragQualFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragQualFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
