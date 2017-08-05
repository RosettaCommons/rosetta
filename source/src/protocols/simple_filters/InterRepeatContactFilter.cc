// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterRepeatContactFilter
/// @brief filter structures by InterRepeatContacts
/// @details
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/InterRepeatContactFilter.hh>
#include <protocols/simple_filters/InterRepeatContactFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>

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
static THREAD_LOCAL basic::Tracer tr("protocols.filters.InterRepeatContactFilter");

namespace protocols {
namespace simple_filters {

// @brief default constructor
InterRepeatContactFilter::InterRepeatContactFilter():
	Filter( "InterRepeatContactsPerResidue" ),
	filtered_value_( 0.0 ),
	numbRepeats_(4),
	sequenceSep_(6),
	distThresh_(10)
{}

// @brief copy constructor
InterRepeatContactFilter::InterRepeatContactFilter( InterRepeatContactFilter const & )= default;

// @brief set filtered value
void InterRepeatContactFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
InterRepeatContactFilter::Real
InterRepeatContactFilter::report_sm( Pose const & pose ) const
{
	return  compute( pose );
}

/// @brief
void
InterRepeatContactFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "InterRepeatContacts: " <<  compute( pose ) << std::endl;
}

/// @brief
InterRepeatContactFilter::Real
InterRepeatContactFilter::compute( Pose const & pose ) const
{
	using numeric::xyzVector;
	Size repeatLength = pose.size()/numbRepeats_;
	Size region1Start = 1;
	Size region1End = repeatLength;
	Size region2Start = repeatLength+1;
	Size region2End = repeatLength*2;
	Size contactCt = 0;
	for ( Size ii=region1Start; ii<=region1End; ++ii ) {
		for ( Size jj=region2Start; jj<=region2End; ++jj ) {
			if ( (jj-ii)>sequenceSep_ ) {
				xyzVector<double> a = pose.xyz(core::id::NamedAtomID("CA", ii));
				xyzVector<double> b = pose.xyz(core::id::NamedAtomID("CA", jj));
				if ( a.distance(b)<distThresh_ ) {
					contactCt+=1;
				}
			}
		}
	}
	return((Real(contactCt))/Real(repeatLength*2));
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool InterRepeatContactFilter::apply( Pose const & pose ) const
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
InterRepeatContactFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", -1.0 );
	numbRepeats_ = tag->getOption<Size>("numb_repeats",4);
	sequenceSep_ = tag->getOption<Size>("sequenceSeperation",6);
	distThresh_ = tag->getOption<Real>("distanceThreshold",10.0);
	tr << "Structures which have InterRepeatContacts less than " << filtered_value_ << " will be filtered." << std::endl;
}

// XRW TEMP filters::FilterOP
// XRW TEMP InterRepeatContactFilterCreator::create_filter() const { return protocols::filters::FilterOP(new InterRepeatContactFilter); }

// XRW TEMP std::string
// XRW TEMP InterRepeatContactFilterCreator::keyname() const { return "InterRepeatContactsPerResidue"; }

std::string InterRepeatContactFilter::name() const {
	return class_name();
}

std::string InterRepeatContactFilter::class_name() {
	return "InterRepeatContactsPerResidue";
}

void InterRepeatContactFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "filter threshold", "-1.0")
		+ XMLSchemaAttribute::attribute_w_default("numb_repeats", xsct_non_negative_integer, "number of repeats", "4")
		+ XMLSchemaAttribute::attribute_w_default("sequenceSeperation", xsct_non_negative_integer, "seperation in the sequence", "6")
		+ XMLSchemaAttribute::attribute_w_default("distanceThreshold", xsct_real, "max distance between two CA atoms being compared", "10.0");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "filter based upon the inter repeat contacts", attlist );
}

std::string InterRepeatContactFilterCreator::keyname() const {
	return InterRepeatContactFilter::class_name();
}

protocols::filters::FilterOP
InterRepeatContactFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new InterRepeatContactFilter );
}

void InterRepeatContactFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterRepeatContactFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
