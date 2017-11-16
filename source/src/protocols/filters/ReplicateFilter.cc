// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/ReplicateFilter.cc
/// @brief Repeat a subfilter multiple times, and pass a value based on the aggregate results
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#include <protocols/filters/ReplicateFilter.hh>
#include <protocols/filters/ReplicateFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector1.hh>
#include <numeric/util.hh>

#include <algorithm>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace filters {

static basic::Tracer TR( "protocols.filters.ReplicateFilter" );

/// @brief default ctor
ReplicateFilter::ReplicateFilter() :
	Filter( "ReplicateFilter" ),
	replicates_(1),
	upper_trim_(0),
	lower_trim_(0),
	median_(false),
	threshold_(0)
{}

ReplicateFilter::ReplicateFilter(protocols::filters::FilterOP subfilter, core::Size replicates, core::Size upper_trim, core::Size lower_trim) :
	Filter( "ReplicateFilter" ),
	subfilter_(std::move(subfilter)),
	replicates_(replicates),
	upper_trim_(upper_trim),
	lower_trim_(lower_trim),
	median_(false),
	threshold_(0)
{}

bool
ReplicateFilter::apply(core::pose::Pose const & pose) const
{
	core::Real value = compute(pose);
	return value <= threshold_;
}

core::Real
ReplicateFilter::report_sm( core::pose::Pose const & pose) const
{
	return( compute(pose) );
}

void
ReplicateFilter::report( std::ostream & out, core::pose::Pose const & pose) const
{
	out<<"ReplicateFilter returns "<< compute(pose) <<std::endl;
}

core::Real
ReplicateFilter::compute(core::pose::Pose const & pose) const {
	runtime_assert( subfilter_ != nullptr );
	if ( lower_trim_ + upper_trim_ >= replicates_ ) {
		TR.Warning << "Replicate Filter trims off all values - returning 0." << std::endl;
	}
	utility::vector1<core::Real> values;
	for ( core::Size ii(1); ii <= replicates_; ++ii ) {
		values.push_back( subfilter_->report_sm( pose ) );
	}
	std::sort( values.begin(), values.end() );
	if ( values[1] == values[replicates_] && replicates_ != 1 ) { // Exact comparison of Real values is intentional
		TR.Warning << "Replicate Filter used, but all " << replicates_ << " replicates are identical! " << std::endl;
	}
	utility::vector1<core::Real> value_trim;
	TR << "Replicate filter: ";
	for ( core::Size jj(1); jj <= replicates_; ++jj ) {
		if ( jj == replicates_ - upper_trim_ + 1 ) { TR << "| "; }
		TR << values[jj] << " ";
		if ( jj == lower_trim_ ) { TR << "| "; }
		if ( jj >= lower_trim_ + 1 && jj <= replicates_ - upper_trim_ ) {
			value_trim.push_back( values[jj] );
		}
	}
	core::Real value;
	if ( median_ ) {
		TR << "=> median ";
		value = numeric::median( value_trim );
	} else {
		TR << "=> mean ";
		value = numeric::mean( value_trim );
	}
	TR << "=> " << value << std::endl;
	return value;
}

void
ReplicateFilter::parse_my_tag( utility::tag::TagCOP tag_ptr,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &)
{
	subfilter_ = protocols::rosetta_scripts::parse_filter( tag_ptr->getOption<std::string>( "filter_name" ) , filters );
	replicates_ = tag_ptr->getOption<core::Size>( "replicates", 1 );
	core::Real upper_cut = tag_ptr->getOption<core::Real>( "upper_cut", 0 );
	core::Real lower_cut = tag_ptr->getOption<core::Real>( "lower_cut", 0 );
	lower_trim_ = core::Size(lower_cut);
	if ( lower_cut < 1.0 ) {
		lower_trim_ = core::Size(lower_cut * replicates_); // Truncation is desired
	}
	upper_trim_ = core::Size(upper_cut);
	if ( upper_cut < 1.0 ) {
		upper_trim_ = core::Size(upper_cut * replicates_); // Truncation is desired
	}
	if ( lower_trim_ + upper_trim_ >= replicates_ ) {
		utility_exit_with_message("In ReplicateFilter, the number of items removed by upper and lower cuts exceed the total number of items.");
	}
	median_ = tag_ptr->getOption< bool >( "median", false );
	threshold_ = tag_ptr->getOption<core::Real>( "threshold", 0.0 );
}

// XRW TEMP FilterOP
// XRW TEMP ReplicateFilterCreator::create_filter() const { return FilterOP( new ReplicateFilter ); }

// XRW TEMP std::string
// XRW TEMP ReplicateFilterCreator::keyname() const { return "ReplicateFilter"; }

std::string ReplicateFilter::name() const {
	return class_name();
}

std::string ReplicateFilter::class_name() {
	return "ReplicateFilter";
}

void ReplicateFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute(
		"filter_name", xs_string,
		"Name of the filter to be applied multiple times." )
		+ XMLSchemaAttribute::attribute_w_default(
		"replicates", xsct_non_negative_integer,
		"Number of replicates.",
		"1" )
		+ XMLSchemaAttribute::attribute_w_default(
		"upper_cut", xsct_real,
		"Trim of the highest value bye upper_cut.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"lower_cut", xsct_real,
		"Trim of the lowest value bye lower_cut.",
		"0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"median", xsct_rosetta_bool,
		"If true, calculate the median instead of the average.",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"threshold", xsct_real,
		"Filter returns true of calculated value is less than threshold.",
		"0" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Repeat a filter multiple times and average.",
		attlist );
}

std::string ReplicateFilterCreator::keyname() const {
	return ReplicateFilter::class_name();
}

protocols::filters::FilterOP
ReplicateFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ReplicateFilter );
}

void ReplicateFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReplicateFilter::provide_xml_schema( xsd );
}


} // filters
} // protocols
