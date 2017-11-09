// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueDistanceFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/ResidueDistanceFilter.hh>
#include <protocols/simple_filters/ResidueDistanceFilterCreator.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer residue_distance_filter_tracer( "protocols.simple_filters.ResidueDistanceFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ResidueDistanceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueDistanceFilter ); }

// XRW TEMP std::string
// XRW TEMP ResidueDistanceFilterCreator::keyname() const { return "ResidueDistance"; }

ResidueDistanceFilter::~ResidueDistanceFilter()= default;

void
ResidueDistanceFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	res1_ = core::pose::get_resnum_string( tag, "res1_" );
	res2_ = core::pose::get_resnum_string( tag, "res2_" );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );

	residue_distance_filter_tracer<<"ResidueDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<<res1_<<" and "<<res2_<<std::endl;
}

bool
ResidueDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	residue_distance_filter_tracer<<"Distance between residues "<<res1_<<" and "<<res2_<<" is "<<distance<<std::endl;
	return( distance<=distance_threshold_ );
}

void
ResidueDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	out<<"Distance between residues "<<res1_<<" and "<<res2_<<" is "<<distance<<'\n';
}

core::Real
ResidueDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	return( distance );
}
core::Real
ResidueDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::Size res1( core::pose::parse_resnum( res1_, pose ) );
	core::Size res2( core::pose::parse_resnum( res2_, pose ) );
	core::conformation::Residue const res_res1( pose.conformation().residue( res1 ) );
	core::conformation::Residue const res_res2( pose.conformation().residue( res2 ) );
	core::Real const distance( res_res1.xyz( res_res1.nbr_atom() ).distance( res_res2.xyz( res_res2.nbr_atom() ) ) );
	return( distance );
}

std::string ResidueDistanceFilter::name() const {
	return class_name();
}

std::string ResidueDistanceFilter::class_name() {
	return "ResidueDistance";
}

void ResidueDistanceFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"distance", xsct_real,
		"Maximal distance between residue 1 and residue 2");

	core::pose::attributes_for_get_resnum_string(attlist, "res1_");
	core::pose::attributes_for_get_resnum_string(attlist, "res2_");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Returns true if the distance between 2 residues does not "
		"exceed the threshold",
		attlist );
}

std::string ResidueDistanceFilterCreator::keyname() const {
	return ResidueDistanceFilter::class_name();
}

protocols::filters::FilterOP
ResidueDistanceFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ResidueDistanceFilter );
}

void ResidueDistanceFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueDistanceFilter::provide_xml_schema( xsd );
}


}
}
