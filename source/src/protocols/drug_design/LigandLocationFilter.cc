// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/LigandLocationFilter.cc
/// @brief
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

//Unit Headers
#include <protocols/drug_design/LigandLocationFilter.hh>
#include <protocols/drug_design/LigandLocationFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/selection.hh>
#include <core/types.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <fstream>

namespace protocols {
namespace drug_design {

static basic::Tracer TR( "protocols.drug_design.LigandLocationFilter" );

protocols::filters::FilterOP
LigandLocationFilterCreator::create_filter() const { return protocols::filters::FilterOP( new LigandLocationFilter ); }

std::string
LigandLocationFilterCreator::keyname() const { return LigandLocationFilter::class_name(); }

void LigandLocationFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	LigandLocationFilter::provide_xml_schema( xsd );
}

std::string LigandLocationFilter::class_name() {
	return "LigandLocationFilter";
}

void LigandLocationFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "chain", xs_string, "Chain ID" )
		+ XMLSchemaAttribute( "radius", xsct_real,
		"Filter is false if the distance between ligand centroid and center is greater than this value" );

	// subelements
	XMLSchemaSimpleSubelementList subelement_list;
	AttributeList Coordinates_attributes;
	Coordinates_attributes
	+ XMLSchemaAttribute::required_attribute("x", xsct_real, "x coordinate of the center.")
	+ XMLSchemaAttribute::required_attribute("y", xsct_real, "y coordinate of the center.")
	+ XMLSchemaAttribute::required_attribute("z", xsct_real, "z coordinate of the center.");
	subelement_list.add_simple_subelement("center", Coordinates_attributes, "Set x,y and z coordinates, center of filter.");

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(),
		"Test whether the ligand is within the pocket region.",
		attlist, subelement_list );
}

void
LigandLocationFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	chain_ = tag->getOption< std::string >( "chain" );
	radius_ = tag->getOption<core::Real>( "radius", radius_ );

	for ( utility::tag::TagCOP child_tag : tag->getTags() ) {
		std::string name= child_tag->getName();
		if ( name == "center" ) {
			center_ = core::Vector(
					child_tag->getOption<core::Real>("x"),
					child_tag->getOption<core::Real>("y"),
					child_tag->getOption<core::Real>("z")
			);
		}
	}
}

bool
LigandLocationFilter::apply( core::pose::Pose const & pose ) const {
	core::Real distance = compute_distance( pose );
	if ( distance > radius_ ) {
		TR << "Failing LigandLocation filter with distance " << distance << " between centroid and center, which is outside the radius of " << radius_ << std::endl;
		return false;
	}
	return true;
}

void
LigandLocationFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const distance( compute_distance( pose ) );
	out << "The distance between the centroid of the ligand and the set center is " << distance << std::endl;
}

core::Real
LigandLocationFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute_distance( pose ) );
}

/*
 * Compute the distance between the centroid of the ligand and a set center
 */
core::Real
LigandLocationFilter::compute_distance( core::pose::Pose const &pose ) const {
	core::Size chain_id = core::pose::get_chain_id_from_chain( chain_, pose );
	core::Vector const centroid = protocols::geometry::centroid_by_chain( pose, chain_id );
	core::Real const distance = ( center_ - centroid ).length();
	return distance;
}

}
}
