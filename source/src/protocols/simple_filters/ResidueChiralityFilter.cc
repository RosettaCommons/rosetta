// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueChiralityFilter.cc
/// @brief checks the chirality of a specific residues, whether it is D or L
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
// Project Headers

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/simple_filters/ResidueChiralityFilter.hh>
#include <protocols/simple_filters/ResidueChiralityFilterCreator.hh>
#include <string>
#include <utility/tag/Tag.hh>

#include <fstream>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue::ResidueType

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.ResidueChiralityFilter" );

protocols::filters::FilterOP
ResidueChiralityFilterCreator::create_filter() const { return utility::pointer::make_shared< ResidueChiralityFilter >(); }

std::string
ResidueChiralityFilterCreator::keyname() const { return "ResidueChirality"; }

ResidueChiralityFilter::~ResidueChiralityFilter(){}

void
ResidueChiralityFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {
	required_type_ = tag->getOption< std::string >( "residue_type", "L" );
	res_num_ = tag->getOption< core::Size >( "res_num", 1 );
}

bool
ResidueChiralityFilter::apply( core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	bool const status = type_found == required_type_;
	return status;
}

void
ResidueChiralityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	out << "found residue " << res_num_ << " type to be " << type_found << std::endl;
}

core::Real
ResidueChiralityFilter::report_sm( core::pose::Pose const & pose ) const {
	std::string type_found = compute( pose );
	return( type_found == required_type_ );
}

std::string
ResidueChiralityFilter::compute( core::pose::Pose const & pose ) const {
	bool is_d = pose.residue(res_num_).type().is_d_aa();
	std::string residue_type = "L";
	if ( is_d ) residue_type = "D";
	return( residue_type );
}

std::string ResidueChiralityFilter::name() const {
	return class_name();
}

std::string ResidueChiralityFilter::class_name() {
	return "ResidueChirality";
}

void ResidueChiralityFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "residue_type" , xs_string , "the type you want to pass" , "L" )
		+ XMLSchemaAttribute::attribute_w_default( "res_num" , xs_integer , "which residue to test" , "1" ) ;
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "fails if the res_num residue is in the wrong chirality (D / L)", attlist );
}

void ResidueChiralityFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueChiralityFilter::provide_xml_schema( xsd );
}

}
}
