// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacySegmentLengthRequirement.cc
///
/// @brief
/// @author Tim Jacobs, Frank Teets

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacySegmentLengthRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacySegmentLengthRequirementCreator.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//Utility headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.legacy_sewing.sampling.requirements.LegacySegmentLengthRequirement");

//////Creator methods/////
LegacyIntraSegmentRequirementOP
LegacySegmentLengthRequirementCreator::create_requirement() const{
	return LegacyIntraSegmentRequirementOP( new LegacySegmentLengthRequirement() );
}

std::string
LegacySegmentLengthRequirementCreator::type_name() const {
	return "LegacySegmentLengthRequirement";
}
//////End Creator methods/////

LegacySegmentLengthRequirement::LegacySegmentLengthRequirement():
	min_length_(1),
	max_length_(1)
{}

LegacySegmentLengthRequirement::LegacySegmentLengthRequirement(
	core::Size min_length,
	core::Size max_length
):
	min_length_(min_length),
	max_length_(max_length)
{}

//@brief segments do not change; if the segment does not violate the requirement it must satisfy it
bool
LegacySegmentLengthRequirement::satisfies(
	SewSegment segment
) const {
	return !violates(segment);
}

//@brief check if the segment falls outside the range parameters
bool
LegacySegmentLengthRequirement::violates(
	SewSegment segment
) const {
	core::Size length = segment.residues_.size();
	if ( length < min_length_ ) {
		//if ( TR.Debug.visible() ) { TR.Debug << length << " Violation! Too Short!" << std::endl; }
		if ( TR.Debug.visible() ) { TR.Debug << " Violation! " << length << " is shorter than " << min_length_ << std::endl; }
		return true;
	} else if ( length > max_length_ ) {
		//if ( TR.Debug.visible() ) { TR.Debug << length << " Violation! Too Long!" << std::endl; }
		if ( TR.Debug.visible() ) { TR.Debug << " Violation! " << length << " is larger than " << max_length_ << std::endl; }
		return true;
	}
	return false;
}

void
LegacySegmentLengthRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
) {
	if ( !tag->hasOption("min_length") && !tag->hasOption("max_length") ) {
		utility_exit_with_message("You must provide either max_length, min_length, or both to the LegacySegmentLengthRequirement tag!");
	} else if ( tag->hasOption("min_length") && !tag->hasOption("max_length") ) {
		min_length_ = tag->getOption<core::Size>("min_length");
		max_length_ = -1;/*max size*/
	} else if ( !tag->hasOption("min_length") && tag->hasOption("max_length") ) {
		max_length_ = tag->getOption<core::Size>("max_length");
		min_length_ = 0;
	} else {
		min_length_ = tag->getOption<core::Size>("min_length");
		max_length_ = tag->getOption<core::Size>("max_length");
	}
}

void
LegacySegmentLengthRequirement::show(
	std::ostream & out
) const {
	out << "/////// LegacySegmentLengthRequirement - Segment must contain between "
		<< min_length_ << " and " << max_length_ << " residues" << std::endl;
}

void
LegacySegmentLengthRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LegacySegmentLengthRequirement::provide_xml_schema( xsd );
}


std::string
LegacySegmentLengthRequirement::class_name(){
	return "LegacySegmentLengthRequirement";
}

void
LegacySegmentLengthRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	//min_length and max_length
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "max_length", xsct_non_negative_integer, "Maximum length of segment in residues" )
		+ XMLSchemaAttribute( "min_length", xsct_non_negative_integer, "Minimum length of segment in residues" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyRequirementFactory::legacy_intra_segment_requirements_ct_namer )
		.element_name( class_name() )
		.description( "Sets limits on minimum and maximum length of a segment" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}


} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace
