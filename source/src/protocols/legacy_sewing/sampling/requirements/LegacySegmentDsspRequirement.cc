// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacySegmentDsspRequirement.cc
///
/// @brief
/// @author Tim Jacobs, Frank Teets

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacySegmentDsspRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacySegmentDsspRequirementCreator.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//Utility headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <string>
#include <basic/Tracer.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

static THREAD_LOCAL basic::Tracer TR("protocols.legacy_sewing.sampling.requirements.LegacySegmentDsspRequirement");

//////Creator methods/////
LegacyIntraSegmentRequirementOP
LegacySegmentDsspRequirementCreator::create_requirement() const{
	return LegacyIntraSegmentRequirementOP( new LegacySegmentDsspRequirement() );
}

std::string
LegacySegmentDsspRequirementCreator::type_name() const {
	return "LegacySegmentDsspRequirement";
}
//////End Creator methods/////

LegacySegmentDsspRequirement::LegacySegmentDsspRequirement(){}

LegacySegmentDsspRequirement::LegacySegmentDsspRequirement(
	std::set<std::string> valid_dssp_codes
):
	valid_dssp_codes_(valid_dssp_codes)
{}

bool
LegacySegmentDsspRequirement::satisfies(
	SewSegment segment
) const {
	return !violates(segment);
}

bool
LegacySegmentDsspRequirement::violates(
	SewSegment segment
) const {
	if ( valid_dssp_codes_.find(std::string(1, segment.dssp_)) == valid_dssp_codes_.end() ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Violation! " << segment.dssp_<< " is invalid!" << std::endl; }
		return true;
	}
	return false;
}


void
LegacySegmentDsspRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
) {
	if ( !tag->hasOption("dssp") ) {
		utility_exit_with_message("LegacySegmentDsspRequiremnt must be given a whitespace delimited 'dssp' option");
	}
	valid_dssp_codes_ = utility::split_to_set(tag->getOption<std::string>("dssp"));
}

void
LegacySegmentDsspRequirement::show(
	std::ostream & out
) const {
	out << "/////// LegacySegmentDsspRequirement - Segment must have on of the following DSSP codes: ";
	std::set<std::string>::const_iterator it = valid_dssp_codes_.begin();
	std::set<std::string>::const_iterator it_end = valid_dssp_codes_.end();
	for ( ; it != it_end; ++it ) {
		out << *it << " ";
	}
	out << std::endl;
}

void
LegacySegmentDsspRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LegacySegmentDsspRequirement::provide_xml_schema( xsd );
}


std::string
LegacySegmentDsspRequirement::class_name(){
	return "LegacySegmentDsspRequirement";
}

void
LegacySegmentDsspRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	//min_length and max_length
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "dssp", xs_string, "Whitespace-delimited list of valid DSSP characters for this segment" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyRequirementFactory::legacy_intra_segment_requirements_ct_namer )
		.element_name( class_name() )
		.description( "Sets allowed DSSP for a segment" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}










} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace
