// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyGlobalLengthRequirement.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalLengthRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalLengthRequirementCreator.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Utility headers
#include <utility>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.legacy_sewing.sampling.requirements.LegacyGlobalLengthRequirement");

//////Creator methods/////
LegacyGlobalRequirementOP
LegacyGlobalLengthRequirementCreator::create_requirement() const{
	return LegacyGlobalRequirementOP( new LegacyGlobalLengthRequirement() );
}

std::string
LegacyGlobalLengthRequirementCreator::type_name() const {
	return "LegacyGlobalLengthRequirement";
}
//////End Creator methods/////

LegacyGlobalLengthRequirement::LegacyGlobalLengthRequirement():
	min_length_(0),
	max_length_(-1)
{}

LegacyGlobalLengthRequirement::LegacyGlobalLengthRequirement(
	std::set<std::string> const & valid_dssp_codes,
	core::Size min_length,
	core::Size max_length
):
	dssp_codes_(valid_dssp_codes),
	min_length_(min_length),
	max_length_(max_length)
{}

bool
LegacyGlobalLengthRequirement::satisfies(
	AssemblyCOP assembly
) const {
	return !violates(assembly);
}

bool
LegacyGlobalLengthRequirement::violates(
	AssemblyCOP assembly
) const {
	utility::vector1<SewSegment> const & segments = assembly->segments();
	for ( core::Size i = 1; i <= segments.size(); ++i ) {
		if ( dssp_codes_.size() == 0 || dssp_codes_.find(std::string(1, segments[i].dssp_)) != dssp_codes_.end() ) {
			if ( segments[i].residues_.size() < min_length_ || segments[i].residues_.size() > max_length_ ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "Segment " << i << ", with dssp " << segments[i].dssp_ << " and "
						<< segments[i].residues_.size() << " residues fails LegacyGlobalLengthRequirement" << std::endl;
				}
				return true;
			}
		}
	}
	return false;
}

void
LegacyGlobalLengthRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	if ( tag->hasOption("dssp") ) {
		dssp_codes_ = utility::split_to_set(tag->getOption<std::string>("dssp"));
	}

	if ( tag->hasOption("min_length") ) {
		min_length_ = tag->getOption<core::Size>("min_length");
	}

	if ( tag->hasOption("max_length") ) {
		max_length_ = tag->getOption<core::Size>("max_length");
	}
}


void
LegacyGlobalLengthRequirement::show(
	std::ostream & out
) const {
	out << "/////// LegacyGlobalLengthRequirement - ";
	if ( dssp_codes_.size() > 0 ) {
		out << "Segments with DSSP codes ";
		auto it = dssp_codes_.begin();
		auto it_end = dssp_codes_.end();
		for ( ; it != it_end; ++it ) {
			out << *it << ",";
		}
	} else {
		out << "All segments";
	}
	out << " must have between " << min_length_ << " and " << max_length_ << " residues." << std::endl;
}

std::string
LegacyGlobalLengthRequirement::class_name(){
	return "LegacyGlobalLengthRequirement";
}

void
LegacyGlobalLengthRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LegacyGlobalLengthRequirement::provide_xml_schema( xsd );
}

void
LegacyGlobalLengthRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	//min_length and max_length
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "max_length", xsct_non_negative_integer, "Maximum length of each segment in residues" )
		+ XMLSchemaAttribute( "min_length", xsct_non_negative_integer, "Minimum length of each segment in residues" )
		+ XMLSchemaAttribute( "dssp", xs_string, "Whitespace-delimited list of valid DSSP characters for which to apply this requirement" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyRequirementFactory::legacy_global_requirements_ct_namer )
		.element_name( class_name() )
		.description( "Sets limits on minimum and maximum length of each segment with a given DSSP" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}


} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace
