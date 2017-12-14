// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyResidueRetentionRequirement.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyResidueRetentionRequirement.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyResidueRetentionRequirementCreator.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Utility headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.legacy_sewing.sampling.requirements.LegacyResidueRetentionRequirement");

//////Creator methods/////
LegacyGlobalRequirementOP
LegacyResidueRetentionRequirementCreator::create_requirement() const{
	return LegacyGlobalRequirementOP( new LegacyResidueRetentionRequirement() );
}

std::string
LegacyResidueRetentionRequirementCreator::type_name() const {
	return "LegacyResidueRetentionRequirement";
}
//////End Creator methods/////

LegacyResidueRetentionRequirement::LegacyResidueRetentionRequirement():
	model_id_(0)
{}

LegacyResidueRetentionRequirement::LegacyResidueRetentionRequirement(
	int model_id
):
	model_id_(model_id)
{}

void
LegacyResidueRetentionRequirement::model_id(
	int model_id
){
	model_id_ = model_id;
}

void
LegacyResidueRetentionRequirement::required_resnums(
	std::set<core::Size> required_resnums
){
	required_resnums_ = required_resnums;
}

void
LegacyResidueRetentionRequirement::add_resnum(
	core::Size resnum
){
	required_resnums_.insert(resnum);
}


bool
LegacyResidueRetentionRequirement::satisfies(
	AssemblyCOP assembly
) const {
	return !violates(assembly);
}

bool
LegacyResidueRetentionRequirement::violates(
	AssemblyCOP assembly
) const {
	auto it = required_resnums_.begin();
	auto it_end = required_resnums_.end();
	for ( ; it != it_end; ++it ) {

		bool found = false;
		utility::vector1<SewSegment> const & segments = assembly->segments();
		for ( core::Size i = 1; i <= segments.size(); ++i ) {
			if ( segments[i].model_id_ == model_id_ ) {
				for ( core::Size j = 1; j <= segments[i].residues_.size(); ++j ) {
					if ( segments[i].residues_[j].resnum_ == *it ) {
						found = true;
						break;
					}
				}
			}
		}

		if ( !found ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Failed to find resnum " << *it << " for model " << model_id_ << std::endl;
			}
			return true;
		}
	}
	return false;
}

void
LegacyResidueRetentionRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){

	if ( !tag->hasOption("model_id") ) {
		utility_exit_with_message("LegacyResidueRetenetionRequirement tag requires a 'model_id' option");
	}
	model_id_ = tag->getOption<int>("model_id");

	if ( !tag->hasOption("required_resnums") ) {
		utility_exit_with_message("LegacyResidueRetenetionRequirement tag requires a whitespace delimited 'resnums' option");
	}
	utility::vector1<std::string> required_resnums_strings = utility::split(tag->getOption<std::string>("resquired_resnums"));
	for ( core::Size i=1; i<=required_resnums_strings.size(); ++i ) {
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );
		core::Size resnum;
		ss << required_resnums_strings[i];
		ss >> resnum;
		required_resnums_.insert(resnum);
	}
}


void
LegacyResidueRetentionRequirement::show(
	std::ostream & out
) const {
	out << "/////// LegacyResidueRetentionRequirement - Assembly must have a segment with ID " << model_id_
		<< " and residues ";
	auto it = required_resnums_.begin();
	auto it_end = required_resnums_.end();
	for ( ; it != it_end; ++it ) {
		out << " " << *it;
	}
	out << std::endl;
}

void
LegacyResidueRetentionRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LegacyResidueRetentionRequirement::provide_xml_schema( xsd );
}


std::string
LegacyResidueRetentionRequirement::class_name(){
	return "LegacyResidueRetentionRequirement";
}

void
LegacyResidueRetentionRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	//min_length and max_length
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "model_id", xs_integer, "Model for which residues are being specified" )
		+ XMLSchemaAttribute::required_attribute( "required_resnums", xsct_nnegative_int_wsslist, "Whitespace delimited list of required residue numbers" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & LegacyRequirementFactory::legacy_global_requirements_ct_namer )
		.element_name( class_name() )
		.description( "Specifies residues that must be retained during assembly" )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}

} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace
